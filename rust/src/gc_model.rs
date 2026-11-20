//! The GC model (GC redesign, item 1): one definition of fragment GC, one table of weights with
//! its provenance, and a strict loader.
//!
//! **Not used by any feature yet.** Every GC-weighted feature still reads
//! [`crate::gc_correction::CorrectionFactors`] until spec Step 2 moves the GC tally and the
//! lookups onto the shared fragment record. Nothing here changes an output.
//!
//! ## The unit
//!
//! A cell is an **exact fragment length** and an **integer GC count** over the fragment's insert
//! with `trim_k` bp dropped from each end ([`gc_key`]). The GC count is never rounded to a
//! percent: at 40 bp a percent is not a whole base, and rounding is what left GC holes at short
//! lengths in the old table. A length `L` has the cells `gc = 0..=L - 2k`.
//!
//! ## Determinism
//!
//! Every table is a dense array indexed by (length, gc) ([`GcGrid`]), and every sum runs in index
//! order. Nothing iterates a `HashMap`, and nothing reduces across threads, so a rerun on the same
//! platform gives a bit-identical table. (Across platforms `f64::exp` comes from the platform
//! libm, so bit-identity there is not claimed.)

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use anyhow::{anyhow, Context, Result};
use arrow::array::{Array, ArrayRef, BooleanArray, Float64Array, UInt16Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use log::info;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use serde::{Deserialize, Serialize};

/// The version of the GC *definition*: [`gc_key`], the trimmed window and the N rule. A model or
/// expected table written under another version is refused on load, never reinterpreted.
pub const GC_DEFINITION_VERSION: u32 = 1;

/// The shortest trimmed window that has a GC key. A fragment shorter than `2k + 20` bp has none.
pub const MIN_TRIMMED_WINDOW: u16 = 20;

/// The fragment-length range of the model, in bp (inclusive). At `trim_k = k` the shortest length
/// with a key is `max(LENGTH_MIN, 2k + MIN_TRIMMED_WINDOW)`; see [`GcGrid::lengths_for`].
pub const LENGTH_MIN: u16 = 20;
pub const LENGTH_MAX: u16 = 600;

/// The trim used until spec Step 3 chooses one from {0, 2, 5, 10} by held-out residual.
pub const DEFAULT_TRIM_K: u16 = 10;

// ---------------------------------------------------------------------------------------------
// The GC key
// ---------------------------------------------------------------------------------------------

/// One cell: an exact fragment length and the GC count over its trimmed window.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GcKey {
    pub length: u16,
    pub gc: u16,
}

/// The GC key of a fragment insert (its reference sequence, 5' to 3', `insert.len()` bp).
///
/// The window is the insert without its first and last `trim_k` bases. A, C, G and T count in
/// either case; **any other base in the window** (N, or an IUPAC code) makes the key `None`, as
/// does an insert shorter than `2 * trim_k + MIN_TRIMMED_WINDOW`. The trimmed-off ends are not
/// read at all. This is the one definition: the observed tally, the expected counts and the lookup
/// all use it, so they cannot disagree about what a cell is.
pub fn gc_key(insert: &[u8], trim_k: u16) -> Option<GcKey> {
    let k = trim_k as usize;
    let length = u16::try_from(insert.len()).ok()?;
    if insert.len() < 2 * k + MIN_TRIMMED_WINDOW as usize {
        return None;
    }
    let mut gc = 0u16;
    for &b in &insert[k..insert.len() - k] {
        match base_class(b) {
            (true, is_gc) => gc += is_gc as u16,
            (false, _) => return None,
        }
    }
    Some(GcKey { length, gc })
}

/// `(counts, is_gc)` for one base: whether it counts toward the window (A/C/G/T, either case),
/// and whether it is G or C. The expected-count builder's prefix sums use this too, so the two
/// cannot disagree about a base.
#[inline]
pub(crate) fn base_class(b: u8) -> (bool, bool) {
    match b {
        b'G' | b'g' | b'C' | b'c' => (true, true),
        b'A' | b'a' | b'T' | b't' => (true, false),
        _ => (false, false),
    }
}

// ---------------------------------------------------------------------------------------------
// The dense grid every table uses
// ---------------------------------------------------------------------------------------------

/// A dense table over every cell of a length range at one trim: length `L` holds `gc = 0..=L-2k`.
#[derive(Debug, Clone, PartialEq)]
pub struct GcGrid<T> {
    min_len: u16,
    max_len: u16,
    trim_k: u16,
    offsets: Vec<usize>,
    cells: Vec<T>,
}

impl<T: Clone> GcGrid<T> {
    /// A grid over `min_len..=max_len` at `trim_k`, every cell set to `fill`. Refuses a range whose
    /// shortest length has no key (`min_len < 2k + MIN_TRIMMED_WINDOW`): a length that can never
    /// hold a fragment would read as one with no support.
    pub fn new(min_len: u16, max_len: u16, trim_k: u16, fill: T) -> Result<Self> {
        let shortest = 2 * trim_k + MIN_TRIMMED_WINDOW;
        if min_len < shortest {
            return Err(anyhow!(
                "GC grid: lengths start at {min_len} bp, but at trim {trim_k} no fragment shorter \
                 than {shortest} bp has a GC key"
            ));
        }
        if min_len > max_len {
            return Err(anyhow!("GC grid: empty length range {min_len}..={max_len}"));
        }
        let mut offsets = Vec::with_capacity((max_len - min_len) as usize + 2);
        let mut n = 0usize;
        for len in min_len..=max_len {
            offsets.push(n);
            n += (len - 2 * trim_k) as usize + 1;
        }
        offsets.push(n);
        Ok(Self { min_len, max_len, trim_k, offsets, cells: vec![fill; n] })
    }
}

impl<T> GcGrid<T> {
    /// The default length range at a trim: `max(LENGTH_MIN, 2k + 20)..=LENGTH_MAX`.
    pub fn lengths_for(trim_k: u16) -> (u16, u16) {
        (LENGTH_MIN.max(2 * trim_k + MIN_TRIMMED_WINDOW), LENGTH_MAX)
    }

    pub fn min_len(&self) -> u16 {
        self.min_len
    }
    pub fn max_len(&self) -> u16 {
        self.max_len
    }
    pub fn trim_k(&self) -> u16 {
        self.trim_k
    }
    pub fn len(&self) -> usize {
        self.cells.len()
    }
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// The largest GC count a length can have: its trimmed window, `L - 2k`.
    pub fn gc_max(&self, length: u16) -> u16 {
        length - 2 * self.trim_k
    }

    /// The flat index of a cell, or `None` outside the grid.
    pub fn index(&self, length: u16, gc: u16) -> Option<usize> {
        if length < self.min_len || length > self.max_len || gc > self.gc_max(length) {
            return None;
        }
        Some(self.offsets[(length - self.min_len) as usize] + gc as usize)
    }

    pub fn get(&self, length: u16, gc: u16) -> Option<&T> {
        self.index(length, gc).map(|i| &self.cells[i])
    }

    /// One length's cells, `gc = 0..=L-2k`, in order.
    pub fn row(&self, length: u16) -> &[T] {
        let i = (length - self.min_len) as usize;
        &self.cells[self.offsets[i]..self.offsets[i + 1]]
    }

    pub fn row_mut(&mut self, length: u16) -> &mut [T] {
        let i = (length - self.min_len) as usize;
        let (a, b) = (self.offsets[i], self.offsets[i + 1]);
        &mut self.cells[a..b]
    }

    pub fn cells(&self) -> &[T] {
        &self.cells
    }

    pub fn cells_mut(&mut self) -> &mut [T] {
        &mut self.cells
    }

    /// Every (length, gc) in index order.
    pub fn keys(&self) -> impl Iterator<Item = GcKey> + '_ {
        (self.min_len..=self.max_len)
            .flat_map(move |length| (0..=self.gc_max(length)).map(move |gc| GcKey { length, gc }))
    }

    /// Same length range and trim, so cells line up index for index.
    pub fn same_shape<U>(&self, other: &GcGrid<U>) -> bool {
        self.min_len == other.min_len && self.max_len == other.max_len && self.trim_k == other.trim_k
    }

    pub(crate) fn shape_string(&self) -> String {
        format!("lengths={}..={} trim_k={}", self.min_len, self.max_len, self.trim_k)
    }
}

// ---------------------------------------------------------------------------------------------
// Support, fit configuration, QC and provenance
// ---------------------------------------------------------------------------------------------

/// The cells a model is fitted on. Outside it a length takes its nearest in-support GC value.
#[derive(Debug, Clone, PartialEq)]
pub struct CellSet {
    cells: GcGrid<bool>,
    share: f64,
}

impl CellSet {
    pub(crate) fn from_grid(cells: GcGrid<bool>, share: f64) -> Self {
        Self { cells, share }
    }

    pub fn contains(&self, length: u16, gc: u16) -> bool {
        self.cells.get(length, gc).copied().unwrap_or(false)
    }

    pub fn grid(&self) -> &GcGrid<bool> {
        &self.cells
    }

    /// The share of expected windows each length's support holds.
    pub fn share(&self) -> f64 {
        self.share
    }

    /// SHA-256 over the shape, the share and every cell's membership in index order.
    pub fn digest(&self) -> String {
        let mut bytes = format!("{} share={:?}\n", self.cells.shape_string(), self.share).into_bytes();
        bytes.extend(self.cells.cells().iter().map(|&b| if b { b'1' } else { b'0' }));
        sha256_hex(&bytes)
    }
}

/// How a model is fitted. Every field is in the compatibility hash ([`Provenance::digest`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FitConfig {
    /// Gaussian kernel SD along fragment length, in bp.
    pub bandwidth_length_bp: f64,
    /// Gaussian kernel SD along GC fraction (GC count over trimmed window).
    pub bandwidth_gc: f64,
    /// Raw weights are capped to `[1/cap, cap]` before normalisation.
    pub cap: f64,
    /// Fewer fragments in the tally than this, and the fit refuses.
    pub min_fragments: u64,
}

/// What the fit did that a reader of the weights cannot see. Reported, never hidden
/// (AGENTS invariant #3): a weight taken from the support's edge or held at the cap is not a
/// measurement.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FitQc {
    /// Fragments in the tally the model was fitted on.
    pub fragments: u64,
    /// Share of the fitted fragments' total weight in cells outside the support.
    pub clamped_weight_share: f64,
    /// Share of the fitted fragments' total weight in cells whose raw weight hit the cap.
    pub capped_weight_share: f64,
    /// Estimated relative smoothing bias, `(h^2 / 2) |s''| / s` summed over both axes:
    /// fragment-weighted mean and maximum over the support.
    pub smoothing_bias_mean: f64,
    pub smoothing_bias_max: f64,
    /// Lengths with no fragment in the tally: their weights are smoothed from neighbours but not
    /// normalised (there is nothing to normalise against).
    pub lengths_unobserved: u32,
    /// Lengths no fragment within the kernel's reach informs: they have no weight (`None`).
    pub lengths_unfitted: u32,
}

/// Where a model came from. [`Self::digest`] is the compatibility hash: two models, or a model and
/// a PON, are comparable only if their digests are equal.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    /// The reference assembly's label (for example `hg19`).
    pub genome: String,
    /// SHA-256 of the territory: every interval's FASTA contig name, contig length, start, end.
    pub territory_sha256: String,
    /// SHA-256 of the support ([`CellSet::digest`]).
    pub support_sha256: String,
    pub trim_k: u16,
    pub gc_definition_version: u32,
    pub fit_config: FitConfig,
    /// Recorded, **not** in the digest: every patch release would otherwise invalidate every PON.
    pub krewlyzer_version: String,
}

impl Provenance {
    /// The compatibility hash: genome, territory, support, trim, GC definition version and fit
    /// configuration. Not the krewlyzer version.
    pub fn digest(&self) -> String {
        let c = &self.fit_config;
        let canonical = format!(
            "genome={}\nterritory_sha256={}\nsupport_sha256={}\ntrim_k={}\ngc_definition_version={}\n\
             bandwidth_length_bp={:?}\nbandwidth_gc={:?}\ncap={:?}\nmin_fragments={}\n",
            self.genome, self.territory_sha256, self.support_sha256, self.trim_k,
            self.gc_definition_version, c.bandwidth_length_bp, c.bandwidth_gc, c.cap, c.min_fragments
        );
        sha256_hex(canonical.as_bytes())
    }
}

// ---------------------------------------------------------------------------------------------
// The model
// ---------------------------------------------------------------------------------------------

/// A fitted GC model: one weight per cell, with the support it was fitted on, its provenance and
/// what the fit could not do.
#[derive(Debug, Clone, PartialEq)]
pub struct GcModel {
    weights: GcGrid<Option<f64>>,
    capped: GcGrid<bool>,
    support: CellSet,
    provenance: Provenance,
    qc: FitQc,
}

const META_PROVENANCE: &str = "krewlyzer.gc_model.provenance";
const META_QC: &str = "krewlyzer.gc_model.qc";
const META_SHAPE: &str = "krewlyzer.gc_model.shape";

#[derive(Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
struct Shape {
    min_len: u16,
    max_len: u16,
    trim_k: u16,
}

impl GcModel {
    pub(crate) fn from_parts(
        weights: GcGrid<Option<f64>>,
        capped: GcGrid<bool>,
        support: CellSet,
        provenance: Provenance,
        qc: FitQc,
    ) -> Result<Self> {
        if !weights.same_shape(support.grid()) || !weights.same_shape(&capped) {
            return Err(anyhow!("GC model: weights, cap flags and support have different shapes"));
        }
        if weights.trim_k() != provenance.trim_k {
            return Err(anyhow!(
                "GC model: grid trim {} differs from provenance trim {}",
                weights.trim_k(), provenance.trim_k
            ));
        }
        if support.digest() != provenance.support_sha256 {
            return Err(anyhow!("GC model: the support does not match its digest in the provenance"));
        }
        Ok(Self { weights, capped, support, provenance, qc })
    }

    /// The weight of a fragment of `length` bp whose key has GC count `gc`.
    ///
    /// `None` when the GC is unresolved (`gc` is `None`), when the length is outside the model or
    /// had nothing to fit from ([`FitQc::lengths_unfitted`]), or when `gc` exceeds the length's
    /// trimmed window. What a consumer does with `None` is Step 3's policy (fragment-model spec
    /// §2d), not this lookup's: it never answers 1.0.
    pub fn weight(&self, length: u16, gc: Option<u16>) -> Option<f64> {
        let gc = gc?;
        self.weights.get(length, gc).copied().flatten()
    }

    pub fn weights(&self) -> &GcGrid<Option<f64>> {
        &self.weights
    }
    pub fn capped(&self) -> &GcGrid<bool> {
        &self.capped
    }
    pub fn support(&self) -> &CellSet {
        &self.support
    }
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }
    pub fn qc(&self) -> &FitQc {
        &self.qc
    }

    /// Write the model as Parquet: one row per cell (`length`, `gc`, `weight` as float64 or null,
    /// `in_support`, `capped`), with the provenance, QC and shape in the file's key-value metadata.
    pub fn write_parquet<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let keys: Vec<GcKey> = self.weights.keys().collect();
        let schema = Arc::new(Schema::new(vec![
            Field::new("length", DataType::UInt16, false),
            Field::new("gc", DataType::UInt16, false),
            Field::new("weight", DataType::Float64, true),
            Field::new("in_support", DataType::Boolean, false),
            Field::new("capped", DataType::Boolean, false),
        ]));
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.length))),
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.gc))),
            Arc::new(Float64Array::from(self.weights.cells().to_vec())),
            Arc::new(BooleanArray::from(self.support.grid().cells().to_vec())),
            Arc::new(BooleanArray::from(self.capped.cells().to_vec())),
        ];
        let shape = Shape {
            min_len: self.weights.min_len(),
            max_len: self.weights.max_len(),
            trim_k: self.weights.trim_k(),
        };
        let meta = vec![
            (META_PROVENANCE, serde_json::to_string(&self.provenance)?),
            (META_QC, serde_json::to_string(&self.qc)?),
            (META_SHAPE, serde_json::to_string(&shape)?),
            (META_SUPPORT_SHARE, serde_json::to_string(&self.support.share())?),
        ];
        write_parquet_with_metadata(path, schema, arrays, meta)?;
        info!("gc_model: wrote {} cells -> {:?}", keys.len(), path);
        Ok(())
    }

    /// Load a model written by [`Self::write_parquet`]. **Refuses**, rather than reinterprets, a
    /// table with no provenance, another GC definition version, a shape that is not the one its
    /// metadata states, a missing, repeated or foreign cell, a non-finite or non-positive weight,
    /// or a support that does not match its digest.
    pub fn load_parquet<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let (meta, batches) = read_parquet_with_metadata(path)?;
        let provenance: Provenance = parse_meta(&meta, META_PROVENANCE, path)?;
        if provenance.gc_definition_version != GC_DEFINITION_VERSION {
            return Err(anyhow!(
                "GC model {:?} was written under GC definition version {}; this krewlyzer reads \
                 version {} only. Refit the model; it is not reinterpreted.",
                path, provenance.gc_definition_version, GC_DEFINITION_VERSION
            ));
        }
        let qc: FitQc = parse_meta(&meta, META_QC, path)?;
        let shape: Shape = parse_meta(&meta, META_SHAPE, path)?;
        if shape.trim_k != provenance.trim_k {
            return Err(anyhow!("GC model {:?}: shape trim {} differs from provenance trim {}",
                path, shape.trim_k, provenance.trim_k));
        }
        let mut weights = GcGrid::new(shape.min_len, shape.max_len, shape.trim_k, None)
            .with_context(|| format!("GC model {:?}", path))?;
        let mut support = GcGrid::new(shape.min_len, shape.max_len, shape.trim_k, false)?;
        let mut capped = GcGrid::new(shape.min_len, shape.max_len, shape.trim_k, false)?;
        let mut seen = vec![false; weights.len()];
        for batch in &batches {
            let length = column::<UInt16Array>(batch, "length", path)?;
            let gc = column::<UInt16Array>(batch, "gc", path)?;
            let weight = column::<Float64Array>(batch, "weight", path)?;
            let in_support = column::<BooleanArray>(batch, "in_support", path)?;
            let cap = column::<BooleanArray>(batch, "capped", path)?;
            for r in 0..batch.num_rows() {
                let i = claim_cell(&weights, &mut seen, length.value(r), gc.value(r), path)?;
                if in_support.is_null(r) || cap.is_null(r) {
                    return Err(anyhow!("GC model {:?}: null flag at length {} gc {}",
                        path, length.value(r), gc.value(r)));
                }
                if !weight.is_null(r) {
                    let w = weight.value(r);
                    if !(w.is_finite() && w > 0.0) {
                        return Err(anyhow!("GC model {:?}: weight {} at length {} gc {} is not a \
                            finite positive number", path, w, length.value(r), gc.value(r)));
                    }
                    weights.cells_mut()[i] = Some(w);
                }
                support.cells_mut()[i] = in_support.value(r);
                capped.cells_mut()[i] = cap.value(r);
            }
        }
        require_every_cell(&seen, path)?;
        let support = CellSet::from_grid(support, read_share(&meta, path)?);
        Self::from_parts(weights, capped, support, provenance, qc)
            .with_context(|| format!("GC model {:?}", path))
    }
}

/// The support share is stored in the provenance's support digest only, so the loader needs it
/// separately to recompute that digest.
const META_SUPPORT_SHARE: &str = "krewlyzer.gc_model.support_share";

fn read_share(meta: &HashMap<String, String>, path: &Path) -> Result<f64> {
    parse_meta(meta, META_SUPPORT_SHARE, path)
}

// ---------------------------------------------------------------------------------------------
// Parquet helpers shared by the model and the expected table
// ---------------------------------------------------------------------------------------------

pub(crate) fn write_parquet_with_metadata(
    path: &Path,
    schema: Arc<Schema>,
    arrays: Vec<ArrayRef>,
    meta: Vec<(&str, String)>,
) -> Result<()> {
    let batch = RecordBatch::try_new(schema.clone(), arrays).context("building GC table batch")?;
    let kv = meta.into_iter().map(|(k, v)| KeyValue::new(k.to_string(), Some(v))).collect();
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_key_value_metadata(Some(kv))
        .build();
    let file = File::create(path).with_context(|| format!("creating {:?}", path))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

pub(crate) fn read_parquet_with_metadata(path: &Path) -> Result<(HashMap<String, String>, Vec<RecordBatch>)> {
    let file = File::open(path).with_context(|| format!("opening GC table {:?}", path))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .with_context(|| format!("GC table {:?}: not a readable Parquet file", path))?;
    let mut meta = HashMap::new();
    if let Some(kvs) = builder.metadata().file_metadata().key_value_metadata() {
        for kv in kvs {
            if let Some(v) = &kv.value {
                meta.insert(kv.key.clone(), v.clone());
            }
        }
    }
    let batches = builder.build()?.collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("GC table {:?}: reading rows", path))?;
    Ok((meta, batches))
}

pub(crate) fn parse_meta<T: for<'de> Deserialize<'de>>(
    meta: &HashMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<T> {
    let raw = meta.get(key).ok_or_else(|| anyhow!(
        "GC table {:?} has no `{}` metadata: it carries no provenance and is refused", path, key))?;
    serde_json::from_str(raw).with_context(|| format!("GC table {:?}: malformed `{}` metadata", path, key))
}

/// A column by name with exactly the type the writer uses; anything else is an error, never a
/// default (architecture.md, "three traps").
pub(crate) fn column<'a, A: Array + 'static>(batch: &'a RecordBatch, name: &str, path: &Path) -> Result<&'a A> {
    let col = batch.column_by_name(name)
        .ok_or_else(|| anyhow!("GC table {:?} has no `{}` column", path, name))?;
    col.as_any().downcast_ref::<A>().ok_or_else(|| anyhow!(
        "GC table {:?}: column `{}` is {:?}, not the type this krewlyzer writes", path, name, col.data_type()))
}

/// Mark one row's cell as seen, refusing a cell outside the shape or a repeated one.
pub(crate) fn claim_cell<T>(grid: &GcGrid<T>, seen: &mut [bool], length: u16, gc: u16, path: &Path) -> Result<usize> {
    let i = grid.index(length, gc).ok_or_else(|| anyhow!(
        "GC table {:?}: cell length {} gc {} is outside its own shape ({})", path, length, gc, grid.shape_string()))?;
    if std::mem::replace(&mut seen[i], true) {
        return Err(anyhow!("GC table {:?}: cell length {} gc {} appears twice", path, length, gc));
    }
    Ok(i)
}

pub(crate) fn require_every_cell(seen: &[bool], path: &Path) -> Result<()> {
    let missing = seen.iter().filter(|s| !**s).count();
    if missing > 0 {
        return Err(anyhow!("GC table {:?}: {} cells of its shape are missing", path, missing));
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// SHA-256 (FIPS 180-4), in-house: the plan allows no new crate dependency, and the digests here
// are identities, not secrets.
// ---------------------------------------------------------------------------------------------

pub(crate) fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    ];
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&((data.len() as u64).wrapping_mul(8)).to_be_bytes());
    for chunk in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for (i, word) in chunk.chunks_exact(4).enumerate() {
            w[i] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
        }
        for (hi, vi) in h.iter_mut().zip(v) {
            *hi = hi.wrapping_add(vi);
        }
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    #[test]
    fn sha256_matches_the_fips_vectors() {
        assert_eq!(sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
        assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        assert_eq!(
            sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    // -- gc_key ---------------------------------------------------------------------------------

    #[test]
    fn gc_key_counts_gc_over_the_trimmed_window_in_either_case() {
        // 24 bp, k = 2: the window is bases 2..22, 20 bp. The ends (GGs) are not read.
        let insert = b"GGacgtACGTacgtACGTacgtGG";
        assert_eq!(gc_key(insert, 2), Some(GcKey { length: 24, gc: 10 }));
        // The same window with AT ends: the same key. A different window: a different key.
        assert_eq!(gc_key(b"ATacgtACGTacgtACGTacgtTA", 2), Some(GcKey { length: 24, gc: 10 }));
        assert_eq!(gc_key(b"GGaaaaACGTacgtACGTacgtGG", 2), Some(GcKey { length: 24, gc: 8 }));
    }

    #[test]
    fn any_base_other_than_acgt_in_the_window_has_no_key() {
        for &bad in b"NnRY-" {
            let mut insert = b"ACGTACGTACGTACGTACGTACGT".to_vec();
            insert[12] = bad;
            assert_eq!(gc_key(&insert, 2), None, "base {:?}", bad as char);
        }
        // An N in the trimmed-off end is not read.
        assert!(gc_key(b"NNACGTACGTACGTACGTACGTNN", 2).is_some());
    }

    #[test]
    fn a_fragment_shorter_than_2k_plus_20_has_no_key() {
        let s = vec![b'A'; 39];
        assert_eq!(gc_key(&s, 10), None);
        assert_eq!(gc_key(&[b'A'; 40], 10), Some(GcKey { length: 40, gc: 0 }));
        assert_eq!(gc_key(&[b'A'; 19], 0), None);
    }

    // -- the model and its file ------------------------------------------------------------------

    pub(crate) fn config() -> FitConfig {
        FitConfig { bandwidth_length_bp: 5.0, bandwidth_gc: 0.02, cap: 10.0, min_fragments: 100_000 }
    }

    fn provenance(territory: &str, support: &CellSet) -> Provenance {
        Provenance {
            genome: "synthetic".into(),
            territory_sha256: sha256_hex(territory.as_bytes()),
            support_sha256: support.digest(),
            trim_k: 2,
            gc_definition_version: GC_DEFINITION_VERSION,
            fit_config: config(),
            krewlyzer_version: "0.0.0-test".into(),
        }
    }

    /// A small model with a non-trivial weight surface, a hole in the support and an unfitted
    /// length, built by hand: the file format is under test here, not the fit.
    fn small_model(territory: &str) -> GcModel {
        let mut weights = GcGrid::new(24, 30, 2, None).unwrap();
        let mut support = GcGrid::new(24, 30, 2, false).unwrap();
        let mut capped = GcGrid::new(24, 30, 2, false).unwrap();
        let keys: Vec<GcKey> = weights.keys().collect();
        for (i, k) in keys.iter().enumerate() {
            if k.length != 27 {
                weights.cells_mut()[i] = Some(1.0 / 3.0 + (k.gc as f64).sqrt() * 0.1 + k.length as f64 * 1e-3);
            }
            support.cells_mut()[i] = (5..=15).contains(&k.gc) && k.gc != 9;
            capped.cells_mut()[i] = k.gc == 0;
        }
        let support = CellSet::from_grid(support, 0.995);
        let qc = FitQc { fragments: 123_456, clamped_weight_share: 0.0123, capped_weight_share: 0.004,
                         smoothing_bias_mean: 0.001, smoothing_bias_max: 0.01,
                         lengths_unobserved: 0, lengths_unfitted: 1 };
        let p = provenance(territory, &support);
        GcModel::from_parts(weights, capped, support, p, qc).unwrap()
    }

    #[test]
    fn write_and_load_round_trip_bit_identically() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("m.gc_model.parquet");
        let m = small_model("1\t1000\t0\t1000\n");
        m.write_parquet(&p).unwrap();
        let back = GcModel::load_parquet(&p).unwrap();
        let bits = |g: &GcGrid<Option<f64>>| g.cells().iter().map(|w| w.map(f64::to_bits)).collect::<Vec<_>>();
        assert_eq!(bits(back.weights()), bits(m.weights()));
        assert_eq!(back, m);
        assert_eq!(back.weight(27, Some(10)), None, "an unfitted length stays unfitted");
    }

    #[test]
    fn a_table_without_provenance_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("bare.parquet");
        let m = small_model("t");
        // Same rows, no metadata.
        let keys: Vec<GcKey> = m.weights().keys().collect();
        let schema = Arc::new(Schema::new(vec![
            Field::new("length", DataType::UInt16, false),
            Field::new("gc", DataType::UInt16, false),
            Field::new("weight", DataType::Float64, true),
            Field::new("in_support", DataType::Boolean, false),
            Field::new("capped", DataType::Boolean, false),
        ]));
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.length))),
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.gc))),
            Arc::new(Float64Array::from(m.weights().cells().to_vec())),
            Arc::new(BooleanArray::from(m.support().grid().cells().to_vec())),
            Arc::new(BooleanArray::from(m.capped().cells().to_vec())),
        ];
        crate::output_utils::write_parquet_batch(&p, schema, arrays).unwrap();
        let err = format!("{:#}", GcModel::load_parquet(&p).unwrap_err());
        assert!(err.contains("provenance"), "{err}");
    }

    #[test]
    fn a_model_from_another_gc_definition_version_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("v2.parquet");
        let mut m = small_model("t");
        m.provenance.gc_definition_version = GC_DEFINITION_VERSION + 1;
        m.write_parquet(&p).unwrap();
        let err = format!("{:#}", GcModel::load_parquet(&p).unwrap_err());
        assert!(err.contains("definition version"), "{err}");
    }

    #[test]
    fn a_tampered_support_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("t.parquet");
        let mut m = small_model("t");
        let i = m.support.cells.index(25, 9).unwrap();
        m.support.cells.cells_mut()[i] = true; // the digest in the provenance no longer matches
        m.write_parquet(&p).unwrap();
        let err = format!("{:#}", GcModel::load_parquet(&p).unwrap_err());
        assert!(err.contains("support"), "{err}");
    }

    #[test]
    fn two_territories_give_two_digests_and_the_krewlyzer_version_gives_one() {
        let a = small_model("1\t1000\t0\t1000\n");
        let b = small_model("1\t1000\t0\t999\n");
        assert_ne!(a.provenance().digest(), b.provenance().digest());
        let mut c = a.provenance().clone();
        c.krewlyzer_version = "9.9.9".into();
        assert_eq!(a.provenance().digest(), c.digest(), "krewlyzer_version is recorded, not hashed");
        let mut d = a.provenance().clone();
        d.trim_k = 5;
        assert_ne!(a.provenance().digest(), d.digest());
        let mut e = a.provenance().clone();
        e.fit_config.cap = 20.0;
        assert_ne!(a.provenance().digest(), e.digest());
    }

    #[test]
    fn weight_is_none_for_an_unresolved_gc_and_outside_the_model() {
        let m = small_model("t");
        assert_eq!(m.weight(25, None), None);
        assert_eq!(m.weight(23, Some(5)), None, "below the length range");
        assert_eq!(m.weight(31, Some(5)), None, "above the length range");
        assert_eq!(m.weight(25, Some(22)), None, "above the length's trimmed window");
        let w = m.weight(25, Some(5)).unwrap();
        assert_ne!(Some(w), m.weight(25, Some(12)), "two cells, two weights");
    }
}
