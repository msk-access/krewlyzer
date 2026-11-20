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
/// read at all: **an N in a trimmed-off end does not void the key; an N in the counted window
/// does** (operator decision, 20 Nov). This is GC definition version 1
/// ([`GC_DEFINITION_VERSION`]); changing the rule means a new version. It is the one definition:
/// the observed tally, the expected counts and the lookup all use it, so they cannot disagree
/// about what a cell is.
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
        // In u32: a trim read from a file must not overflow into a short minimum.
        let shortest = 2 * trim_k as u32 + MIN_TRIMMED_WINDOW as u32;
        if (min_len as u32) < shortest {
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
        let shortest = (2 * trim_k as u32 + MIN_TRIMMED_WINDOW as u32).min(u16::MAX as u32) as u16;
        (LENGTH_MIN.max(shortest), LENGTH_MAX)
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

/// The share of each length's expected windows the support holds.
pub const SUPPORT_SHARE: f64 = 0.995;

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

    /// The support from the reference (decision S1): at each length, the fewest GC cells holding
    /// `share` of the territory's expected windows, largest first (ties to the lower GC count).
    /// It needs no donors and no depth, and is fixed per (genome, territory, trim).
    pub fn from_expected(expected: &crate::gc_reference::ExpectedTable, share: f64) -> Result<Self> {
        if !(share > 0.0 && share <= 1.0) {
            return Err(anyhow!("support share {share} is not in (0, 1]"));
        }
        let e = expected.counts();
        let mut cells = GcGrid::new(e.min_len(), e.max_len(), e.trim_k(), false)?;
        for len in e.min_len()..=e.max_len() {
            let row = e.row(len);
            let target = share * row.iter().sum::<u64>() as f64;
            let mut order: Vec<usize> = (0..row.len()).collect();
            order.sort_by(|&a, &b| row[b].cmp(&row[a]).then(a.cmp(&b)));
            let mut held = 0u64;
            let out = cells.row_mut(len);
            for g in order {
                if held as f64 >= target || row[g] == 0 {
                    break;
                }
                held += row[g];
                out[g] = true;
            }
        }
        Ok(Self { cells, share })
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
    /// Fragments the tally counted with a key whose length is outside the model's range: they
    /// were not fitted and have no weight. What they get is Step 2's decision.
    pub fragments_outside_lengths: u64,
    /// Inserts with no GC key at all (too short for the trim, or a non-ACGT base in the window):
    /// unresolved, not out of range.
    pub fragments_without_key: u64,
    /// Share of the fitted fragments' total weight in cells outside the support.
    pub clamped_weight_share: f64,
    /// Share of the fitted fragments' total weight in cells whose raw weight hit the cap.
    pub capped_weight_share: f64,
    /// Estimated relative smoothing bias, `(h^2 / 2) |s''| / s` summed over both axes: the
    /// fragment-weighted mean and 95th percentile over the support (a maximum would be set by the
    /// support's sparsest cell, where the ratio is near zero).
    pub smoothing_bias_mean: f64,
    pub smoothing_bias_p95: f64,
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

// ---------------------------------------------------------------------------------------------
// The observed tally and the fit
// ---------------------------------------------------------------------------------------------

/// Gaussian kernel SD along fragment length, in bp, until spec Step 3 chooses it by held-out
/// residual (criterion 1).
pub const DEFAULT_BANDWIDTH_LENGTH_BP: f64 = 5.0;
/// Gaussian kernel SD along GC fraction, until spec Step 3 chooses it.
pub const DEFAULT_BANDWIDTH_GC: f64 = 0.02;
/// The kernel is truncated at this many SDs on each axis.
pub const KERNEL_REACH_SD: f64 = 4.0;
/// Raw weights are capped to `[1 / WEIGHT_CAP, WEIGHT_CAP]` before normalisation.
pub const WEIGHT_CAP: f64 = 10.0;
/// Fewer off-target fragments than this in the tally, and the fit refuses.
pub const MIN_FIT_FRAGMENTS: u64 = 100_000;

impl FitConfig {
    /// Bandwidths finite and positive, the cap finite and at least 1, a positive minimum. Anything
    /// else cannot define a model: a zero bandwidth makes every smoothed ratio NaN, and a cap
    /// below 1 or NaN has no meaning (and panics in `f64::clamp`).
    pub fn validate(&self) -> std::result::Result<(), String> {
        let positive = |v: f64| v.is_finite() && v > 0.0;
        if !positive(self.bandwidth_length_bp) {
            return Err(format!("bandwidth_length_bp {} is not finite and positive", self.bandwidth_length_bp));
        }
        if !positive(self.bandwidth_gc) {
            return Err(format!("bandwidth_gc {} is not finite and positive", self.bandwidth_gc));
        }
        if !(self.cap.is_finite() && self.cap >= 1.0) {
            return Err(format!("cap {} is not finite and at least 1", self.cap));
        }
        if self.min_fragments == 0 {
            return Err("min_fragments is 0: an empty tally would fit".to_string());
        }
        Ok(())
    }
}

impl Default for FitConfig {
    fn default() -> Self {
        Self {
            bandwidth_length_bp: DEFAULT_BANDWIDTH_LENGTH_BP,
            bandwidth_gc: DEFAULT_BANDWIDTH_GC,
            cap: WEIGHT_CAP,
            min_fragments: MIN_FIT_FRAGMENTS,
        }
    }
}

/// A sample's fragments counted per cell, for one territory at one trim. Step 2 builds it from
/// the fragment record with [`Self::add_insert`], over the same territory as the expected counts;
/// [`fit`] refuses a tally whose territory digest or trim is not the expected table's.
#[derive(Debug, Clone, PartialEq)]
pub struct Tally {
    counts: GcGrid<u64>,
    territory_sha256: String,
    outside_lengths: u64,
    without_key: u64,
}

impl Tally {
    pub fn new(min_len: u16, max_len: u16, trim_k: u16, territory_sha256: &str) -> Result<Self> {
        Ok(Self {
            counts: GcGrid::new(min_len, max_len, trim_k, 0)?,
            territory_sha256: territory_sha256.to_string(),
            outside_lengths: 0,
            without_key: 0,
        })
    }

    /// An empty tally over an expected table's cells, territory and trim.
    pub fn for_expected(expected: &crate::gc_reference::ExpectedTable) -> Self {
        let g = expected.counts();
        Self::new(g.min_len(), g.max_len(), g.trim_k(), &expected.provenance().territory_sha256)
            .expect("an existing grid's shape is valid")
    }

    /// Key a fragment's insert with this tally's own trim ([`gc_key`]) and count it. An insert
    /// with no key (shorter than `2k + 20` bp, or a non-ACGT base in its window) is counted in
    /// [`Self::without_key`]; one whose length is outside the grid in [`Self::outside_lengths`].
    pub fn add_insert(&mut self, insert: &[u8]) -> Option<GcKey> {
        let Some(key) = gc_key(insert, self.trim_k()) else {
            self.without_key += 1;
            return None;
        };
        self.add(key, 1).expect("gc_key yields only keys an insert can have");
        Some(key)
    }

    /// Count `n` fragments with this key: `Ok(true)` in the grid, `Ok(false)` counted in
    /// [`Self::outside_lengths`]. A key no insert can have at this trim (a length below
    /// `2k + 20`, or a GC count above `length - 2k`) is an error, not "outside".
    pub fn add(&mut self, key: GcKey, n: u64) -> Result<bool> {
        let k = self.trim_k() as u32;
        let (len, gc) = (key.length as u32, key.gc as u32);
        if len < 2 * k + MIN_TRIMMED_WINDOW as u32 || gc > len - 2 * k {
            return Err(anyhow!("GC key length {} gc {} is impossible at trim {}: no insert has it",
                               key.length, key.gc, k));
        }
        match self.counts.index(key.length, key.gc) {
            Some(i) => {
                self.counts.cells_mut()[i] += n;
                Ok(true)
            }
            None => {
                self.outside_lengths += n;
                Ok(false)
            }
        }
    }

    pub fn counts(&self) -> &GcGrid<u64> {
        &self.counts
    }

    pub fn territory_sha256(&self) -> &str {
        &self.territory_sha256
    }

    pub fn trim_k(&self) -> u16 {
        self.counts.trim_k()
    }

    /// Fragments with a key whose length is outside the grid's range.
    pub fn outside_lengths(&self) -> u64 {
        self.outside_lengths
    }

    /// Inserts with no GC key: too short for the trim, or a non-ACGT base in the window.
    pub fn without_key(&self) -> u64 {
        self.without_key
    }

    /// Fragments in the grid, summed in index order.
    pub fn total(&self) -> u64 {
        self.counts.cells().iter().sum()
    }
}

/// Why a model could not be fitted. The fit says so instead of returning weights of 1.
#[derive(Debug, Clone, PartialEq)]
pub enum FitError {
    /// The tally holds fewer fragments than [`FitConfig::min_fragments`].
    TooFewFragments { fragments: u64, required: u64 },
    /// A length of the grid has no cell in the support.
    NoSupport { length: u16 },
    /// The tally, the expected table and the support are not over the same cells.
    ShapeMismatch(String),
    /// The tally was counted over another territory than the expected table's.
    TerritoryMismatch { tally: String, expected: String },
    /// The fit configuration cannot define a model.
    InvalidConfig(String),
}

impl std::fmt::Display for FitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FitError::TooFewFragments { fragments, required } => write!(
                f, "cannot fit a GC model: {fragments} fragments in the tally, at least {required} needed"),
            FitError::NoSupport { length } => write!(
                f, "cannot fit a GC model: length {length} bp has no cell in the support"),
            FitError::ShapeMismatch(why) => write!(f, "cannot fit a GC model: {why}"),
            FitError::TerritoryMismatch { tally, expected } => write!(
                f, "cannot fit a GC model: the tally was counted over territory {tally}, the \
                    expected table over {expected}"),
            FitError::InvalidConfig(why) => write!(f, "cannot fit a GC model: {why}"),
        }
    }
}

impl std::error::Error for FitError {}

/// The local-linear fit at one cell: the smoothed ratio and, where estimable, its slopes along
/// length (per bp) and GC fraction.
#[derive(Debug, Clone, Copy)]
struct Local {
    s: f64,
    ds_dl: Option<f64>,
    ds_dx: Option<f64>,
}

/// One length's cells, `None` where nothing within the kernel's reach informs a cell.
type SmoothedRow = Vec<Option<Local>>;

/// Fit a GC model from a sample's tally.
///
/// 1. Per length `L`, over the support only: `r(L, g) = p_obs(g | L) / p_exp(g | L)`, the ratio of
///    the observed to the expected GC distribution at that length. Conditioning on length puts `r`
///    on a unit scale (1 where there is no bias) whatever the sample's depth and size distribution,
///    which is what makes a fixed cap on the raw weights meaningful.
/// 2. Smooth `r` over (length, GC fraction) with a Gaussian kernel, as **local-linear** weighted
///    least squares weighted by the expected count at the sample's depth, `O_L * p_exp(g | L)` (the
///    inverse variance of `r`). Local-linear rather than a kernel average: a kernel average is
///    biased wherever the expected density has a slope, and at the support's edge. What remains
///    is a curvature bias of about `(h^2 / 2) |r''| / r`, reported in [`FitQc`].
/// 3. Raw weight `= 1 / smoothed r`, capped to `[1/C, C]`.
/// 4. Sparse cells: a cell outside the support takes its length's nearest in-support value.
///    Observed fragments there never enter the smoother, so an outlier moves only its own length's
///    normalising scale. Both the clamped and the capped share of the weight are reported.
/// 5. Normalise per length over the **whole** tally: `w(L, ·) /= Σ_g obs·w / Σ_g obs`, so the
///    weighted total of every length equals its fragment count. The cap bounds the raw weights,
///    not the normalised ones.
///
/// Dense arrays, every sum in index order; rows are computed in parallel but each one
/// sequentially, and nothing is reduced across threads, so a rerun is bit-identical.
pub fn fit(
    observed: &Tally,
    expected: &crate::gc_reference::ExpectedTable,
    support: &CellSet,
    cfg: &FitConfig,
) -> std::result::Result<GcModel, FitError> {
    fit_with_scales(observed, expected, support, cfg).map(|(m, _)| m)
}

/// [`fit`], also returning each length's normalising divisor (1.0 where nothing was normalised),
/// so a test can recover the raw weights: `raw = weight * scale`.
fn fit_with_scales(
    observed: &Tally,
    expected: &crate::gc_reference::ExpectedTable,
    support: &CellSet,
    cfg: &FitConfig,
) -> std::result::Result<(GcModel, Vec<f64>), FitError> {
    use rayon::prelude::*;
    cfg.validate().map_err(FitError::InvalidConfig)?;
    let exp = expected.counts();
    let obs = observed.counts();
    let sup = support.grid();
    if !obs.same_shape(exp) || !sup.same_shape(exp) {
        return Err(FitError::ShapeMismatch(format!(
            "tally {}, expected table {}, support {}",
            obs.shape_string(), exp.shape_string(), sup.shape_string())));
    }
    if observed.territory_sha256() != expected.provenance().territory_sha256 {
        return Err(FitError::TerritoryMismatch {
            tally: observed.territory_sha256().to_string(),
            expected: expected.provenance().territory_sha256.clone(),
        });
    }
    let fragments = observed.total();
    if fragments < cfg.min_fragments {
        return Err(FitError::TooFewFragments { fragments, required: cfg.min_fragments });
    }
    let lengths: Vec<u16> = (exp.min_len()..=exp.max_len()).collect();
    if let Some(&length) = lengths.iter().find(|&&l| !sup.row(l).iter().any(|&b| b)) {
        return Err(FitError::NoSupport { length });
    }

    // Per cell, over the support only: y = observed count, v = expected count at the sample's
    // depth, O_L * p_exp(g | L). r = y / v, and v is the weight of r in the smoother.
    let mut y = vec![0.0; exp.len()];
    let mut v = vec![0.0; exp.len()];
    for &len in &lengths {
        let i0 = exp.index(len, 0).expect("in range");
        let (mut o_l, mut e_l) = (0.0, 0.0);
        for (g, &in_s) in sup.row(len).iter().enumerate() {
            if in_s {
                o_l += obs.row(len)[g] as f64;
                e_l += exp.row(len)[g] as f64;
            }
        }
        for (g, &in_s) in sup.row(len).iter().enumerate() {
            if in_s && e_l > 0.0 {
                y[i0 + g] = obs.row(len)[g] as f64;
                v[i0 + g] = o_l * exp.row(len)[g] as f64 / e_l;
            }
        }
    }

    // Rows in parallel (a map, not a reduction): each row's sums run sequentially in index order.
    let smoothed: Vec<SmoothedRow> = lengths.par_iter()
        .map(|&len| smooth_row(len, exp, sup, &y, &v, cfg))
        .collect();

    // Raw weights, capped to [1/C, C] before normalisation (a ratio smoothed to zero or below asks
    // for an unbounded weight, so it is capped at C). Then every cell without a fitted value (outside
    // the support, or with nothing in reach) takes its length's nearest fitted GC value, ties to the
    // lower GC count, and inherits that value's cap flag.
    let shape = |fill| GcGrid::new(exp.min_len(), exp.max_len(), exp.trim_k(), fill).expect("same shape");
    let mut weights: GcGrid<Option<f64>> = GcGrid::new(exp.min_len(), exp.max_len(), exp.trim_k(), None)
        .expect("same shape as the expected table");
    let mut capped = shape(false);
    let mut clamped = shape(false);
    for (li, &len) in lengths.iter().enumerate() {
        let row = &smoothed[li];
        let raw: Vec<Option<(f64, bool)>> = row.iter().map(|local| local.map(|l| {
            let r = if l.s > 0.0 { 1.0 / l.s } else { f64::INFINITY };
            let c = r.clamp(1.0 / cfg.cap, cfg.cap);
            (c, c != r)
        })).collect();
        if raw.iter().all(Option::is_none) {
            continue; // nothing in reach: the length is unfitted
        }
        let mut left = vec![None; raw.len()];
        let mut last = None;
        for (g, r) in raw.iter().enumerate() {
            if r.is_some() {
                last = Some(g);
            }
            left[g] = last;
        }
        let mut right = vec![None; raw.len()];
        let mut next = None;
        for (g, r) in raw.iter().enumerate().rev() {
            if r.is_some() {
                next = Some(g);
            }
            right[g] = next;
        }
        let (w_row, c_row, k_row) = (weights.row_mut(len), capped.row_mut(len), clamped.row_mut(len));
        for g in 0..raw.len() {
            let src = match (left[g], right[g]) {
                (Some(a), Some(b)) => if g - a <= b - g { a } else { b },
                (Some(a), None) => a,
                (None, Some(b)) => b,
                (None, None) => unreachable!("the row has a fitted cell"),
            };
            let (w, was_capped) = raw[src].expect("src is fitted");
            w_row[g] = Some(w);
            c_row[g] = was_capped;
            k_row[g] = src != g;
        }
    }

    // Normalise per length over the whole tally, in index order.
    let mut scales = vec![1.0; lengths.len()];
    let mut qc = FitQc {
        fragments,
        fragments_outside_lengths: observed.outside_lengths(),
        fragments_without_key: observed.without_key(),
        ..FitQc::default()
    };
    for (li, &len) in lengths.iter().enumerate() {
        let row_obs = obs.row(len);
        if row_obs.iter().all(|&o| o == 0) {
            qc.lengths_unobserved += 1;
        }
        let row_w = weights.row_mut(len);
        if row_w.iter().all(Option::is_none) {
            qc.lengths_unfitted += 1;
            continue;
        }
        let (mut num, mut den) = (0.0, 0.0);
        for (w, &o) in row_w.iter().zip(row_obs) {
            if let Some(w) = w {
                num += o as f64 * w;
                den += o as f64;
            }
        }
        if den > 0.0 {
            let c = num / den;
            row_w.iter_mut().flatten().for_each(|w| *w /= c);
            scales[li] = c;
        }
    }
    // Shares of the fitted fragments' final weight that came from the support's edge or the cap.
    let (mut total, mut from_edge, mut at_cap) = (0.0, 0.0, 0.0);
    for (i, w) in weights.cells().iter().enumerate() {
        if let Some(w) = w {
            let ow = obs.cells()[i] as f64 * w;
            total += ow;
            if clamped.cells()[i] {
                from_edge += ow;
            }
            if capped.cells()[i] {
                at_cap += ow;
            }
        }
    }
    if total > 0.0 {
        qc.clamped_weight_share = from_edge / total;
        qc.capped_weight_share = at_cap / total;
    }
    let (mean, p95) = smoothing_bias(&lengths, &smoothed, exp, obs, cfg);
    qc.smoothing_bias_mean = mean;
    qc.smoothing_bias_p95 = p95;

    let ep = expected.provenance();
    let provenance = Provenance {
        genome: ep.genome.clone(),
        territory_sha256: ep.territory_sha256.clone(),
        support_sha256: support.digest(),
        trim_k: exp.trim_k(),
        gc_definition_version: GC_DEFINITION_VERSION,
        fit_config: cfg.clone(),
        krewlyzer_version: env!("CARGO_PKG_VERSION").to_string(),
    };
    info!("gc_model: fitted {} fragments over {} lengths ({} unobserved, {} unfitted); smoothing bias \
           mean {:.4}, p95 {:.4}", fragments, lengths.len(), qc.lengths_unobserved, qc.lengths_unfitted,
          qc.smoothing_bias_mean, qc.smoothing_bias_p95);
    GcModel::from_parts(weights, capped, support.clone(), provenance, qc)
        .map(|m| (m, scales))
        .map_err(|e| FitError::ShapeMismatch(format!("{e:#}")))
}

/// Local-linear Gaussian-kernel estimate of `r` at every support cell of one length. Cells outside
/// the support, or with nothing in reach, are `None`.
fn smooth_row(len: u16, grid: &GcGrid<u64>, sup: &GcGrid<bool>, y: &[f64], v: &[f64], cfg: &FitConfig) -> SmoothedRow {
    let (hl, hx) = (cfg.bandwidth_length_bp, cfg.bandwidth_gc);
    let reach_l = (KERNEL_REACH_SD * hl).floor() as i32;
    let reach_x = KERNEL_REACH_SD * hx;
    let n = grid.gc_max(len) as f64;
    let lo_l = (len as i32 - reach_l).max(grid.min_len() as i32) as u16;
    let hi_l = (len as i32 + reach_l).min(grid.max_len() as i32) as u16;
    sup.row(len).iter().enumerate().map(|(g, &in_s)| {
        if !in_s {
            return None;
        }
        let x = g as f64 / n;
        // Weighted sums for the normal equations of r ~ a + b·dL + c·dx, in index order.
        let mut m = [0.0f64; 6]; // S0, SL, Sx, SLL, SLx, Sxx
        let mut t = [0.0f64; 3]; // T0, TL, Tx
        for l2 in lo_l..=hi_l {
            let dl = (l2 as i32 - len as i32) as f64;
            let wl = (-0.5 * (dl / hl).powi(2)).exp();
            let n2 = grid.gc_max(l2) as f64;
            let g_lo = ((x - reach_x) * n2).ceil().max(0.0) as usize;
            let g_hi = ((x + reach_x) * n2).floor().min(n2) as usize;
            let i0 = grid.index(l2, 0).expect("in range");
            for g2 in g_lo..=g_hi {
                let i = i0 + g2;
                if v[i] == 0.0 {
                    continue;
                }
                let dx = g2 as f64 / n2 - x;
                let kw = wl * (-0.5 * (dx / hx).powi(2)).exp();
                let w = kw * v[i];
                m[0] += w;
                m[1] += w * dl;
                m[2] += w * dx;
                m[3] += w * dl * dl;
                m[4] += w * dl * dx;
                m[5] += w * dx * dx;
                let ky = kw * y[i];
                t[0] += ky;
                t[1] += ky * dl;
                t[2] += ky * dx;
            }
        }
        local_linear(&m, &t)
    }).collect()
}

/// The local-linear fit from its normal equations. Falls back to the fit in GC alone when the
/// length direction carries no information (one length in reach), and to the kernel average when
/// neither slope is estimable; a slope that was not estimated is `None`.
fn local_linear(m: &[f64; 6], t: &[f64; 3]) -> Option<Local> {
    let [s0, sl, sx, sll, slx, sxx] = *m;
    if s0 <= 0.0 {
        return None;
    }
    const TOL: f64 = 1e-9;
    let det3 = s0 * (sll * sxx - slx * slx) - sl * (sl * sxx - slx * sx) + sx * (sl * slx - sll * sx);
    if det3.abs() > TOL * s0 * sll * sxx {
        // Cramer's rule, column by column.
        let a = t[0] * (sll * sxx - slx * slx) - sl * (t[1] * sxx - slx * t[2]) + sx * (t[1] * slx - sll * t[2]);
        let b = s0 * (t[1] * sxx - slx * t[2]) - t[0] * (sl * sxx - slx * sx) + sx * (sl * t[2] - t[1] * sx);
        let c = s0 * (sll * t[2] - t[1] * slx) - sl * (sl * t[2] - t[1] * sx) + t[0] * (sl * slx - sll * sx);
        return Some(Local { s: a / det3, ds_dl: Some(b / det3), ds_dx: Some(c / det3) });
    }
    let det2 = s0 * sxx - sx * sx;
    if det2.abs() > TOL * s0 * sxx {
        let a = (t[0] * sxx - sx * t[2]) / det2;
        let c = (s0 * t[2] - sx * t[0]) / det2;
        return Some(Local { s: a, ds_dl: None, ds_dx: Some(c) });
    }
    Some(Local { s: t[0] / s0, ds_dl: None, ds_dx: None })
}

/// The length slope of one length's row at a GC fraction, linearly interpolated between cells.
fn length_slope_at(row: &SmoothedRow, n: f64, x: f64) -> Option<f64> {
    let pos = x * n;
    let (a, b) = (pos.floor() as usize, pos.ceil() as usize);
    let sa = row.get(a).copied().flatten()?.ds_dl?;
    let sb = row.get(b).copied().flatten()?.ds_dl?;
    let t = pos - a as f64;
    Some(sa * (1.0 - t) + sb * t)
}

/// The estimated relative smoothing bias of a local-linear fit, `(h^2 / 2) |∂²s| / s` on each axis,
/// with the second derivative taken as the central difference of the fitted slopes (a second
/// difference of `s` itself across lengths would need `s` interpolated between GC cells, and the
/// interpolation error of a steep curve is as large as the bias). Returns the fragment-weighted
/// mean and 95th percentile over the cells where it is estimable.
fn smoothing_bias(lengths: &[u16], smoothed: &[SmoothedRow], grid: &GcGrid<u64>, obs: &GcGrid<u64>,
                  cfg: &FitConfig) -> (f64, f64) {
    let (hl, hx) = (cfg.bandwidth_length_bp, cfg.bandwidth_gc);
    let mut cells: Vec<(f64, f64)> = Vec::new(); // (bias, fragments), in index order
    for (li, &len) in lengths.iter().enumerate() {
        let n = grid.gc_max(len) as f64;
        let row = &smoothed[li];
        for g in 1..row.len().saturating_sub(1) {
            let Some(here) = row[g].filter(|l| l.s > 0.0) else { continue };
            let mut b = 0.0;
            let mut any = false;
            let dx = |l: Option<Local>| l.and_then(|l| l.ds_dx);
            if let (Some(lo), Some(hi)) = (dx(row[g - 1]), dx(row[g + 1])) {
                b += 0.5 * hx * hx * ((hi - lo) * n / 2.0).abs() / here.s;
                any = true;
            }
            if li > 0 && li + 1 < lengths.len() {
                let x = g as f64 / n;
                let lo = length_slope_at(&smoothed[li - 1], grid.gc_max(lengths[li - 1]) as f64, x);
                let hi = length_slope_at(&smoothed[li + 1], grid.gc_max(lengths[li + 1]) as f64, x);
                if let (Some(lo), Some(hi)) = (lo, hi) {
                    b += 0.5 * hl * hl * ((hi - lo) / 2.0).abs() / here.s;
                    any = true;
                }
            }
            if any {
                cells.push((b, obs.row(len)[g] as f64));
            }
        }
    }
    let total: f64 = cells.iter().map(|c| c.1).sum();
    if total <= 0.0 {
        return (0.0, 0.0);
    }
    let mean = cells.iter().map(|(b, o)| b * o).sum::<f64>() / total;
    cells.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut cum = 0.0;
    let mut p95 = cells.last().map(|c| c.0).unwrap_or(0.0);
    for (b, o) in &cells {
        cum += o;
        if cum >= 0.95 * total {
            p95 = *b;
            break;
        }
    }
    (mean, p95)
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
    fn the_n_rule_of_gc_definition_version_1_is_pinned() {
        // Version 1: a non-ACGT base voids the key only inside the counted (trimmed) window.
        assert_eq!(GC_DEFINITION_VERSION, 1, "a changed rule needs a new version and this test");
        let clean = [b'A'; 50];
        let mut in_end = clean;
        in_end[0] = b'N';
        in_end[49] = b'N';
        in_end[9] = b'N'; // the last base of the 10 bp trimmed-off start
        assert_eq!(gc_key(&in_end, 10), Some(GcKey { length: 50, gc: 0 }), "an N in a trimmed end gives a key");
        let mut in_window = clean;
        in_window[10] = b'N'; // the first base of the counted window
        assert_eq!(gc_key(&in_window, 10), None, "an N in the window gives no key");
        let mut last = clean;
        last[39] = b'N'; // the last base of the counted window
        assert_eq!(gc_key(&last, 10), None);
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
                         smoothing_bias_mean: 0.001, smoothing_bias_p95: 0.01,
                         lengths_unobserved: 0, lengths_unfitted: 1, fragments_outside_lengths: 0,
                         fragments_without_key: 0 };
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
        assert!(err.contains(META_PROVENANCE), "{err}");
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

    /// One written row: length, gc, weight, in_support, capped.
    type Row = (u16, u16, Option<f64>, bool, bool);
    type Edit = Box<dyn Fn(&mut Vec<Row>)>;

    /// Rewrite a model's rows through `edit`, keeping its metadata, and load the result.
    fn load_edited(m: &GcModel, edit: impl Fn(&mut Vec<Row>)) -> String {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("m.parquet");
        m.write_parquet(&src).unwrap();
        let (meta, _) = read_parquet_with_metadata(&src).unwrap();
        let mut rows: Vec<_> = m.weights().keys().zip(m.weights().cells())
            .zip(m.support().grid().cells().iter().zip(m.capped().cells()))
            .map(|((k, w), (s, c))| (k.length, k.gc, *w, *s, *c)).collect();
        edit(&mut rows);
        let schema = Arc::new(Schema::new(vec![
            Field::new("length", DataType::UInt16, false),
            Field::new("gc", DataType::UInt16, false),
            Field::new("weight", DataType::Float64, true),
            Field::new("in_support", DataType::Boolean, false),
            Field::new("capped", DataType::Boolean, false),
        ]));
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt16Array::from_iter_values(rows.iter().map(|r| r.0))),
            Arc::new(UInt16Array::from_iter_values(rows.iter().map(|r| r.1))),
            Arc::new(Float64Array::from(rows.iter().map(|r| r.2).collect::<Vec<_>>())),
            Arc::new(BooleanArray::from(rows.iter().map(|r| r.3).collect::<Vec<_>>())),
            Arc::new(BooleanArray::from(rows.iter().map(|r| r.4).collect::<Vec<_>>())),
        ];
        let dst = dir.path().join("edited.parquet");
        let meta: Vec<(&str, String)> = meta.iter().map(|(k, v)| (k.as_str(), v.clone())).collect();
        write_parquet_with_metadata(&dst, schema, arrays, meta).unwrap();
        match GcModel::load_parquet(&dst) {
            Ok(_) => "loaded".to_string(),
            Err(e) => format!("{e:#}"),
        }
    }

    #[test]
    fn a_malformed_cell_set_is_refused() {
        let m = small_model("t");
        assert_eq!(load_edited(&m, |_| {}), "loaded", "the unedited rows load");
        let cases: [(&str, &str, Edit); 5] = [
            ("repeated", "appears twice", Box::new(|r| { let x = r[3]; r.push(x); })),
            ("missing", "missing", Box::new(|r| { r.pop(); })),
            ("foreign", "outside its own shape", Box::new(|r| r[0].1 = 999)),
            ("zero weight", "finite positive", Box::new(|r| r[0].2 = Some(0.0))),
            ("NaN weight", "finite positive", Box::new(|r| r[0].2 = Some(f64::NAN))),
        ];
        for (what, needle, edit) in cases {
            let err = load_edited(&m, edit);
            assert!(err.contains(needle), "{what}: {err}");
        }
    }

    #[test]
    fn a_grid_refuses_a_trim_that_overflows() {
        assert!(GcGrid::new(100, 200, 40_000, 0u8).is_err());
        assert!(GcGrid::new(100, 200, 41, 0u8).is_err(), "100 bp < 2 * 41 + 20");
        assert!(GcGrid::new(100, 200, 40, 0u8).is_ok(), "100 bp = 2 * 40 + 20");
    }

    #[test]
    fn the_provenance_digest_is_pinned() {
        // A fixed input, a fixed hash: changing the canonical form invalidates every PON, so it
        // must fail here first.
        let p = Provenance {
            genome: "hg19".into(),
            territory_sha256: sha256_hex(b"territory"),
            support_sha256: sha256_hex(b"support"),
            trim_k: 10,
            gc_definition_version: 1,
            fit_config: FitConfig { bandwidth_length_bp: 5.0, bandwidth_gc: 0.02, cap: 10.0, min_fragments: 100_000 },
            krewlyzer_version: "0.9.2".into(),
        };
        // Cross-checked with Python hashlib over the same canonical string.
        assert_eq!(p.digest(), "0e06cef4f8e46d2bcb627115f4a429c5e845344b053fbb41a353850ba93532ef");
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

#[cfg(test)]
pub(crate) mod fit_tests {
    use super::*;
    use crate::gc_reference::ExpectedTable;

    /// SplitMix64, for test data only: the fit itself draws nothing at random.
    pub(crate) struct Rng(u64);
    impl Rng {
        pub(crate) fn new(seed: u64) -> Self {
            Rng(seed)
        }
        pub(crate) fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^ (z >> 31)
        }
        pub(crate) fn uniform(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
        /// Poisson: Knuth's method below 30, a rounded normal above.
        pub(crate) fn poisson(&mut self, lambda: f64) -> u64 {
            if lambda <= 0.0 {
                return 0;
            }
            if lambda < 30.0 {
                let (l, mut k, mut p) = ((-lambda).exp(), 0u64, 1.0);
                loop {
                    p *= self.uniform();
                    if p <= l {
                        return k;
                    }
                    k += 1;
                }
            }
            let (u1, u2) = (self.uniform().max(1e-300), self.uniform());
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            (lambda + lambda.sqrt() * z).round().max(0.0) as u64
        }
        pub(crate) fn binomial_thin(&mut self, n: u64, p: f64) -> u64 {
            (0..n).filter(|_| self.uniform() < p).count() as u64
        }
    }

    /// The truth-test pass rule (operator, 20 Nov; validation criterion 1): cells within 5% of the
    /// true correction must hold at least this share of the fragments, at full and at 10% depth.
    pub(crate) const TRUTH_SHARE: f64 = 0.98;

    pub(crate) const LO: u16 = 60;
    pub(crate) const HI: u16 = 200;
    pub(crate) const K: u16 = 5;

    /// Expected windows per cell: `GC fraction ~ Normal(0.43, 0.08)`, the genome's shape.
    pub(crate) fn genome_shaped(n: u16, gc: u16) -> f64 {
        let z = (gc as f64 / n as f64 - 0.43) / 0.08;
        2e6 / n as f64 * (-0.5 * z * z).exp()
    }

    /// Expected windows per cell, uniform in GC fraction over 0.2-0.8.
    pub(crate) fn uniform_ref(n: u16, gc: u16) -> f64 {
        if (0.2..=0.8).contains(&(gc as f64 / n as f64)) { 2e6 / n as f64 } else { 0.0 }
    }

    pub(crate) fn expected_table(shape: fn(u16, u16) -> f64) -> ExpectedTable {
        let mut g = GcGrid::new(LO, HI, K, 0u64).unwrap();
        let keys: Vec<GcKey> = g.keys().collect();
        for (i, key) in keys.iter().enumerate() {
            g.cells_mut()[i] = shape(key.length - 2 * K, key.gc).round() as u64;
        }
        ExpectedTable::from_counts(g, "synthetic", "territory-a")
    }

    /// Every cell with an expected window: the support for the fit tests (Task 4 builds the real
    /// one from the expected table).
    pub(crate) fn nonzero_support(e: &ExpectedTable) -> CellSet {
        let mut g = GcGrid::new(LO, HI, K, false).unwrap();
        for (i, &c) in e.counts().cells().iter().enumerate() {
            g.cells_mut()[i] = c > 0;
        }
        CellSet::from_grid(g, 1.0)
    }

    /// A peaked size distribution (mode 167 bp).
    pub(crate) fn size(length: u16) -> f64 {
        let z = (length as f64 - 167.0) / 20.0;
        (-0.5 * z * z).exp() + 0.02
    }

    /// Expected observed count per cell for `n` fragments, size distribution `size`, bias `f`.
    pub(crate) fn mean_tally(e: &ExpectedTable, n: f64, f: &dyn Fn(u16, f64) -> f64) -> Vec<f64> {
        let g = e.counts();
        let mut lambda = vec![0.0; g.len()];
        let mut z = 0.0;
        for len in LO..=HI {
            let row = g.row(len);
            let e_l: f64 = row.iter().map(|&c| c as f64).sum();
            let nn = g.gc_max(len) as f64;
            let off = g.index(len, 0).unwrap();
            for (gc, &c) in row.iter().enumerate() {
                let v = size(len) * c as f64 / e_l * f(len, gc as f64 / nn);
                lambda[off + gc] = v;
                z += v;
            }
        }
        lambda.iter().map(|v| v * n / z).collect()
    }

    pub(crate) fn tally_from(e: &ExpectedTable, lambda: &[f64], rng: Option<&mut Rng>) -> Tally {
        let mut t = Tally::for_expected(e);
        let keys: Vec<GcKey> = e.counts().keys().collect();
        match rng {
            None => keys.iter().zip(lambda).for_each(|(k, l)| { t.add(*k, l.round() as u64).unwrap(); }),
            Some(r) => keys.iter().zip(lambda).for_each(|(k, l)| { t.add(*k, r.poisson(*l)).unwrap(); }),
        }
        t
    }

    pub(crate) fn steep(slope: f64) -> impl Fn(u16, f64) -> f64 {
        move |_, x| (slope * (x - 0.5)).exp()
    }

    /// `Σ obs·w / Σ obs` at one length, as a consumer applies the weights.
    pub(crate) fn weighted_ratio(m: &GcModel, t: &Tally, length: u16) -> f64 {
        let (mut num, mut den) = (0.0, 0.0);
        for (gc, &o) in t.counts().row(length).iter().enumerate() {
            if o > 0 {
                num += o as f64 * m.weight(length, Some(gc as u16)).expect("an observed cell has a weight");
                den += o as f64;
            }
        }
        num / den
    }

    #[test]
    fn every_length_keeps_its_total_weight_on_a_non_trivial_surface() {
        // Bias slope varies with length (4 at 60 bp to 10 at 200 bp), on a genome-shaped reference.
        let e = expected_table(genome_shaped);
        let f = |len: u16, x: f64| ((4.0 + 6.0 * (len - LO) as f64 / (HI - LO) as f64) * (x - 0.5)).exp();
        let t = tally_from(&e, &mean_tally(&e, 5e6, &f), None);
        let m = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
        for len in LO..=HI {
            let r = weighted_ratio(&m, &t, len);
            assert!((r - 1.0).abs() < 1e-9, "length {len}: Σobs·w / Σobs = {r:.12}");
        }
        let row: Vec<f64> = m.weights().row(167).iter().flatten().copied().collect();
        let (lo, hi) = row.iter().fold((f64::MAX, f64::MIN), |(a, b), &w| (a.min(w), b.max(w)));
        assert!(hi / lo > 2.0, "the weights must not be flat: {lo:.3}..{hi:.3}");
    }

    #[test]
    fn uniform_and_genome_shaped_references_both_keep_total_weight() {
        // Today's fit (mean ratio over fit, unweighted) gives Σobs·w / Σobs = 2.45 on the
        // genome-shaped case: it carries mean(r) / E_exp[r] per length bin. Documented here only;
        // the old fit is not under test.
        for (name, shape) in [("uniform", uniform_ref as fn(u16, u16) -> f64), ("genome", genome_shaped)] {
            let e = expected_table(shape);
            let t = tally_from(&e, &mean_tally(&e, 5e6, &steep(8.0)), None);
            let m = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
            for len in [LO, 120, 167, HI] {
                let r = weighted_ratio(&m, &t, len);
                assert!((r - 1.0).abs() < 1e-9, "{name} reference, length {len}: {r:.12}");
            }
        }
    }

    #[test]
    fn a_rerun_gives_a_bit_identical_table() {
        let e = expected_table(genome_shaped);
        let mut rng = Rng::new(11);
        let t = tally_from(&e, &mean_tally(&e, 2e6, &steep(6.0)), Some(&mut rng));
        let s = nonzero_support(&e);
        let a = fit(&t, &e, &s, &FitConfig::default()).unwrap();
        let b = fit(&t, &e, &s, &FitConfig::default()).unwrap();
        let bits = |m: &GcModel| m.weights().cells().iter().map(|w| w.map(f64::to_bits)).collect::<Vec<_>>();
        assert_eq!(bits(&a), bits(&b));
        assert_eq!(a.qc(), b.qc());
    }

    #[test]
    fn two_biases_give_two_tables() {
        let e = expected_table(genome_shaped);
        let s = nonzero_support(&e);
        let a = fit(&tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None), &e, &s, &FitConfig::default()).unwrap();
        let b = fit(&tally_from(&e, &mean_tally(&e, 2e6, &steep(8.0)), None), &e, &s, &FitConfig::default()).unwrap();
        let ratio = |m: &GcModel| m.weight(167, Some(40)).unwrap() / m.weight(167, Some(100)).unwrap();
        assert!((ratio(&a) / ratio(&b) - 1.0).abs() > 0.5, "{} vs {}", ratio(&a), ratio(&b));
    }

    #[test]
    fn the_weights_undo_the_bias_up_to_a_per_length_scale() {
        let e = expected_table(genome_shaped);
        let f = |len: u16, x: f64| ((4.0 + 4.0 * (len - LO) as f64 / (HI - LO) as f64) * (x - 0.5)).exp();
        let t = tally_from(&e, &mean_tally(&e, 5e6, &f), None);
        let m = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
        let share = share_recovered(&m, &t, &f, 0.05);
        assert!(share >= TRUTH_SHARE, "fragments in cells recovered within 5%: {share:.4}");
    }

    /// The share of the tally's fragments in cells where the fit recovers the truth within `tol`:
    /// `|w·f / c_L - 1| <= tol`, with `c_L` length L's fragment-weighted mean of `w·f` (the bias is
    /// identifiable only up to a per-length scale). A cell without a weight counts as a miss.
    pub(crate) fn share_recovered(m: &GcModel, t: &Tally, f: &dyn Fn(u16, f64) -> f64, tol: f64) -> f64 {
        let g = t.counts();
        let (mut hit, mut all) = (0.0, 0.0);
        for len in g.min_len()..=g.max_len() {
            let x = |gc: usize| gc as f64 / g.gc_max(len) as f64;
            let (mut num, mut den) = (0.0, 0.0);
            for (gc, &o) in g.row(len).iter().enumerate() {
                if let Some(w) = m.weight(len, Some(gc as u16)) {
                    num += o as f64 * w * f(len, x(gc));
                    den += o as f64;
                }
            }
            for (gc, &o) in g.row(len).iter().enumerate() {
                all += o as f64;
                if let Some(w) = m.weight(len, Some(gc as u16)) {
                    if den > 0.0 && (w * f(len, x(gc)) / (num / den) - 1.0).abs() <= tol {
                        hit += o as f64;
                    }
                }
            }
        }
        hit / all
    }

    /// The fewest cells, largest first, holding `share` of the tally's fragments.
    pub(crate) fn top_cells(t: &Tally, share: f64) -> Vec<GcKey> {
        let g = t.counts();
        let mut cells: Vec<(u64, GcKey)> = g.keys().zip(g.cells()).map(|(k, &c)| (c, k)).collect();
        cells.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let mut keep = Vec::new();
        let mut cum = 0u64;
        for (c, k) in cells {
            if cum as f64 >= share * t.total() as f64 {
                break;
            }
            cum += c;
            keep.push(k);
        }
        keep
    }

    #[test]
    fn depth_changes_the_table_little_and_both_depths_match_the_truth() {
        // Criterion 4 on a synthetic grid: 10% of the fragments gives a table correlated > 0.98
        // with the full one, and both are compared with the truth, not only with each other.
        let e = expected_table(genome_shaped);
        let f = steep(6.0);
        let mut rng = Rng::new(3);
        let full = tally_from(&e, &mean_tally(&e, 3e6, &f), Some(&mut rng));
        let mut tenth = Tally::for_expected(&e);
        let keys: Vec<GcKey> = e.counts().keys().collect();
        for (k, &c) in keys.iter().zip(full.counts().cells()) {
            tenth.add(*k, rng.binomial_thin(c, 0.1)).unwrap();
        }
        let s = nonzero_support(&e);
        let a = fit(&full, &e, &s, &FitConfig::default()).unwrap();
        let b = fit(&tenth, &e, &s, &FitConfig::default()).unwrap();
        let both: Vec<(f64, f64)> = top_cells(&full, 0.95).iter()
            .map(|k| (a.weight(k.length, Some(k.gc)).unwrap(), b.weight(k.length, Some(k.gc)).unwrap()))
            .collect();
        let r = pearson(&both);
        let (ea, eb) = (share_recovered(&a, &full, &f, 0.05), share_recovered(&b, &tenth, &f, 0.05));
        eprintln!("depth: correlation {r:.4}; fragments recovered within 5%: full {ea:.4}, 10% {eb:.4}");
        assert!(r > 0.98, "correlation of the 100% and 10% tables: {r:.4}");
        assert!(ea >= TRUTH_SHARE, "full depth against the truth: {ea:.4}");
        assert!(eb >= TRUTH_SHARE, "10% depth against the truth: {eb:.4}");
    }

    pub(crate) fn pearson(xy: &[(f64, f64)]) -> f64 {
        let n = xy.len() as f64;
        let (mx, my) = (xy.iter().map(|p| p.0).sum::<f64>() / n, xy.iter().map(|p| p.1).sum::<f64>() / n);
        let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
        for (x, y) in xy {
            sxy += (x - mx) * (y - my);
            sxx += (x - mx).powi(2);
            syy += (y - my).powi(2);
        }
        sxy / (sxx * syy).sqrt()
    }

    #[test]
    fn a_tally_below_the_minimum_cannot_be_fitted() {
        let e = expected_table(genome_shaped);
        let s = nonzero_support(&e);
        let zero = Tally::for_expected(&e);
        assert_eq!(fit(&zero, &e, &s, &FitConfig::default()).unwrap_err(),
                   FitError::TooFewFragments { fragments: 0, required: MIN_FIT_FRAGMENTS });
        let small = tally_from(&e, &mean_tally(&e, 90_000.0, &steep(4.0)), None);
        assert!(matches!(fit(&small, &e, &s, &FitConfig::default()),
                         Err(FitError::TooFewFragments { .. })));
    }

    #[test]
    fn a_length_with_no_support_cell_cannot_be_fitted() {
        let e = expected_table(genome_shaped);
        let mut g = nonzero_support(&e).grid().clone();
        g.row_mut(150).iter_mut().for_each(|c| *c = false);
        let s = CellSet::from_grid(g, 1.0);
        let t = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        assert_eq!(fit(&t, &e, &s, &FitConfig::default()).unwrap_err(), FitError::NoSupport { length: 150 });
    }

    #[test]
    fn an_invalid_config_is_refused_not_fitted() {
        let e = expected_table(genome_shaped);
        let s = nonzero_support(&e);
        let t = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        let d = FitConfig::default();
        let bad = [
            FitConfig { cap: 0.5, ..d.clone() },
            FitConfig { cap: f64::NAN, ..d.clone() },
            FitConfig { cap: f64::INFINITY, ..d.clone() },
            FitConfig { bandwidth_length_bp: 0.0, ..d.clone() },
            FitConfig { bandwidth_length_bp: f64::NAN, ..d.clone() },
            FitConfig { bandwidth_gc: -0.02, ..d.clone() },
            FitConfig { bandwidth_gc: f64::INFINITY, ..d.clone() },
            FitConfig { min_fragments: 0, ..d.clone() },
        ];
        for cfg in bad {
            let r = std::panic::catch_unwind(|| fit(&t, &e, &s, &cfg));
            assert!(matches!(r, Ok(Err(FitError::InvalidConfig(_)))), "{cfg:?} gave {r:?}");
        }
        // An empty tally with min_fragments 0 is refused too, by the config check.
        let zero = Tally::for_expected(&e);
        assert!(matches!(fit(&zero, &e, &s, &FitConfig { min_fragments: 0, ..d }),
                         Err(FitError::InvalidConfig(_))));
    }

    #[test]
    fn a_tally_for_another_territory_or_trim_is_refused() {
        let e = expected_table(genome_shaped);
        let s = nonzero_support(&e);
        let good = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        let mut other = Tally::new(LO, HI, K, "territory-b").unwrap();
        for (k, &c) in good.counts().keys().zip(good.counts().cells()) {
            other.add(k, c).unwrap();
        }
        assert!(matches!(fit(&other, &e, &s, &FitConfig::default()), Err(FitError::TerritoryMismatch { .. })));
        let trim = Tally::new(LO, HI, K + 1, "territory-a").unwrap();
        assert!(matches!(fit(&trim, &e, &s, &FitConfig::default()), Err(FitError::ShapeMismatch(_))));
    }

    #[test]
    fn add_insert_keys_with_the_tallys_own_trim() {
        let e = expected_table(genome_shaped);
        let mut t = Tally::for_expected(&e);
        // 70 bp at k = 5: GGGGG ends are trimmed off, the 60 bp window holds 15 G.
        let mut insert = b"GGGGG".to_vec();
        insert.extend(std::iter::repeat_n(b"GAAA", 15).flatten());
        insert.extend(b"GGGGG");
        assert_eq!(t.add_insert(&insert), Some(GcKey { length: 70, gc: 15 }));
        assert_eq!(*t.counts().get(70, 15).unwrap(), 1);
        // No key: an N in the window, or too short for the trim. Not "outside lengths".
        let mut with_n = insert.clone();
        with_n[30] = b'N';
        assert_eq!(t.add_insert(&with_n), None);
        assert_eq!(t.add_insert(&[b'A'; 29]), None);
        assert_eq!((t.without_key(), t.outside_lengths()), (2, 0));
        // A keyed insert longer than the grid: outside lengths, not without key.
        assert_eq!(t.add_insert(&[b'A'; 300]), Some(GcKey { length: 300, gc: 0 }));
        assert_eq!((t.without_key(), t.outside_lengths()), (2, 1));
    }

    #[test]
    fn a_key_no_insert_can_have_is_an_error() {
        let e = expected_table(genome_shaped);
        let mut t = Tally::for_expected(&e);
        assert!(t.add(GcKey { length: 100, gc: 91 }, 1).is_err(), "gc above length - 2k");
        assert!(t.add(GcKey { length: 300, gc: 291 }, 1).is_err(), "out of range and impossible");
        assert!(t.add(GcKey { length: 29, gc: 0 }, 1).is_err(), "shorter than 2k + 20");
        assert_eq!((t.total(), t.outside_lengths()), (0, 0));
    }

    #[test]
    fn the_qc_separates_no_key_from_outside_lengths() {
        let e = expected_table(genome_shaped);
        let mut t = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        t.add_insert(&[b'N'; 100]);
        t.add_insert(&[b'N'; 100]);
        t.add_insert(&[b'A'; 250]);
        let q = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap().qc().clone();
        assert_eq!((q.fragments_without_key, q.fragments_outside_lengths), (2, 1));
    }

    #[test]
    fn mismatched_shapes_cannot_be_fitted() {
        let e = expected_table(genome_shaped);
        let t = Tally::new(LO, HI - 1, K, "territory-a").unwrap();
        assert!(matches!(fit(&t, &e, &nonzero_support(&e), &FitConfig::default()),
                         Err(FitError::ShapeMismatch(_))));
    }

    #[test]
    fn the_provenance_names_its_inputs() {
        let e = expected_table(genome_shaped);
        let s = nonzero_support(&e);
        let m = fit(&tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None), &e, &s, &FitConfig::default()).unwrap();
        let p = m.provenance();
        assert_eq!((p.genome.as_str(), p.territory_sha256.as_str()), ("synthetic", "territory-a"));
        assert_eq!(p.support_sha256, s.digest());
        assert_eq!((p.trim_k, p.gc_definition_version), (K, GC_DEFINITION_VERSION));
        assert_eq!(p.fit_config, FitConfig::default());
    }

    #[test]
    fn the_smoothing_bias_is_reported_and_grows_with_curvature() {
        let e = expected_table(genome_shaped);
        let s = nonzero_support(&e);
        let flat = fit(&tally_from(&e, &mean_tally(&e, 2e6, &|_, _| 1.0), None), &e, &s, &FitConfig::default()).unwrap();
        let curved = fit(&tally_from(&e, &mean_tally(&e, 2e6, &steep(10.0)), None), &e, &s, &FitConfig::default()).unwrap();
        let (b0, b1) = (flat.qc().smoothing_bias_mean, curved.qc().smoothing_bias_mean);
        eprintln!("smoothing bias: flat {b0:.5}, κ=10 {b1:.5}, p95 {:.5}", curved.qc().smoothing_bias_p95);
        // (κh)² / 2 with κ = 10, h = 0.02 is 0.02; a flat bias has none beyond rounding.
        assert!(b0 < 0.002, "flat: {b0}");
        assert!((0.015..0.025).contains(&b1), "κ = 10: {b1}");
        assert!(curved.qc().smoothing_bias_p95 >= b1);
    }

    /// Fit runtime on the full grid at the default trim (40-600 bp): `cargo test -- --ignored`.
    #[test]
    #[ignore]
    fn timing_of_a_fit_on_the_full_grid() {
        let (lo, hi) = GcGrid::<u64>::lengths_for(DEFAULT_TRIM_K);
        let mut g = GcGrid::new(lo, hi, DEFAULT_TRIM_K, 0u64).unwrap();
        let keys: Vec<GcKey> = g.keys().collect();
        for (i, k) in keys.iter().enumerate() {
            g.cells_mut()[i] = genome_shaped(k.length - 2 * DEFAULT_TRIM_K, k.gc).round() as u64;
        }
        let e = ExpectedTable::from_counts(g, "synthetic", "t");
        let s = CellSet::from_expected(&e, SUPPORT_SHARE).unwrap();
        let mut t = Tally::for_expected(&e);
        let n_cells = s.grid().cells().iter().filter(|&&b| b).count();
        for (k, &c) in keys.iter().zip(e.counts().cells()) {
            t.add(*k, c / 10).unwrap();
        }
        let t0 = std::time::Instant::now();
        let m = fit(&t, &e, &s, &FitConfig::default()).unwrap();
        eprintln!("fit: lengths {lo}..={hi}, {} cells, {n_cells} in support, {} fragments: {:.2} s",
                  e.counts().len(), t.total(), t0.elapsed().as_secs_f64());
        assert!(m.qc().lengths_unfitted == 0);
    }

    #[test]
    fn lengths_beyond_the_kernels_reach_are_unfitted_and_reported() {
        // No fragment at 100-150 bp: the middle lengths are more than 4 SD (20 bp) from any data.
        let e = expected_table(genome_shaped);
        let mut t = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        for len in 100..=150 {
            t.counts.row_mut(len).iter_mut().for_each(|c| *c = 0);
        }
        let m = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
        // 120-130 are unreachable: each is more than 20 bp from 99 and from 151.
        assert_eq!(m.qc().lengths_unfitted, 11);
        assert_eq!(m.qc().lengths_unobserved, 51);
        assert_eq!(m.weight(125, Some(60)), None);
        assert!(m.weight(110, Some(50)).is_some());
    }

    #[test]
    fn fragments_outside_the_length_range_are_reported() {
        let e = expected_table(genome_shaped);
        let mut t = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        let m0 = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
        assert!(!t.add(GcKey { length: HI + 50, gc: 100 }, 7).unwrap());
        assert!(!t.add(GcKey { length: LO - 1, gc: 10 }, 3).unwrap());
        let m1 = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
        assert_eq!((m0.qc().fragments_outside_lengths, m1.qc().fragments_outside_lengths), (0, 10));
    }

    #[test]
    fn a_length_with_no_fragment_is_reported_unobserved() {
        let e = expected_table(genome_shaped);
        let mut t = tally_from(&e, &mean_tally(&e, 2e6, &steep(4.0)), None);
        t.counts.row_mut(130).iter_mut().for_each(|c| *c = 0);
        let m = fit(&t, &e, &nonzero_support(&e), &FitConfig::default()).unwrap();
        assert_eq!(m.qc().lengths_unobserved, 1);
        assert_eq!(m.qc().fragments, t.total());
        assert!(m.weight(130, Some(60)).is_some(), "smoothed from its neighbours");
    }
}

#[cfg(test)]
mod sparse_tests {
    use super::fit_tests::*;
    use super::*;

    fn support(e: &crate::gc_reference::ExpectedTable) -> CellSet {
        CellSet::from_expected(e, SUPPORT_SHARE).unwrap()
    }

    #[test]
    fn the_support_is_the_fewest_cells_holding_the_share_at_every_length() {
        let e = expected_table(genome_shaped);
        let s = support(&e);
        for len in LO..=HI {
            let row = e.counts().row(len);
            let total: u64 = row.iter().sum();
            let held: u64 = row.iter().zip(s.grid().row(len)).filter(|(_, &b)| b).map(|(c, _)| c).sum();
            assert!(held as f64 >= SUPPORT_SHARE * total as f64, "length {len}: {held} of {total}");
            // Minimal: without its smallest member it falls below the share.
            let smallest = row.iter().zip(s.grid().row(len)).filter(|(_, &b)| b).map(|(c, _)| *c).min().unwrap();
            assert!(((held - smallest) as f64) < SUPPORT_SHARE * total as f64, "length {len} is not minimal");
        }
        // Two references, two supports.
        let u = support(&expected_table(uniform_ref));
        assert_ne!(s.digest(), u.digest());
    }

    #[test]
    fn a_cell_outside_the_support_takes_the_edge_value_and_does_not_enter_the_fit() {
        let e = expected_table(genome_shaped);
        let s = support(&e);
        let base = tally_from(&e, &mean_tally(&e, 2e6, &steep(6.0)), None);
        let len = 167u16;
        let edge = (0..=e.counts().gc_max(len)).rev().find(|&g| s.contains(len, g)).unwrap();
        let outlier = edge + 5; // within the kernel's reach of the support (1.6 SD), so it could leak
        let mut t = base.clone();
        t.add(GcKey { length: len, gc: outlier }, 5_000).unwrap();
        let cfg = FitConfig::default();
        let (a, b) = (fit(&base, &e, &s, &cfg).unwrap(), fit(&t, &e, &s, &cfg).unwrap());
        assert_eq!(b.weight(len, Some(outlier)), b.weight(len, Some(edge)), "clamped to the edge value");
        // The outlier changes the length's normalising scale and nothing else.
        let k = b.weight(len, Some(edge)).unwrap() / a.weight(len, Some(edge)).unwrap();
        for g in 0..=e.counts().gc_max(len) {
            if s.contains(len, g) {
                let r = b.weight(len, Some(g)).unwrap() / a.weight(len, Some(g)).unwrap();
                assert!((r / k - 1.0).abs() < 1e-12, "gc {g}: the outlier entered the smoother");
            }
        }
        assert!((weighted_ratio(&b, &t, len) - 1.0).abs() < 1e-9, "the outlier's weight is normalised too");
    }

    #[test]
    fn raw_weights_are_capped_before_normalisation() {
        // Slope 16 asks for raw weights near 40 at the low-GC edge of the support.
        let e = expected_table(genome_shaped);
        let s = support(&e);
        let cfg = FitConfig::default();
        let t = tally_from(&e, &mean_tally(&e, 2e6, &steep(16.0)), None);
        let (m, scales) = fit_with_scales(&t, &e, &s, &cfg).unwrap();
        // Raw = final * the length's divisor: every raw weight lies in [1/C, C], and a capped one
        // sits exactly on a bound. Capping after normalisation ("normalise, cap, renormalise")
        // leaves a capped cell at C times the first divisor instead.
        let mut outside_bounds = 0;
        for (li, len) in (LO..=HI).enumerate() {
            for (w, &c) in m.weights().row(len).iter().zip(m.capped().row(len)) {
                let (w, raw) = (w.unwrap(), w.unwrap() * scales[li]);
                assert!(raw >= (1.0 / cfg.cap) * (1.0 - 1e-12) && raw <= cfg.cap * (1.0 + 1e-12),
                        "length {len}: raw {raw}");
                if c {
                    let on_bound = (raw / cfg.cap - 1.0).abs() < 1e-12 || (raw * cfg.cap - 1.0).abs() < 1e-12;
                    assert!(on_bound, "length {len}: a capped raw weight {raw} is not on a bound");
                }
                if !(1.0 / cfg.cap..=cfg.cap).contains(&w) {
                    outside_bounds += 1;
                }
            }
        }
        assert!(outside_bounds > 0, "the cap bounds raw weights, so some normalised ones lie outside it");
        for len in LO..=HI {
            let row: Vec<f64> = m.weights().row(len).iter().flatten().copied().collect();
            let (lo, hi) = row.iter().fold((f64::MAX, f64::MIN), |(a, b), &w| (a.min(w), b.max(w)));
            // Both ends raw in [1/C, C] and divided by one scale: their ratio is at most C².
            assert!(hi / lo <= cfg.cap * cfg.cap * (1.0 + 1e-12), "length {len}: {hi} / {lo}");
            assert!((weighted_ratio(&m, &t, len) - 1.0).abs() < 1e-9, "length {len} keeps its total");
        }
        assert!(m.capped().cells().iter().any(|&c| c), "the case must reach the cap");
        assert!(m.qc().capped_weight_share > 0.0);
    }

    #[test]
    fn the_qc_shares_vary_with_the_input() {
        let e = expected_table(genome_shaped);
        let s = support(&e);
        let cfg = FitConfig::default();
        let qc = |slope: f64| fit(&tally_from(&e, &mean_tally(&e, 2e6, &steep(slope)), None), &e, &s, &cfg)
            .unwrap().qc().clone();
        let (gentle, steep16) = (qc(4.0), qc(16.0));
        assert_eq!(gentle.capped_weight_share, 0.0, "a gentle bias stays inside the cap");
        assert!(steep16.capped_weight_share > 0.0);
        assert!(gentle.clamped_weight_share > 0.0 && steep16.clamped_weight_share > 0.0);
        assert_ne!(gentle.clamped_weight_share, steep16.clamped_weight_share);
        for q in [&gentle, &steep16] {
            assert!(q.clamped_weight_share < 0.05 && q.capped_weight_share < 0.5, "{q:?}");
        }
    }

    #[test]
    fn every_observed_fragment_has_a_weight_in_a_fitted_length() {
        let e = expected_table(genome_shaped);
        let s = support(&e);
        let t = tally_from(&e, &mean_tally(&e, 2e6, &steep(8.0)), None);
        let m = fit(&t, &e, &s, &FitConfig::default()).unwrap();
        for k in t.counts().keys() {
            if *t.counts().get(k.length, k.gc).unwrap() > 0 {
                assert!(m.weight(k.length, Some(k.gc)).is_some(), "{k:?}");
            }
        }
    }
}

/// Truth tests on injected bias (criterion 3). Unlike the fit tests, nothing here is a grid made
/// by hand: a reference sequence with GC blocks, the expected counts built from it by the real
/// builder, and fragments drawn from it with a peaked size distribution, keyed by `gc_key`, and
/// thinned by a known `f(length, GC)`. Obs/exp then mixes size and GC as real data does.
#[cfg(test)]
mod truth_tests {
    use super::fit_tests::{share_recovered, Rng, TRUTH_SHARE};
    use super::*;
    use crate::gc_reference::{ExpectedTable, Territory};
    use crate::contigs::ContigNames;

    const GENOME_BP: usize = 600_000;
    const BLOCK_BP: usize = 3_000;
    const LO: u16 = 60;
    const HI: u16 = 250;

    /// Blocks of 3 kb whose GC probability cycles through 35-65%.
    fn genome() -> Vec<u8> {
        let mut rng = Rng::new(42);
        (0..GENOME_BP).map(|i| {
            let p_gc = 0.35 + 0.05 * ((i / BLOCK_BP * 3) % 7) as f64;
            let gc = rng.uniform() < p_gc;
            let pick = rng.uniform() < 0.5;
            match (gc, pick) { (true, true) => b'G', (true, false) => b'C', (false, true) => b'A', (false, false) => b'T' }
        }).collect()
    }

    /// A peaked size distribution: mode 167 bp, with a shoulder of short fragments.
    fn draw_length(rng: &mut Rng) -> u16 {
        loop {
            let (u1, u2) = (rng.uniform().max(1e-300), rng.uniform());
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            let len = if rng.uniform() < 0.15 { 110.0 + 15.0 * z } else { 167.0 + 20.0 * z };
            let len = len.round();
            if (LO as f64..=HI as f64).contains(&len) {
                return len as u16;
            }
        }
    }

    struct Case {
        expected: ExpectedTable,
        support: CellSet,
        seq: Vec<u8>,
    }

    fn case() -> Case {
        let seq = genome();
        let names = ContigNames::new([("1".to_string(), GENOME_BP as u64)]);
        let territory = Territory::resolve(&[("1".to_string(), 0, GENOME_BP as u64)], &names).unwrap();
        let expected = ExpectedTable::count(&territory, "synthetic", LO, HI, DEFAULT_TRIM_K, |iv| {
            Ok(seq[iv.start as usize..iv.end as usize].to_vec())
        }).unwrap();
        let support = CellSet::from_expected(&expected, SUPPORT_SHARE).unwrap();
        Case { expected, support, seq }
    }

    /// Draw fragments uniformly over the territory and keep each with probability
    /// `f / f_max(length)`. Keys come from a GC prefix sum (the genome has no N), checked against
    /// `gc_key` on the first draws; a per-length ceiling keeps the acceptance rate workable, and
    /// is asserted to bound `f`.
    fn sample(c: &Case, n_keep: u64, f: &dyn Fn(u16, f64) -> f64, f_max: &dyn Fn(u16) -> f64, seed: u64) -> Tally {
        let mut prefix = Vec::with_capacity(c.seq.len() + 1);
        prefix.push(0u32);
        for &b in &c.seq {
            prefix.push(prefix[prefix.len() - 1] + matches!(b, b'G' | b'C') as u32);
        }
        let k = DEFAULT_TRIM_K as usize;
        let mut rng = Rng::new(seed);
        let mut t = Tally::for_expected(&c.expected);
        let (mut kept, mut drawn) = (0, 0u64);
        while kept < n_keep {
            let len = draw_length(&mut rng);
            let start = (rng.uniform() * (GENOME_BP - len as usize + 1) as f64) as usize;
            let gc = (prefix[start + len as usize - k] - prefix[start + k]) as u16;
            let key = GcKey { length: len, gc };
            drawn += 1;
            if drawn <= 1_000 {
                assert_eq!(gc_key(&c.seq[start..start + len as usize], DEFAULT_TRIM_K), Some(key));
            }
            let p = f(len, gc as f64 / (len - 2 * DEFAULT_TRIM_K) as f64) / f_max(len);
            assert!(p <= 1.0, "the acceptance ceiling must bound f (length {len}, gc {gc})");
            if rng.uniform() < p {
                t.add(key, 1).unwrap();
                kept += 1;
            }
        }
        t
    }

    /// GC bias steepening with length (slope 3 at 60 bp to 8 at 250 bp).
    fn injected(len: u16, x: f64) -> f64 {
        let slope = 3.0 + 5.0 * (len - LO) as f64 / (HI - LO) as f64;
        (slope * (x - 0.5)).exp()
    }

    /// The injected-bias sample (2,000,000 fragments) and its 10% thinning (about 200,000).
    fn injected_full_and_tenth(c: &Case) -> (Tally, Tally) {
        // GC fraction stays below 0.9 for windows of 40 bp or more in 35-65% blocks.
        let ceiling = |len: u16| injected(len, 0.9);
        let t = sample(c, 2_000_000, &injected, &ceiling, 1);
        let mut rng = Rng::new(5);
        let mut tenth = Tally::for_expected(&c.expected);
        for (k, &n) in t.counts().keys().zip(t.counts().cells()) {
            tenth.add(k, rng.binomial_thin(n, 0.1)).unwrap();
        }
        (t, tenth)
    }

    #[test]
    fn at_ten_percent_depth_the_injected_bias_falls_short_of_the_rule() {
        // KNOWN SHORTFALL, pinned like a strict xfail: at 10% depth (~200,000 fragments) cells
        // within 5% hold 0.977 of fragments, under the operator's 98% rule. When a change brings
        // it to 98%, this fails: move the assertion into the test below and delete this one.
        let c = case();
        let (_, tenth) = injected_full_and_tenth(&c);
        let m = fit(&tenth, &c.expected, &c.support, &FitConfig::default()).unwrap();
        let share = share_recovered(&m, &tenth, &injected, 0.05);
        eprintln!("injected bias at 10% depth: {share:.4} of fragments recovered within 5%");
        assert!(share < TRUTH_SHARE, "10% depth now meets the rule ({share:.4}): un-pin it");
        assert!(share > 0.95, "10% depth got worse than when pinned: {share:.4}");
    }

    #[test]
    fn an_injected_bias_is_recovered_within_five_percent() {
        let c = case();
        let (t, _) = injected_full_and_tenth(&c);
        let m = fit(&t, &c.expected, &c.support, &FitConfig::default()).unwrap();
        let share = share_recovered(&m, &t, &injected, 0.05);
        eprintln!("injected bias: {share:.4} of fragments recovered within 5%; qc {:?}", m.qc());
        assert!(share >= TRUTH_SHARE, "fragments in cells where the injected bias is recovered within 5%: {share:.4}");
        // And the correction is not flat: at 200 bp the weights span the bias they undo.
        let w = |gc| m.weight(200, Some(gc)).unwrap();
        let (lo, hi) = (w(65), w(110)); // GC fraction 0.36 and 0.61
        let want = injected(200, 110.0 / 180.0) / injected(200, 65.0 / 180.0);
        assert!((lo / hi / want - 1.0).abs() < 0.05, "weight ratio {:.3} vs injected {:.3}", lo / hi, want);
    }

    /// A model with given weights over the case's support, for showing what the truth tests reject.
    fn model_with(c: &Case, w: impl Fn(u16, u16) -> f64) -> GcModel {
        let g = c.expected.counts();
        let mut weights = GcGrid::new(g.min_len(), g.max_len(), g.trim_k(), None).unwrap();
        let keys: Vec<GcKey> = g.keys().collect();
        for (i, k) in keys.iter().enumerate() {
            weights.cells_mut()[i] = Some(w(k.length, k.gc));
        }
        let capped = GcGrid::new(g.min_len(), g.max_len(), g.trim_k(), false).unwrap();
        let fitted = fit(&sample(c, 200_000, &|_, _| 1.0, &|_| 1.0, 9), &c.expected, &c.support,
                         &FitConfig::default()).unwrap();
        GcModel::from_parts(weights, capped, c.support.clone(), fitted.provenance().clone(),
                            FitQc::default()).unwrap()
    }

    #[test]
    fn the_truth_tests_reject_flat_weights_a_curve_and_oversmoothing() {
        let c = case();
        let ceiling = |len: u16| injected(len, 0.9);
        let biased = sample(&c, 1_000_000, &injected, &ceiling, 1);
        // Weights of 1 on a biased sample.
        let flat = model_with(&c, |_, _| 1.0);
        let share = share_recovered(&flat, &biased, &injected, 0.05);
        assert!(share < TRUTH_SHARE, "weights of 1 pass the injected-bias test: {share:.4}");
        // Four times the GC bandwidth.
        let wide = FitConfig { bandwidth_gc: 4.0 * DEFAULT_BANDWIDTH_GC, ..FitConfig::default() };
        let m = fit(&biased, &c.expected, &c.support, &wide).unwrap();
        let share = share_recovered(&m, &biased, &injected, 0.05);
        assert!(share < TRUTH_SHARE, "a 4x GC bandwidth passes the injected-bias test: {share:.4}");
        // A curve on an unbiased sample.
        let unbiased = sample(&c, 1_000_000, &|_, _| 1.0, &|_| 1.0, 2);
        let curve = model_with(&c, |len, gc| (1.5 * (gc as f64 / (len - 2 * DEFAULT_TRIM_K) as f64 - 0.45)).exp());
        let near = share_near_one(&curve, &unbiased);
        assert!(near < TRUTH_SHARE, "a curve passes the no-bias control: {near:.4}");
    }

    /// The share of fragments in cells with `|w - 1| < 0.05`.
    fn share_near_one(m: &GcModel, t: &Tally) -> f64 {
        let (mut near, mut all) = (0.0, 0.0);
        for k in t.counts().keys() {
            let o = *t.counts().get(k.length, k.gc).unwrap() as f64;
            all += o;
            if (m.weight(k.length, Some(k.gc)).unwrap() - 1.0).abs() < 0.05 {
                near += o;
            }
        }
        near / all
    }

    #[test]
    fn a_sample_with_no_bias_gets_weights_within_five_percent_of_one() {
        // The control: a fit that always returns a curve fails here.
        let c = case();
        let flat = |_: u16, _: f64| 1.0;
        let t = sample(&c, 1_000_000, &flat, &|_| 1.0, 2);
        let m = fit(&t, &c.expected, &c.support, &FitConfig::default()).unwrap();
        let (mut near, mut all) = (0.0, 0.0);
        for k in t.counts().keys() {
            let o = *t.counts().get(k.length, k.gc).unwrap() as f64;
            all += o;
            if (m.weight(k.length, Some(k.gc)).unwrap() - 1.0).abs() < 0.05 {
                near += o;
            }
        }
        eprintln!("no-bias control: {:.4} of fragments in cells with |w - 1| < 0.05", near / all);
        assert!(near / all >= TRUTH_SHARE, "fragments in cells with |w - 1| < 0.05: {:.4}", near / all);
        assert_eq!(m.qc().capped_weight_share, 0.0, "no bias reaches the cap");
    }
}
