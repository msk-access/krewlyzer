//! GC Reference Asset Generation
//!
//! Generates pre-computed GC correction assets:
//! 1. valid_regions.bed - Curated 100kb bins excluding problematic regions
//! 2. ref_genome_GC.parquet - Expected fragment counts per (length, GC)
//!
//! These assets are generated once per reference genome and shipped with krewlyzer.
//!
//! 3. [`ExpectedTable`] - exhaustive expected counts per exact (length, GC count) over a
//!    [`Territory`], for the redesigned GC model (`gc_model.rs`). Additive: the two assets above,
//!    [`LENGTH_BINS`] and [`get_length_bin_index`] stay until spec Step 2 moves their consumers.

use std::path::Path;
use std::io::{BufRead, BufReader, Write};
use std::fs::File;
use std::collections::HashMap;

use anyhow::{Result, Context};
use log::info;
use pyo3::prelude::*;
use rust_htslib::faidx;
use crate::contigs::{norm_chrom, ContigNames};

/// Length bins for GC correction (5bp bins from 60-1000bp = 188 bins)
/// 
/// Extended range captures ultra-long fragments for:
/// - Necrosis-derived cfDNA
/// - Fetal cfDNA (placental origin)
/// - Apoptosis late-stage markers
/// 
/// 5bp granularity provides finer-grained GC bias modeling for ML features.
/// Bin 0 = 60-64bp, Bin 1 = 65-69bp, ..., Bin 187 = 995-999bp
/// 
/// Note: This matches LengthBin in gc_correction.rs for consistency.
pub const LENGTH_BINS: [(u32, u32); 188] = [
    // 60-100bp (Ultra-short: bins 0-7)
    (60, 65), (65, 70), (70, 75), (75, 80), (80, 85), (85, 90), (90, 95), (95, 100),
    // 100-150bp (Core-short: bins 8-17)
    (100, 105), (105, 110), (110, 115), (115, 120), (120, 125), (125, 130), (130, 135), (135, 140),
    (140, 145), (145, 150),
    // 150-220bp (Mono-nucleosomal: bins 18-31)
    (150, 155), (155, 160), (160, 165), (165, 170), (170, 175), (175, 180),
    (180, 185), (185, 190), (190, 195), (195, 200), (200, 205), (205, 210), (210, 215), (215, 220),
    // 220-260bp (Di-nucleosomal: bins 32-39)
    (220, 225), (225, 230), (230, 235), (235, 240), (240, 245), (245, 250), (250, 255), (255, 260),
    // 260-400bp (Long: bins 40-67)
    (260, 265), (265, 270), (270, 275), (275, 280), (280, 285), (285, 290), (290, 295), (295, 300),
    (300, 305), (305, 310), (310, 315), (315, 320), (320, 325), (325, 330), (330, 335), (335, 340),
    (340, 345), (345, 350), (350, 355), (355, 360), (360, 365), (365, 370), (370, 375), (375, 380),
    (380, 385), (385, 390), (390, 395), (395, 400),
    // 400-600bp (Ultra-long part 1: bins 68-107)
    (400, 405), (405, 410), (410, 415), (415, 420), (420, 425), (425, 430), (430, 435), (435, 440),
    (440, 445), (445, 450), (450, 455), (455, 460), (460, 465), (465, 470), (470, 475), (475, 480),
    (480, 485), (485, 490), (490, 495), (495, 500), (500, 505), (505, 510), (510, 515), (515, 520),
    (520, 525), (525, 530), (530, 535), (535, 540), (540, 545), (545, 550), (550, 555), (555, 560),
    (560, 565), (565, 570), (570, 575), (575, 580), (580, 585), (585, 590), (590, 595), (595, 600),
    // 600-800bp (Ultra-long part 2: bins 108-147)
    (600, 605), (605, 610), (610, 615), (615, 620), (620, 625), (625, 630), (630, 635), (635, 640),
    (640, 645), (645, 650), (650, 655), (655, 660), (660, 665), (665, 670), (670, 675), (675, 680),
    (680, 685), (685, 690), (690, 695), (695, 700), (700, 705), (705, 710), (710, 715), (715, 720),
    (720, 725), (725, 730), (730, 735), (735, 740), (740, 745), (745, 750), (750, 755), (755, 760),
    (760, 765), (765, 770), (770, 775), (775, 780), (780, 785), (785, 790), (790, 795), (795, 800),
    // 800-1000bp (Ultra-long part 3: bins 148-187)
    (800, 805), (805, 810), (810, 815), (815, 820), (820, 825), (825, 830), (830, 835), (835, 840),
    (840, 845), (845, 850), (850, 855), (855, 860), (860, 865), (865, 870), (870, 875), (875, 880),
    (880, 885), (885, 890), (890, 895), (895, 900), (900, 905), (905, 910), (910, 915), (915, 920),
    (920, 925), (925, 930), (930, 935), (935, 940), (940, 945), (945, 950), (950, 955), (955, 960),
    (960, 965), (965, 970), (970, 975), (975, 980), (980, 985), (985, 990), (990, 995), (995, 1000),
];

/// Number of length bins (5bp bins from 60-1000bp)
pub const NUM_LENGTH_BINS: usize = 188;

/// Get the length bin index for a given fragment length (5bp bins)
/// Returns None if the length is outside the range [60, 1000)
pub fn get_length_bin_index(length: u32) -> Option<usize> {
    if !(60..1000).contains(&length) {
        return None;
    }
    // (length - 60) / 5 -> 0..187
    Some(((length - 60) / 5) as usize)
}

/// Valid genomic regions for GC correction estimation
/// Excludes problematic regions, assembly gaps, and low-mappability areas
#[derive(Debug, Clone)]
pub struct ValidRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub gc_content: f64,
}

/// Load valid regions from a BED file (plain or BGZF compressed)
/// 
/// # Arguments
/// * `bed_path` - Path to the valid_regions.bed or .bed.gz file
/// 
/// # Returns
/// * Vector of ValidRegion structs
pub fn load_valid_regions(bed_path: &Path) -> Result<Vec<ValidRegion>> {
    use noodles::bgzf;
    use std::io::BufRead;
    
    info!("Loading valid regions from {:?}", bed_path);
    
    let file = File::open(bed_path)
        .with_context(|| format!("Failed to open valid regions file: {:?}", bed_path))?;
    
    // Check if file is BGZF compressed based on extension
    let is_bgzf = bed_path.extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false);
    
    let mut regions = Vec::new();
    
    // Create appropriate reader based on compression
    if is_bgzf {
        let mut reader = bgzf::io::Reader::new(file);
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            if let Some(region) = parse_bed_line(&line) {
                regions.push(region);
            }
            line.clear();
        }
    } else {
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            if let Some(region) = parse_bed_line(&line) {
                regions.push(region);
            }
        }
    }
    
    info!("Loaded {} valid regions", regions.len());
    Ok(regions)
}

/// Parse a single BED line into a ValidRegion
fn parse_bed_line(line: &str) -> Option<ValidRegion> {
    let line = line.trim();
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 3 {
        return None;
    }
    
    let chrom = norm_chrom(fields[0]).to_string();
    let start: u64 = fields[1].parse().ok()?;
    let end: u64 = fields[2].parse().ok()?;
    
    // GC content may be in 4th column if pre-computed
    let gc_content = if fields.len() >= 4 {
        fields[3].parse().unwrap_or(0.0)
    } else {
        0.0
    };
    
    Some(ValidRegion {
        chrom,
        start,
        end,
        gc_content,
    })
}

/// Generate on-target valid regions by intersecting with target regions
/// 
/// # Arguments
/// * `valid_regions_path` - Path to existing valid_regions.bed.gz
/// * `target_regions_path` - Path to target regions BED (e.g., capture baits)
/// * `output_path` - Path to write valid_regions_ontarget.bed.gz
/// 
/// # Returns
/// * Number of on-target valid regions generated
#[pyfunction]
pub fn generate_valid_regions_ontarget(
    valid_regions_path: &str,
    target_regions_path: &str,
    output_path: &str,
) -> PyResult<usize> {
    info!("Generating on-target valid regions...");
    info!("  Valid regions: {}", valid_regions_path);
    info!("  Target regions: {}", target_regions_path);
    
    // Load existing valid regions
    let valid_regions = load_valid_regions(std::path::Path::new(valid_regions_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to load valid regions: {}", e)))?;
    
    // Load target regions  
    let target_regions = load_exclude_regions(std::path::Path::new(target_regions_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to load target regions: {}", e)))?;
    
    // Filter valid regions to those overlapping target regions
    let mut ontarget_regions = Vec::new();
    
    for region in &valid_regions {
        let chrom_key = norm_chrom(&region.chrom).to_string();
        
        if let Some(targets) = target_regions.get(&chrom_key) {
            // Check if this valid region overlaps any target region
            for (t_start, t_end) in targets {
                if region.start < *t_end && region.end > *t_start {
                    // Overlaps - include this valid region
                    ontarget_regions.push(region.clone());
                    break;
                }
            }
        }
    }
    
    info!("  {} of {} valid regions overlap targets", ontarget_regions.len(), valid_regions.len());
    
    // Write output
    let output_file = File::create(output_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to create output file: {}", e)))?;
    
    let mut writer = noodles::bgzf::io::Writer::new(output_file);
    
    for region in &ontarget_regions {
        let line = format!("{}\t{}\t{}\t{:.6}\n", region.chrom, region.start, region.end, region.gc_content);
        writer.write_all(line.as_bytes())
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to write: {}", e)))?;
    }
    
    writer.finish()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to finish BGZF: {}", e)))?;
    
    let count = ontarget_regions.len();
    info!("Generated {} on-target valid regions", count);
    
    Ok(count)
}


/// Generate valid regions BED file by excluding problematic and gap regions
/// 
/// # Arguments
/// * `reference_path` - Path to reference FASTA
/// * `exclude_regions_path` - Path to exclude regions BED (e.g., ENCODE list)
/// * `output_path` - Path to write valid_regions.bed
/// * `bin_size` - Size of each bin (default: 100000)
/// 
/// # Returns
/// * Number of valid regions generated
#[pyfunction]
#[pyo3(signature = (reference_path, exclude_regions_path, output_path, bin_size=100000))]
pub fn generate_valid_regions(
    reference_path: &str,
    exclude_regions_path: &str,
    output_path: &str,
    bin_size: u64,
) -> PyResult<usize> {
    info!("Generating valid regions from reference: {}", reference_path);
    info!("Exclude regions: {}", exclude_regions_path);
    info!("Bin size: {} bp", bin_size);
    
    // Load reference index to get chromosome lengths
    let faidx = faidx::Reader::from_path(reference_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to open reference FASTA: {}. Make sure .fai exists.", e)
        ))?;
    
    // Load exclude regions
    let exclude_regions = load_exclude_regions(Path::new(exclude_regions_path))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load exclude regions: {}", e)
        ))?;
    
    // Valid chromosomes (1-22, X, Y) - try with and without chr prefix
    let valid_chroms: Vec<String> = (1..=22)
        .map(|i| i.to_string())
        .chain(vec!["X".to_string(), "Y".to_string()])
        .collect();
    
    let mut regions: Vec<ValidRegion> = Vec::new();
    
    // Get number of sequences in the reference
    let n_seqs = faidx.n_seqs();
    info!("Reference has {} sequences", n_seqs);
    
    // Build a map of chromosome name -> length from reference
    let mut chrom_lengths: HashMap<String, u64> = HashMap::new();
    for i in 0..n_seqs {
        if let Ok(name) = faidx.seq_name(i as i32) {
            // Normalize chromosome name (remove "chr" prefix if present)
            let normalized = norm_chrom(&name).to_string();
            
            // Only include valid chromosomes
            if valid_chroms.contains(&normalized) {
                // Get sequence length
                let len = faidx.fetch_seq_len(&name);
                chrom_lengths.insert(normalized, len);
            }
        }
    }
    
    info!("Found {} valid chromosomes in reference", chrom_lengths.len());
    
    // Generate bins for each chromosome, excluding problematic regions
    for chrom in &valid_chroms {
        if let Some(&chrom_len) = chrom_lengths.get(chrom) {
            let excluded = exclude_regions.get(chrom).cloned().unwrap_or_default();
            
            // Generate bins of bin_size, excluding those that overlap with excluded regions
            let mut pos: u64 = 0;
            while pos + bin_size <= chrom_len {
                let bin_start = pos;
                let bin_end = pos + bin_size;
                
                // Check if this bin overlaps with any excluded region
                let overlaps_excluded = excluded.iter().any(|(excl_start, excl_end)| {
                    // Overlap if: bin_start < excl_end AND bin_end > excl_start
                    bin_start < *excl_end && bin_end > *excl_start
                });
                
                if !overlaps_excluded {
                    regions.push(ValidRegion {
                        chrom: chrom.clone(),
                        start: bin_start,
                        end: bin_end,
                        gc_content: 0.0,  // Will be computed later
                    });
                }
                
                pos += bin_size;
            }
        }
    }
    
    info!("Generated {} valid regions (excluded {} problematic bins)", 
          regions.len(), 
          chrom_lengths.values().map(|l| l / bin_size).sum::<u64>() as usize - regions.len());
    
    // Write output BED.gz (BGZF compressed)
    use noodles::bgzf;
    use std::io::Write;
    
    let output_file = File::create(output_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    let mut writer = bgzf::io::Writer::new(output_file);
    
    // Write header
    writeln!(writer, "# Valid regions for GC correction estimation")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    writeln!(writer, "# Excludes: problematic regions, assembly gaps, low mappability")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    writeln!(writer, "#chrom\tstart\tend")
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    
    // Write regions
    for region in &regions {
        writeln!(writer, "{}\t{}\t{}", region.chrom, region.start, region.end)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    }
    
    // Finish BGZF stream
    writer.finish()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to finish BGZF stream: {}", e)
        ))?;
    
    let region_count = regions.len();
    
    info!("Generated {} valid regions", region_count);
    Ok(region_count)
}

/// Load exclude regions from BED file (supports plain and gzipped)
fn load_exclude_regions(path: &Path) -> Result<HashMap<String, Vec<(u64, u64)>>> {
    // Use bed::get_reader for transparent BGZF/gzip handling
    let reader = crate::bed::get_reader(path)
        .with_context(|| format!("Failed to open exclude regions file: {:?}", path))?;
    
    let mut exclude_regions: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        
        let chrom = norm_chrom(fields[0]).to_string();
        let start: u64 = fields[1].parse().unwrap_or(0);
        let end: u64 = fields[2].parse().unwrap_or(0);
        
        exclude_regions.entry(chrom).or_default().push((start, end));
    }
    
    // Sort intervals by start position
    for intervals in exclude_regions.values_mut() {
        intervals.sort_by_key(|i| i.0);
    }
    
    info!("Loaded exclude regions with {} chromosomes", exclude_regions.len());
    Ok(exclude_regions)
}

/// Generate reference genome GC expected counts
/// 
/// # Algorithm
/// For each length bin (60-80, 80-100, ..., 380-400):
///   For each GC percent (0-100):
///     Count how many theoretical fragments of that length have that GC
///     
/// This represents the "expected" distribution under uniform sampling.
/// 
/// # Arguments
/// * `reference_path` - Path to reference FASTA
/// * `valid_regions_path` - Path to valid_regions.bed
/// * `output_path` - Path to write ref_genome_GC.parquet
/// 
/// # Returns
/// * Total fragments counted
#[pyfunction]
#[pyo3(signature = (reference_path, valid_regions_path, output_path))]
pub fn generate_ref_genome_gc(
    reference_path: &str,
    valid_regions_path: &str,
    output_path: &str,
) -> PyResult<u64> {
    info!("Generating ref_genome_GC from reference: {}", reference_path);
    info!("Using valid regions: {}", valid_regions_path);
    
    // Load reference
    let faidx = faidx::Reader::from_path(reference_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to open reference FASTA: {}", e)
        ))?;
    let contig_names = ContigNames::from_faidx(&faidx);
    
    // Load valid regions
    let regions = load_valid_regions(Path::new(valid_regions_path))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to load valid regions: {}", e)
        ))?;
    
    // Initialize counts matrix [length_bin][gc_percent]
    // Uses NUM_LENGTH_BINS (188 bins of 5bp spanning 60..999)
    let mut counts: Vec<Vec<u64>> = vec![vec![0; 101]; NUM_LENGTH_BINS];
    let mut total_fragments: u64 = 0;
    
    info!("Processing {} regions...", regions.len());
    
    // Progress tracking
    let total_regions = regions.len();
    let mut processed_regions = 0;
    
    // For each valid region, scan for GC content at each position for each length
    for region in &regions {
        // Valid regions name contigs without `chr`; fetch under the reference's own name
        // (contigs.rs). Every region must resolve: the expected-GC table is only right if each
        // contributes, and a region the reference lacks means the valid-regions file was built
        // from another reference. (The try-both-spellings fallback this replaced never fell
        // through: rust-htslib's `fetch_seq` returns a `Vec` over a null pointer, not `Err`.)
        let fasta_name = contig_names.resolve(&region.chrom, 0).name().map(str::to_string).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "valid region on contig {} is not in the reference FASTA under either spelling; \
                 the valid-regions file was built from another reference", region.chrom))
        })?;
        let seq = faidx.fetch_seq(&fasta_name, region.start as usize, region.end as usize - 1)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!(
                "cannot fetch {}:{}-{} from the reference: {}", fasta_name, region.start, region.end, e)))?;
        
        let seq_len = seq.len();
        
        // For each length bin, sample positions and count GC
        for (bin_idx, &(min_len, max_len)) in LENGTH_BINS.iter().enumerate() {
            // Use the midpoint of each length bin for sampling
            let frag_len = ((min_len + max_len) / 2) as usize;
            
            // Sample every 100bp to speed up (still statistically representative)
            let step = 100;
            
            let mut pos = 0;
            while pos + frag_len <= seq_len {
                // Get the subsequence for this theoretical fragment
                let frag_seq = &seq[pos..pos + frag_len];
                
                // Count GC content
                let gc_count = frag_seq.iter()
                    .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
                    .count();
                
                // Skip if too many Ns (uncertain bases)
                let n_count = frag_seq.iter()
                    .filter(|&&b| b == b'N' || b == b'n')
                    .count();
                
                if n_count * 10 > frag_len {
                    // More than 10% Ns, skip this position
                    pos += step;
                    continue;
                }
                
                // Calculate GC percent (0-100)
                let valid_bases = frag_len - n_count;
                let gc_percent = if valid_bases > 0 {
                    ((gc_count as f64 / valid_bases as f64) * 100.0).round() as usize
                } else {
                    50  // Default to 50% if no valid bases
                };
                
                let gc_idx = gc_percent.min(100);
                
                // Increment count
                counts[bin_idx][gc_idx] += 1;
                total_fragments += 1;
                
                pos += step;
            }
        }
        
        processed_regions += 1;
        
        // Log progress every 1000 regions
        if processed_regions % 1000 == 0 {
            info!("Processed {}/{} regions ({:.1}%), {} fragments so far", 
                  processed_regions, total_regions,
                  (processed_regions as f64 / total_regions as f64) * 100.0,
                  total_fragments);
        }
    }
    
    info!("Completed GC counting: {} total fragments across {} regions", 
          total_fragments, processed_regions);
    
    // Save as Parquet file (consistent with PON format)
    // Schema: length_bin_min (u32), length_bin_max (u32), gc_percent (u8), expected_count (u64)
    use arrow::array::{UInt32Array, UInt8Array, UInt64Array};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;
    
    // Build arrays from counts
    let mut length_bin_min_vec: Vec<u32> = Vec::new();
    let mut length_bin_max_vec: Vec<u32> = Vec::new();
    let mut gc_percent_vec: Vec<u8> = Vec::new();
    let mut expected_count_vec: Vec<u64> = Vec::new();
    
    for (bin_idx, &(min_len, max_len)) in LENGTH_BINS.iter().enumerate() {
        for gc_percent in 0..=100u8 {
            let count = counts[bin_idx][gc_percent as usize];
            if count > 0 {  // Only write non-zero counts to save space
                length_bin_min_vec.push(min_len);
                length_bin_max_vec.push(max_len);
                gc_percent_vec.push(gc_percent);
                expected_count_vec.push(count);
            }
        }
    }
    
    // Create arrow arrays
    let length_bin_min_array = UInt32Array::from(length_bin_min_vec);
    let length_bin_max_array = UInt32Array::from(length_bin_max_vec);
    let gc_percent_array = UInt8Array::from(gc_percent_vec);
    let expected_count_array = UInt64Array::from(expected_count_vec);
    
    // Create schema
    let schema = Arc::new(Schema::new(vec![
        Field::new("length_bin_min", DataType::UInt32, false),
        Field::new("length_bin_max", DataType::UInt32, false),
        Field::new("gc_percent", DataType::UInt8, false),
        Field::new("expected_count", DataType::UInt64, false),
    ]));
    
    // Create record batch
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(length_bin_min_array),
            Arc::new(length_bin_max_array),
            Arc::new(gc_percent_array),
            Arc::new(expected_count_array),
        ],
    ).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
        format!("Failed to create record batch: {}", e)
    ))?;
    
    // Write to Parquet file
    let output_file = std::fs::File::create(output_path)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create output file: {}", e)
        ))?;
    
    let mut writer = ArrowWriter::try_new(output_file, schema, None)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to create Parquet writer: {}", e)
        ))?;
    
    writer.write(&batch).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
        format!("Failed to write Parquet batch: {}", e)
    ))?;
    
    writer.close().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
        format!("Failed to close Parquet file: {}", e)
    ))?;
    
    info!("Generated ref_genome_GC with {} total fragments", total_fragments);
    Ok(total_fragments)
}

// =================================================================================================
// Expected counts over a territory (GC redesign item 1). Additive: nothing above changes.
// =================================================================================================

use crate::gc_model::{base_class, sha256_hex, GcGrid, GC_DEFINITION_VERSION};
use serde::{Deserialize, Serialize};

/// One interval of a [`Territory`], under the reference FASTA's own contig name.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct TerritoryInterval {
    pub contig: String,
    pub contig_len: u64,
    pub start: u64,
    pub end: u64,
}

/// The genomic intervals a GC model is fitted over, resolved against the reference FASTA through
/// [`ContigNames`] (the one contig-naming rule). Sorted, and never overlapping: a window counts
/// only if it lies wholly inside **one** interval, and adjacent intervals are not merged.
#[derive(Debug, Clone, PartialEq)]
pub struct Territory {
    intervals: Vec<TerritoryInterval>,
}

impl Territory {
    /// Resolve BED-style records `(contig, start, end)` (0-based, half-open) against the FASTA.
    /// Refuses a contig the FASTA lacks under either spelling, an empty interval or one past its
    /// contig's end, and overlapping intervals.
    pub fn resolve(records: &[(String, u64, u64)], names: &ContigNames) -> Result<Self> {
        let mut intervals = Vec::with_capacity(records.len());
        for (name, start, end) in records {
            let contig = names.resolve(name, 0).name().map(str::to_string).ok_or_else(|| anyhow::anyhow!(
                "territory contig {name} is not in the reference FASTA under either spelling; the \
                 territory was built for another reference"))?;
            let contig_len = names.length(&contig).expect("a resolved contig has a length");
            if start >= end || *end > contig_len {
                return Err(anyhow::anyhow!(
                    "territory interval {name}:{start}-{end} is empty or runs past the end of {contig} \
                     ({contig_len} bp)"));
            }
            intervals.push(TerritoryInterval { contig, contig_len, start: *start, end: *end });
        }
        intervals.sort();
        for pair in intervals.windows(2) {
            if pair[0].contig == pair[1].contig && pair[1].start < pair[0].end {
                return Err(anyhow::anyhow!(
                    "territory intervals {}:{}-{} and {}:{}-{} overlap; windows in the overlap would \
                     be counted twice", pair[0].contig, pair[0].start, pair[0].end,
                    pair[1].contig, pair[1].start, pair[1].end));
            }
        }
        Ok(Self { intervals })
    }

    pub fn intervals(&self) -> &[TerritoryInterval] {
        &self.intervals
    }

    /// Whether `[start, end)` on `contig` (the FASTA's own name, as [`ContigNames`] resolves it)
    /// lies wholly inside **one** interval: the rule [`count_windows`] applies to the expected
    /// counts, so a tally built with it counts the same fragments. Adjacent intervals are not
    /// merged.
    pub fn contains(&self, contig: &str, start: u64, end: u64) -> bool {
        // The last interval on `contig` starting at or before `start` is the only candidate:
        // intervals are sorted and never overlap.
        let i = self.intervals.partition_point(|iv| (iv.contig.as_str(), iv.start) <= (contig, start));
        i > 0 && {
            let iv = &self.intervals[i - 1];
            iv.contig == contig && start < end && end <= iv.end
        }
    }

    /// SHA-256 over every interval's FASTA contig name, contig length, start and end, in order.
    /// The contig length is in it so that the same coordinates on another assembly differ.
    pub fn digest(&self) -> String {
        let mut text = String::new();
        for iv in &self.intervals {
            text.push_str(&format!("{}\t{}\t{}\t{}\n", iv.contig, iv.contig_len, iv.start, iv.end));
        }
        sha256_hex(text.as_bytes())
    }
}

/// Add every window of `seq` to `counts`: for each length of the grid (outer) and each start
/// position (inner), the window's key over its trimmed window ([`crate::gc_model::gc_key`]),
/// from prefix sums of GC and of bases that are not A/C/G/T. A window containing any such base in
/// its trimmed window is not counted.
pub fn count_windows(seq: &[u8], counts: &mut GcGrid<u64>) {
    // Prefix sums: gc[i] and bad[i] count G/C and non-ACGT bases in seq[..i].
    let n = seq.len();
    let mut gc = Vec::with_capacity(n + 1);
    let mut bad = Vec::with_capacity(n + 1);
    gc.push(0u32);
    bad.push(0u32);
    for &b in seq {
        let (counts_base, is_gc) = base_class(b);
        gc.push(gc[gc.len() - 1] + is_gc as u32);
        bad.push(bad[bad.len() - 1] + !counts_base as u32);
    }
    let k = counts.trim_k() as usize;
    for len in counts.min_len()..=counts.max_len() {
        let l = len as usize;
        if l > n {
            break;
        }
        let row = counts.row_mut(len);
        for p in 0..=n - l {
            let (a, b) = (p + k, p + l - k);
            if bad[b] == bad[a] {
                row[(gc[b] - gc[a]) as usize] += 1;
            }
        }
    }
}

/// Provenance of an expected table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExpectedProvenance {
    pub genome: String,
    pub territory_sha256: String,
    pub trim_k: u16,
    pub gc_definition_version: u32,
    pub krewlyzer_version: String,
}

/// Expected counts: for each exact (length, GC count), the number of windows of that length in the
/// territory with that key.
#[derive(Debug, Clone, PartialEq)]
pub struct ExpectedTable {
    counts: GcGrid<u64>,
    provenance: ExpectedProvenance,
}

impl ExpectedTable {
    /// Count every window of every length in `min_len..=max_len` over the territory, fetching each
    /// interval's sequence with `fetch`.
    pub fn count<F>(territory: &Territory, genome: &str, min_len: u16, max_len: u16, trim_k: u16,
                    mut fetch: F) -> Result<Self>
    where
        F: FnMut(&TerritoryInterval) -> Result<Vec<u8>>,
    {
        let mut counts = GcGrid::new(min_len, max_len, trim_k, 0u64)?;
        let started = std::time::Instant::now();
        for iv in territory.intervals() {
            let seq = fetch(iv)?;
            if seq.len() as u64 != iv.end - iv.start {
                return Err(anyhow::anyhow!(
                    "fetched {} bp for {}:{}-{}, expected {}", seq.len(), iv.contig, iv.start, iv.end,
                    iv.end - iv.start));
            }
            count_windows(&seq, &mut counts);
        }
        info!("gc_reference: counted {} windows over {} intervals, lengths {}..={} at trim {}, in {:.1} s",
              counts.cells().iter().sum::<u64>(), territory.intervals().len(), min_len, max_len, trim_k,
              started.elapsed().as_secs_f64());
        Ok(Self {
            counts,
            provenance: ExpectedProvenance {
                genome: genome.to_string(),
                territory_sha256: territory.digest(),
                trim_k,
                gc_definition_version: GC_DEFINITION_VERSION,
                krewlyzer_version: env!("CARGO_PKG_VERSION").to_string(),
            },
        })
    }

    pub fn counts(&self) -> &GcGrid<u64> {
        &self.counts
    }

    /// A table from counts made by hand, for tests of the fit (which needs no sequence).
    #[cfg(test)]
    pub(crate) fn from_counts(counts: GcGrid<u64>, genome: &str, territory_sha256: &str) -> Self {
        let provenance = ExpectedProvenance {
            genome: genome.to_string(),
            territory_sha256: territory_sha256.to_string(),
            trim_k: counts.trim_k(),
            gc_definition_version: GC_DEFINITION_VERSION,
            krewlyzer_version: env!("CARGO_PKG_VERSION").to_string(),
        };
        Self { counts, provenance }
    }

    pub fn provenance(&self) -> &ExpectedProvenance {
        &self.provenance
    }

    /// SHA-256 of the table: its genome, territory digest, trim, GC definition version, length
    /// range and every count in index order (not the krewlyzer version that built it). A support
    /// records it, so a fit can refuse a support built from another table.
    pub fn digest(&self) -> String {
        let p = &self.provenance;
        let mut bytes = format!(
            "genome={}\nterritory_sha256={}\ntrim_k={}\ngc_definition_version={}\nlengths={}..={}\n",
            p.genome, p.territory_sha256, p.trim_k, p.gc_definition_version,
            self.counts.min_len(), self.counts.max_len()).into_bytes();
        for c in self.counts.cells() {
            bytes.extend_from_slice(&c.to_le_bytes());
        }
        sha256_hex(&bytes)
    }

    /// One row per cell: `length`, `gc` (UInt16), `expected` (UInt64); the provenance and the
    /// length range in the file's key-value metadata.
    pub fn write_parquet<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use arrow::array::{ArrayRef, UInt16Array, UInt64Array};
        use arrow::datatypes::{DataType, Field, Schema};
        use std::sync::Arc;
        let keys: Vec<_> = self.counts.keys().collect();
        let schema = Arc::new(Schema::new(vec![
            Field::new("length", DataType::UInt16, false),
            Field::new("gc", DataType::UInt16, false),
            Field::new("expected", DataType::UInt64, false),
        ]));
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.length))),
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.gc))),
            Arc::new(UInt64Array::from(self.counts.cells().to_vec())),
        ];
        let range = (self.counts.min_len(), self.counts.max_len());
        let meta = vec![
            (EXPECTED_META_PROVENANCE, serde_json::to_string(&self.provenance)?),
            (EXPECTED_META_LENGTHS, serde_json::to_string(&range)?),
        ];
        crate::gc_model::write_parquet_with_metadata(path.as_ref(), schema, arrays, meta)
    }

    /// Load a table written by [`Self::write_parquet`]; refuses one with no provenance, another
    /// GC definition version, or a missing, repeated or foreign cell.
    pub fn load_parquet<P: AsRef<Path>>(path: P) -> Result<Self> {
        use arrow::array::{UInt16Array, UInt64Array};
        use crate::gc_model::{claim_cell, column, parse_meta, read_parquet_with_metadata, require_every_cell};
        let path = path.as_ref();
        let (meta, batches) = read_parquet_with_metadata(path)?;
        let provenance: ExpectedProvenance = parse_meta(&meta, EXPECTED_META_PROVENANCE, path)?;
        if provenance.gc_definition_version != GC_DEFINITION_VERSION {
            return Err(anyhow::anyhow!(
                "expected table {:?} was written under GC definition version {}; this krewlyzer \
                 reads version {} only. Rebuild it; it is not reinterpreted.",
                path, provenance.gc_definition_version, GC_DEFINITION_VERSION));
        }
        let (min_len, max_len): (u16, u16) = parse_meta(&meta, EXPECTED_META_LENGTHS, path)?;
        let mut counts = GcGrid::new(min_len, max_len, provenance.trim_k, 0u64)?;
        let mut seen = vec![false; counts.len()];
        for batch in &batches {
            let length = column::<UInt16Array>(batch, "length", path)?;
            let gc = column::<UInt16Array>(batch, "gc", path)?;
            let expected = column::<UInt64Array>(batch, "expected", path)?;
            for r in 0..batch.num_rows() {
                let i = claim_cell(&counts, &mut seen, length.value(r), gc.value(r), path)?;
                counts.cells_mut()[i] = expected.value(r);
            }
        }
        require_every_cell(&seen, path)?;
        Ok(Self { counts, provenance })
    }
}

/// [`ExpectedTable::count`] over a reference FASTA (with `.fai`): `records` are BED-style
/// intervals, resolved through [`ContigNames`].
pub fn expected_counts(fasta: &Path, genome: &str, records: &[(String, u64, u64)],
                       min_len: u16, max_len: u16, trim_k: u16) -> Result<ExpectedTable> {
    let fa = faidx::Reader::from_path(fasta)
        .map_err(|e| anyhow::anyhow!("cannot open reference FASTA {:?}: {}", fasta, e))?;
    let territory = Territory::resolve(records, &ContigNames::from_faidx(&fa))?;
    // Fetch only under the FASTA's own names (contigs.rs): htslib answers an unknown name with a
    // Vec over a null pointer rather than an error.
    ExpectedTable::count(&territory, genome, min_len, max_len, trim_k, |iv| {
        fa.fetch_seq(&iv.contig, iv.start as usize, iv.end as usize - 1)
            .map(|s| s.to_vec())
            .map_err(|e| anyhow::anyhow!("cannot fetch {}:{}-{}: {}", iv.contig, iv.start, iv.end, e))
    })
}

const EXPECTED_META_PROVENANCE: &str = "krewlyzer.gc_expected.provenance";
const EXPECTED_META_LENGTHS: &str = "krewlyzer.gc_expected.lengths";

#[cfg(test)]
mod expected_tests {
    use super::*;
    use crate::gc_model::gc_key;

    /// A deterministic pseudo-random sequence (SplitMix64) with lowercase bases and a few Ns.
    fn synthetic(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n).map(|_| {
            x = x.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = x;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^= z >> 31;
            match z % 1000 { 0 => b'N', 1 => b'n', 2 => b'R', 3..=150 => b'g', r => b"ACGTacgt"[(r % 8) as usize] }
        }).collect()
    }

    fn names() -> ContigNames {
        ContigNames::new([("1".to_string(), 5_000), ("2".to_string(), 3_000)])
    }

    fn genome() -> HashMap<String, Vec<u8>> {
        HashMap::from([("1".to_string(), synthetic(5_000, 1)), ("2".to_string(), synthetic(3_000, 2))])
    }

    fn table(records: &[(String, u64, u64)], min_len: u16, max_len: u16, k: u16) -> ExpectedTable {
        let g = genome();
        let t = Territory::resolve(records, &names()).unwrap();
        ExpectedTable::count(&t, "synthetic", min_len, max_len, k, |iv| {
            Ok(g[&iv.contig][iv.start as usize..iv.end as usize].to_vec())
        }).unwrap()
    }

    fn recs(r: &[(&str, u64, u64)]) -> Vec<(String, u64, u64)> {
        r.iter().map(|(c, s, e)| (c.to_string(), *s, *e)).collect()
    }

    #[test]
    fn counts_equal_a_brute_force_count_for_every_cell() {
        let r = recs(&[("1", 100, 1_400), ("chr2", 0, 900), ("1", 2_000, 2_150)]);
        let (lo, hi, k) = (30, 160, 5);
        let got = table(&r, lo, hi, k);
        let mut want = GcGrid::new(lo, hi, k, 0u64).unwrap();
        let g = genome();
        for (c, s, e) in &r {
            let seq = &g[norm_chrom(c)][*s as usize..*e as usize];
            for len in lo..=hi {
                for p in 0..seq.len().saturating_sub(len as usize - 1) {
                    if let Some(key) = gc_key(&seq[p..p + len as usize], k) {
                        let i = want.index(key.length, key.gc).unwrap();
                        want.cells_mut()[i] += 1;
                    }
                }
            }
        }
        assert!(want.cells().iter().sum::<u64>() > 100_000, "the case must exercise many windows");
        assert_eq!(got.counts(), &want);
    }

    #[test]
    fn a_window_with_n_in_its_trimmed_window_is_never_counted() {
        // 100 bp of A with one N at 50, windows of 40 bp at k = 0: 61 windows, 40 of them cover it.
        let mut seq = vec![b'A'; 100];
        seq[50] = b'N';
        let mut g = GcGrid::new(40, 40, 0, 0u64).unwrap();
        count_windows(&seq, &mut g);
        assert_eq!(g.cells().iter().sum::<u64>(), 21);
        // At k = 5 a window whose N falls in a trimmed-off end still counts: 61 - 30 = 31.
        let mut g = GcGrid::new(40, 40, 5, 0u64).unwrap();
        count_windows(&seq, &mut g);
        assert_eq!(g.cells().iter().sum::<u64>(), 31);
    }

    #[test]
    fn two_territories_give_two_tables_and_two_checksums() {
        let a = table(&recs(&[("1", 0, 2_000)]), 40, 120, 2);
        let b = table(&recs(&[("1", 2_000, 4_000)]), 40, 120, 2);
        assert_ne!(a.counts(), b.counts());
        assert_ne!(a.provenance().territory_sha256, b.provenance().territory_sha256);
    }

    #[test]
    fn adjacent_intervals_are_not_merged() {
        let whole = table(&recs(&[("1", 0, 200)]), 40, 40, 0);
        let split = table(&recs(&[("1", 0, 100), ("1", 100, 200)]), 40, 40, 0);
        let n = |t: &ExpectedTable| t.counts().cells().iter().sum::<u64>();
        assert!(n(&split) < n(&whole), "split {} vs whole {}", n(&split), n(&whole));
    }

    #[test]
    fn contains_agrees_with_count_windows_window_by_window() {
        // A 300 bp contig; intervals with a gap, two adjacent ones, and one at the contig's end.
        let n = ContigNames::new([("1".to_string(), 300)]);
        let t = Territory::resolve(&recs(&[("1", 10, 60), ("1", 80, 130), ("1", 130, 170), ("1", 250, 300)]), &n)
            .unwrap();
        for len in [20u16, 25, 50, 51] {
            for start in 0..=(300 - len as u64) {
                let end = start + len as u64;
                // A genome of N except this window: the builder counts exactly this window, or none.
                let mut seq = vec![b'N'; 300];
                seq[start as usize..end as usize].fill(b'A');
                let table = ExpectedTable::count(&t, "g", len, len, 0, |iv| {
                    Ok(seq[iv.start as usize..iv.end as usize].to_vec())
                }).unwrap();
                let counted = table.counts().cells().iter().sum::<u64>();
                assert!(counted <= 1);
                assert_eq!(t.contains("1", start, end), counted == 1, "length {len} at {start}");
            }
        }
        assert!(!t.contains("2", 10, 30), "another contig");
        assert!(!t.contains("chr1", 10, 30), "contains takes the FASTA's own name");
    }

    #[test]
    fn territory_contigs_resolve_through_contig_names() {
        let n = names();
        let ucsc = Territory::resolve(&recs(&[("chr1", 10, 20)]), &n).unwrap();
        let plain = Territory::resolve(&recs(&[("1", 10, 20)]), &n).unwrap();
        assert_eq!(ucsc.intervals()[0].contig, "1", "fetch under the FASTA's own name");
        assert_eq!(ucsc.digest(), plain.digest(), "the spelling is not part of the territory");
        for bad in [recs(&[("3", 0, 10)]), recs(&[("1", 10, 10)]), recs(&[("1", 4_990, 5_001)]),
                    recs(&[("1", 0, 100), ("chr1", 50, 150)])] {
            assert!(Territory::resolve(&bad, &n).is_err(), "{bad:?} must be refused");
        }
        // Another assembly with the same coordinates: another digest.
        let other = ContigNames::new([("1".to_string(), 6_000)]);
        assert_ne!(Territory::resolve(&recs(&[("1", 10, 20)]), &other).unwrap().digest(), plain.digest());
        // Order in the BED does not matter.
        let ab = Territory::resolve(&recs(&[("1", 0, 10), ("2", 0, 10)]), &n).unwrap();
        let ba = Territory::resolve(&recs(&[("2", 0, 10), ("1", 0, 10)]), &n).unwrap();
        assert_eq!(ab.digest(), ba.digest());
    }

    #[test]
    fn the_fasta_builder_matches_the_in_memory_count() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let g = genome();
        let mut text = String::new();
        for c in ["1", "2"] {
            text.push_str(&format!(">{c}\n"));
            for line in g[c].chunks(60) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
        }
        std::fs::write(&fa, text).unwrap();
        let r = recs(&[("chr1", 100, 1_400), ("2", 0, 900)]);
        let from_file = expected_counts(&fa, "synthetic", &r, 30, 90, 5).unwrap();
        assert_eq!(from_file, table(&r, 30, 90, 5));
        assert!(expected_counts(&fa, "synthetic", &recs(&[("chr9", 0, 10)]), 30, 90, 5).is_err());
    }

    #[test]
    fn the_expected_table_round_trips_and_refuses_another_version() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("e.parquet");
        let t = table(&recs(&[("1", 0, 1_000)]), 30, 60, 2);
        t.write_parquet(&p).unwrap();
        assert_eq!(ExpectedTable::load_parquet(&p).unwrap(), t);
        let mut v = t.clone();
        v.provenance.gc_definition_version += 1;
        v.write_parquet(&p).unwrap();
        let err = format!("{:#}", ExpectedTable::load_parquet(&p).unwrap_err());
        assert!(err.contains("definition version"), "{err}");
        // The same rows with no metadata: refused.
        use arrow::array::{ArrayRef, UInt16Array, UInt64Array};
        use arrow::datatypes::{DataType, Field, Schema};
        use std::sync::Arc;
        let keys: Vec<_> = t.counts().keys().collect();
        let schema = Arc::new(Schema::new(vec![
            Field::new("length", DataType::UInt16, false),
            Field::new("gc", DataType::UInt16, false),
            Field::new("expected", DataType::UInt64, false),
        ]));
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.length))),
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.gc))),
            Arc::new(UInt64Array::from(t.counts().cells().to_vec())),
        ];
        crate::output_utils::write_parquet_batch(&p, schema, arrays).unwrap();
        let err = format!("{:#}", ExpectedTable::load_parquet(&p).unwrap_err());
        assert!(err.contains(EXPECTED_META_PROVENANCE), "{err}");
    }

    #[test]
    fn a_null_expected_count_is_refused() {
        use arrow::array::{ArrayRef, UInt16Array, UInt64Array};
        use arrow::datatypes::{DataType, Field, Schema};
        use std::sync::Arc;
        let dir = tempfile::tempdir().unwrap();
        let good = dir.path().join("e.parquet");
        let t = table(&recs(&[("1", 0, 1_000)]), 30, 60, 2);
        t.write_parquet(&good).unwrap();
        let (meta, _) = crate::gc_model::read_parquet_with_metadata(&good).unwrap();
        let keys: Vec<_> = t.counts().keys().collect();
        let schema = Arc::new(Schema::new(vec![
            Field::new("length", DataType::UInt16, false),
            Field::new("gc", DataType::UInt16, false),
            Field::new("expected", DataType::UInt64, true),
        ]));
        let expected: Vec<Option<u64>> = t.counts().cells().iter().enumerate()
            .map(|(i, &c)| if i == 5 { None } else { Some(c) }).collect();
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.length))),
            Arc::new(UInt16Array::from_iter_values(keys.iter().map(|k| k.gc))),
            Arc::new(UInt64Array::from(expected)),
        ];
        let bad = dir.path().join("null.parquet");
        let meta: Vec<(&str, String)> = meta.iter().map(|(k, v)| (k.as_str(), v.clone())).collect();
        crate::gc_model::write_parquet_with_metadata(&bad, schema, arrays, meta).unwrap();
        let err = format!("{:#}", ExpectedTable::load_parquet(&bad).unwrap_err());
        assert!(err.contains("`expected`") && err.contains("null"), "{err}");
    }

    /// Synthetic-scale timing for the expected-count builder (`cargo test --release -- --ignored`).
    /// The hg19 measurement needs the FASTA, so it is a cluster job.
    #[test]
    #[ignore]
    fn timing_on_a_synthetic_megabase() {
        let seq = synthetic(1_000_000, 7);
        let (lo, hi) = GcGrid::<u64>::lengths_for(10);
        let mut g = GcGrid::new(lo, hi, 10, 0u64).unwrap();
        let t0 = std::time::Instant::now();
        count_windows(&seq, &mut g);
        let s = t0.elapsed().as_secs_f64();
        let windows: u64 = g.cells().iter().sum();
        eprintln!("count_windows: 1 Mbp x {} lengths in {:.2} s ({} windows counted, {:.2} ns/window \
                   position); extrapolated {:.2} core-hours per Gbp",
                  hi - lo + 1, s, windows, s * 1e9 / (1e6 * (hi - lo + 1) as f64), s * 1000.0 / 3600.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_length_bin_index() {
        // 5bp bins spanning 60..999 -> NUM_LENGTH_BINS (188) bins:
        // 60-64 = bin 0, 65-69 = bin 1, ..., 995-999 = bin 187
        assert_eq!(get_length_bin_index(60), Some(0));
        assert_eq!(get_length_bin_index(64), Some(0));  // Still bin 0
        assert_eq!(get_length_bin_index(65), Some(1));  // Now bin 1
        assert_eq!(get_length_bin_index(79), Some(3));  // (79-60)/5 = 3
        assert_eq!(get_length_bin_index(80), Some(4));  // (80-60)/5 = 4
        assert_eq!(get_length_bin_index(166), Some(21)); // Nucleosome peak ~166bp
        assert_eq!(get_length_bin_index(400), Some(68)); // Multi-nucleosomal range
        assert_eq!(get_length_bin_index(999), Some(NUM_LENGTH_BINS - 1)); // Last bin
        assert_eq!(get_length_bin_index(1000), None);    // Out of range (exclusive)
        assert_eq!(get_length_bin_index(59), None);      // Out of range
    }

}
