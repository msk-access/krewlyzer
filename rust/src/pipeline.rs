//! Unified Feature Pipeline Orchestration
//!
//! Single-pass feature extraction from BED.gz fragments.
//! Coordinates multiple consumers (FSC, WPS, FSD, OCF) in a single read.
//!
//! Entry point: `run_unified_pipeline()` - exposed to Python via PyO3.
//!
//! Features:
//! - GC correction (compute or load factors)
//! - On/off-target splitting for panel data
//! - BaitMask for WPS edge detection

use pyo3::prelude::*;
use std::path::PathBuf;
use anyhow::{Result, Context};
use log::{info, warn};
use crate::bed::{Fragment, ChromosomeMap};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use crate::fsc::FscConsumer;
use crate::wps::{WpsConsumer, WpsBackgroundConsumer};
use crate::fsd::FsdConsumer;
use crate::ocf::OcfConsumer;

/// A consumer that delegates to multiple other consumers
#[derive(Clone)]
pub struct MultiConsumer {
    fsc: Option<FscConsumer>,
    wps: Option<WpsConsumer>,
    wps_background: Option<WpsBackgroundConsumer>,
    fsd: Option<FsdConsumer>,
    ocf: Option<OcfConsumer>,
}

impl MultiConsumer {
    pub fn new(
        fsc: Option<FscConsumer>,
        wps: Option<WpsConsumer>,
        wps_background: Option<WpsBackgroundConsumer>,
        fsd: Option<FsdConsumer>,
        ocf: Option<OcfConsumer>,
    ) -> Self {
        Self { fsc, wps, wps_background, fsd, ocf }
    }
    
    pub fn write_outputs(
        &self,
        fsc_path: Option<PathBuf>,
        wps_parquet_path: Option<PathBuf>,  // Foreground Parquet (TSS/CTCF)
        wps_background_path: Option<PathBuf>,  // Background Parquet (Alu stacking)
        fsd_path: Option<PathBuf>,
        ocf_dir: Option<PathBuf>,
        _empty: bool, // for WPS (unused with Parquet)
    ) -> Result<()> {
        if let Some(c) = &self.fsc {
            if let Some(p) = fsc_path {
                c.write_output(&p).context("Writing FSC output")?;
            }
        }
        
        if let Some(c) = &self.wps {
            // Write foreground Parquet (ML-ready vectors)
            if let Some(p) = wps_parquet_path {
                if let Err(e) = c.write_parquet(&p, None) {
                    // Create detailed error message that will survive PyO3 conversion
                    let full_error = format!("WPS Parquet write failed: {:#}", e);
                    return Err(anyhow::anyhow!(full_error));
                }
            }
        }
        
        if let Some(c) = &self.wps_background {
            // Write background Parquet (Alu stacking)
            if let Some(p) = wps_background_path {
                c.write_parquet(&p).context("Writing WPS Background Parquet output")?;
            }
        }
        
        if let Some(c) = &self.fsd {
            if let Some(p) = fsd_path {
                c.write_output(&p).context("Writing FSD output")?;
            }
        }
        
        if let Some(c) = &self.ocf {
            if let Some(p) = ocf_dir {
                c.write_output(&p).context("Writing OCF output")?;
            }
        }
        
        Ok(())
    }
}

impl FragmentConsumer for MultiConsumer {
    fn name(&self) -> &str {
        "MultiConsumer"
    }

    fn consume(&mut self, fragment: &Fragment) {
        if let Some(c) = &mut self.fsc { c.consume(fragment); }
        if let Some(c) = &mut self.wps { c.consume(fragment); }
        if let Some(c) = &mut self.wps_background { c.consume(fragment); }
        if let Some(c) = &mut self.fsd { c.consume(fragment); }
        if let Some(c) = &mut self.ocf { c.consume(fragment); }
    }

    fn merge(&mut self, other: Self) {
        if let (Some(a), Some(b)) = (&mut self.fsc, other.fsc) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.wps, other.wps) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.wps_background, other.wps_background) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.fsd, other.fsd) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.ocf, other.ocf) { a.merge(b); }
    }
}

#[pyfunction]
#[pyo3(signature = (
    bed_path,
    // GC Correction (Optional)
    gc_ref_path=None, valid_regions_path=None, correction_out_path=None,
    // Pre-computed factors (skip computation if provided)
    correction_input_path=None,
    // FSC Args
    fsc_bins=None, fsc_output=None,
    // WPS Args (Foreground Parquet)
    wps_regions=None, wps_output=None,
    // WPS Background Args (Alu stacking)
    wps_background_regions=None, wps_background_output=None,
    wps_empty=false,
    // FSD Args
    fsd_arms=None, fsd_output=None,
    // OCF Args
    ocf_regions=None, ocf_output=None,
    // Target regions for on/off-target split (panel data)
    target_regions_path=None,
    // Bait edge padding (adaptive safety applies)
    bait_padding=50,
    // Progress control
    silent=false
))]
pub fn run_unified_pipeline(
    _py: Python,
    bed_path: PathBuf,
    // GC Correction
    gc_ref_path: Option<PathBuf>,
    valid_regions_path: Option<PathBuf>,
    correction_out_path: Option<PathBuf>,
    // Pre-computed factors (load instead of compute if provided)
    correction_input_path: Option<PathBuf>,
    // FSC
    fsc_bins: Option<PathBuf>, fsc_output: Option<PathBuf>,
    // WPS Foreground (TSS/CTCF Parquet)
    wps_regions: Option<PathBuf>, wps_output: Option<PathBuf>,
    // WPS Background (Alu stacking)
    wps_background_regions: Option<PathBuf>, wps_background_output: Option<PathBuf>,
    wps_empty: bool,
    // FSD
    fsd_arms: Option<PathBuf>, fsd_output: Option<PathBuf>,
    // OCF
    ocf_regions: Option<PathBuf>, ocf_output: Option<PathBuf>,
    // Target regions for on/off-target split (panel data like MSK-ACCESS)
    target_regions_path: Option<PathBuf>,
    // Bait edge padding in bp (default 50, adaptive safety applies)
    bait_padding: u64,
    // Progress control
    silent: bool,
) -> PyResult<()> {
    use crate::gc_correction::{ReferenceData, ValidRegionFilter, GcObservationConsumer, compute_gcfix_factors, CorrectionFactors};
    use crate::engine;
    use std::sync::Arc;
    
    // ═══════════════════════════════════════════════════════════════════
    // Phase 1: GC Correction Factors (Load or Compute)
    // ═══════════════════════════════════════════════════════════════════
    let factors = if let Some(input_path) = correction_input_path {
        // Load pre-computed factors from CSV
        info!("Phase 1: Loading pre-computed GC correction factors from {:?}", input_path);
        match CorrectionFactors::load_csv(&input_path) {
            Ok(f) => {
                info!("Loaded {} correction factor entries", f.data.len());
                Some(f)
            },
            Err(e) => {
                warn!("Failed to load correction factors: {}. Proceeding without GC correction.", e);
                None
            }
        }
    } else if let (Some(gc_ref), Some(valid_reg), Some(out_factors)) = (gc_ref_path, valid_regions_path, correction_out_path) {
        // Compute factors from scratch
        info!("Phase 1: Computing GC correction factors...");
        
        // 1. Load Reference Data
        let ref_data = ReferenceData::load(&gc_ref)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to load GC reference: {}", e)))?;
            
        // 2. Load Valid Regions (needs mutable ChromosomeMap for loading)
        let mut valid_chrom_map = crate::bed::ChromosomeMap::new();
        let valid_regions = Arc::new(ValidRegionFilter::load(&valid_reg, &mut valid_chrom_map)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to load valid regions: {}", e)))?);
    
        // 3. Observed Counts
        let obs_consumer = GcObservationConsumer::new(valid_regions);
        
        // Run Engine
        let engine = engine::FragmentAnalyzer::new(obs_consumer, 500_000); // 500k batch
        let obs_result = engine.process_file(&bed_path, &mut crate::bed::ChromosomeMap::new(), silent)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Phase 1 failed: {}", e)))?;
            
        // 4. Compute Factors (use .observed field, not .counts)
        let computed = compute_gcfix_factors(&obs_result.observed, &ref_data, None)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Factor computation failed: {}", e)))?;
            
        // 5. Write Factors
        computed.write_csv(&out_factors)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to write factors: {}", e)))?;
            
        Some(computed)
    } else {
        info!("GC correction skipped (no factors provided and no reference files)");
        None
    };


    // ═══════════════════════════════════════════════════════════════════
    // Phase 2: Feature Extraction (with Correction)
    // ═══════════════════════════════════════════════════════════════════
    
    // 1. Initialize Chromosome Map (Shared)
    let mut chrom_map = ChromosomeMap::new();
    
    // 2. Initialize Consumers (with GC correction)
    let factors_arc = factors.map(Arc::new);
    
    let mut fsc_consumer = None;
    if let Some(bins) = fsc_bins {
        // Need to load regions. Reuse verify code?
        // fsc::parse_bin_file is public.
        let regions = crate::fsc::parse_bin_file(&bins).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        // Load target regions if provided (for on/off-target split)
        let target_tree = if let Some(ref target_path) = target_regions_path {
            let tree = crate::fsd::load_target_regions(target_path, &mut chrom_map)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            info!("FSC: Panel mode enabled, on/off-target split active");
            Some(Arc::new(tree))
        } else {
            None
        };
        
        fsc_consumer = Some(FscConsumer::new(&regions, &mut chrom_map, factors_arc.clone(), target_tree));
    }
    
    let mut wps_consumer = None;
    if let Some(regs) = wps_regions {
        // Parse regions using exposed helper
        let regions = crate::wps::parse_regions(&regs)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to parse WPS regions: {}", e)))?;
        
        // Create base consumer
        let mut consumer = WpsConsumer::new(regions, &mut chrom_map, factors_arc.clone());
        
        // Load BaitMask if target regions provided (for panel edge detection)
        if let Some(ref target_path) = target_regions_path {
            let bait_mask = crate::wps::BaitMask::from_bed(target_path, &mut chrom_map, bait_padding)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("Failed to load target regions for WPS: {}", e)))?;
            consumer = consumer.with_bait_mask(bait_mask);
            info!("WPS: BaitMask enabled, padding={}bp (adaptive safety applies)", bait_padding);
        }
        
        wps_consumer = Some(consumer);
    }
    
    let mut fsd_consumer = None;
    if let Some(arms) = fsd_arms {
        // Parse arms
        let regions = crate::fsd::parse_regions_file(&arms)
             .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        
        // Load target regions if provided (for on/off-target split)
        let target_tree = if let Some(ref target_path) = target_regions_path {
            let tree = crate::fsd::load_target_regions(target_path, &mut chrom_map)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            Some(Arc::new(tree))
        } else {
            None
        };
        
        fsd_consumer = Some(FsdConsumer::new(regions, &mut chrom_map, factors_arc.clone(), target_tree));
    }
    
    let mut ocf_consumer = None;
    if let Some(regs) = ocf_regions {
        // Load target regions if provided (for on/off-target split)
        let target_tree = if let Some(ref target_path) = target_regions_path {
            let tree = crate::fsd::load_target_regions(target_path, &mut chrom_map)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            info!("OCF: Panel mode enabled, on/off-target split active");
            Some(Arc::new(tree))
        } else {
            None
        };
        
        // OcfConsumer parses internally from path
        let consumer = OcfConsumer::new(&regs, &mut chrom_map, factors_arc.clone(), target_tree)
              .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        ocf_consumer = Some(consumer);
    }
    
    // WPS Background Consumer (Alu stacking with hierarchical groups)
    let mut wps_background_consumer = None;
    if let Some(bg_regions_path) = wps_background_regions {
        let regions = crate::wps::parse_alu_regions(&bg_regions_path)
             .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        let consumer = crate::wps::WpsBackgroundConsumer::new(regions, &mut chrom_map);
        wps_background_consumer = Some(consumer);
        info!("WPS Background: Enabled hierarchical Alu stacking");
    }
    
    // 3. Create MultiConsumer
    let consumer = MultiConsumer::new(fsc_consumer, wps_consumer, wps_background_consumer, fsd_consumer, ocf_consumer);
    
    // 4. Run Analysis
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map, silent)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    // 5. Write Outputs
    final_consumer.write_outputs(fsc_output, wps_output, wps_background_output, fsd_output, ocf_output, wps_empty)
         .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    
    Ok(())
}
