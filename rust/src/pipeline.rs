use pyo3::prelude::*;
use std::path::PathBuf;
use anyhow::{Result, Context};

use crate::bed::{Fragment, ChromosomeMap};
use crate::engine::{FragmentConsumer, FragmentAnalyzer};
use crate::fsc::FscConsumer;
use crate::wps::WpsConsumer;
use crate::fsd::FsdConsumer;
use crate::ocf::OcfConsumer;

/// A consumer that delegates to multiple other consumers
#[derive(Clone)]
pub struct MultiConsumer {
    fsc: Option<FscConsumer>,
    wps: Option<WpsConsumer>,
    fsd: Option<FsdConsumer>,
    ocf: Option<OcfConsumer>,
}

impl MultiConsumer {
    pub fn new(
        fsc: Option<FscConsumer>,
        wps: Option<WpsConsumer>,
        fsd: Option<FsdConsumer>,
        ocf: Option<OcfConsumer>,
    ) -> Self {
        Self { fsc, wps, fsd, ocf }
    }
    
    pub fn write_outputs(
        &self,
        fsc_path: Option<PathBuf>,
        wps_path: Option<PathBuf>, // Directory or file pattern? WPS takes bed and tsv_file and writes to out_dir relative... no wait. 
        // calculate_wps takes `output_file`, which is a pattern "%s.tsv.gz".
        // wps consumer `write_output` takes `output_pattern`.
        fsd_path: Option<PathBuf>,
        ocf_dir: Option<PathBuf>,
        empty: bool, // for WPS
    ) -> Result<()> {
        if let Some(c) = &self.fsc {
            if let Some(p) = fsc_path {
                c.write_output(&p).context("Writing FSC output")?;
            }
        }
        
        if let Some(c) = &self.wps {
            if let Some(p) = wps_path {
                // p is likely file pattern if logic matches calculate_wps?
                // calculate_wps takes "output_dir" and "file_stem".
                // And constructs filename.
                // But wrapper might pass full path?
                // `run_unified_pipeline` takes `wps_output`.
                // If it is full path, we use it.
                // WpsConsumer::write_output just writes to path.
                c.write_output(&p, None, empty, false, false).context("Writing WPS output")?;
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
        if let Some(c) = &mut self.fsd { c.consume(fragment); }
        if let Some(c) = &mut self.ocf { c.consume(fragment); }
    }

    fn merge(&mut self, other: Self) {
        if let (Some(a), Some(b)) = (&mut self.fsc, other.fsc) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.wps, other.wps) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.fsd, other.fsd) { a.merge(b); }
        if let (Some(a), Some(b)) = (&mut self.ocf, other.ocf) { a.merge(b); }
    }
}

#[pyfunction]
#[pyo3(signature = (
    bed_path,
    // FSC Args
    fsc_bins=None, fsc_output=None,
    // WPS Args
    wps_regions=None, wps_output=None, wps_empty=false,
    // FSD Args
    fsd_arms=None, fsd_output=None,
    // OCF Args
    ocf_regions=None, ocf_output=None
))]
pub fn run_unified_pipeline(
    _py: Python,
    bed_path: PathBuf,
    // FSC
    fsc_bins: Option<PathBuf>, fsc_output: Option<PathBuf>,
    // WPS
    wps_regions: Option<PathBuf>, wps_output: Option<PathBuf>, wps_empty: bool,
    // FSD
    fsd_arms: Option<PathBuf>, fsd_output: Option<PathBuf>,
    // OCF
    ocf_regions: Option<PathBuf>, ocf_output: Option<PathBuf>,
) -> PyResult<()> {
    
    // 1. Initialize Chromosome Map (Shared)
    let mut chrom_map = ChromosomeMap::new();
    
    // 2. Initialize Consumers
    let mut fsc_consumer = None;
    if let Some(bins) = fsc_bins {
        // Need to load regions. Reuse verify code?
        // fsc::parse_bin_file is public.
        let regions = crate::fsc::parse_bin_file(&bins).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        fsc_consumer = Some(FscConsumer::new(&regions, &mut chrom_map));
    }
    
    let mut wps_consumer = None;
    if let Some(regs) = wps_regions {
        // Parse regions using exposed helper
        let regions = crate::wps::parse_regions(&regs)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to parse WPS regions: {}", e)))?;
        wps_consumer = Some(WpsConsumer::new(regions, &mut chrom_map));
    }
    
    let mut fsd_consumer = None;
    if let Some(arms) = fsd_arms {
        // Parse arms
        let regions = crate::fsd::parse_regions_file(&arms)
             .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        fsd_consumer = Some(FsdConsumer::new(regions, &mut chrom_map));
    }
    
    let mut ocf_consumer = None;
    if let Some(regs) = ocf_regions {
        // OcfConsumer parses internally from path
        let consumer = OcfConsumer::new(&regs, &mut chrom_map)
              .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        ocf_consumer = Some(consumer);
    }
    
    // 3. Create MultiConsumer
    let consumer = MultiConsumer::new(fsc_consumer, wps_consumer, fsd_consumer, ocf_consumer);
    
    // 4. Run Analysis
    let analyzer = FragmentAnalyzer::new(consumer, 100_000);
    let final_consumer = analyzer.process_file(&bed_path, &mut chrom_map)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        
    // 5. Write Outputs
    final_consumer.write_outputs(fsc_output, wps_output, fsd_output, ocf_output, wps_empty)
         .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    
    Ok(())
}
