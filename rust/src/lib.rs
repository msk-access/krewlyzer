//! krewlyzer-core: High-performance Rust backend for cfDNA analysis
//! 
//! This crate provides fast implementations of fragment size analysis functions
//! for cell-free DNA analysis, exposed to Python via PyO3.

use pyo3::prelude::*;

// fsc, wps, fsd, ocf, mfsd, uxm maintained for individual tool access via legacy-compatible APIs (using new engines internally)
pub mod filters;
pub mod bed;
pub mod fsc;
pub mod fsd;
pub mod ocf;
pub mod wps;
pub mod mfsd;
pub mod uxm;
pub mod extract_motif;
pub mod engine;
pub mod pipeline;
pub mod gc_correction;

/// Read filtering configuration
#[pyclass]
#[derive(Clone, Debug)]
pub struct ReadFilters {
    #[pyo3(get, set)]
    pub mapq: u8,
    #[pyo3(get, set)]
    pub min_length: u32,
    #[pyo3(get, set)]
    pub max_length: u32,
    #[pyo3(get, set)]
    pub skip_duplicates: bool,
    #[pyo3(get, set)]
    pub require_proper_pair: bool,
}

#[pymethods]
impl ReadFilters {
    #[new]
    #[pyo3(signature = (mapq=20, min_length=65, max_length=400, skip_duplicates=true, require_proper_pair=true))]
    fn new(
        mapq: u8,
        min_length: u32,
        max_length: u32,
        skip_duplicates: bool,
        require_proper_pair: bool,
    ) -> Self {
        Self {
            mapq,
            min_length,
            max_length,
            skip_duplicates,
            require_proper_pair,
        }
    }
}

#[pymodule]
#[pyo3(name = "_core")]
fn krewlyzer_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize logging (bridges Rust log -> Python logging)
    let _ = pyo3_log::init();

    m.add_class::<ReadFilters>()?;
    
    // Thread configuration
    m.add_function(wrap_pyfunction!(configure_threads, m)?)?;
    
    // FSC functions
    m.add_function(wrap_pyfunction!(fsc::count_fragments_by_bins, m)?)?;
    m.add_function(wrap_pyfunction!(fsc::count_fragments_gc_corrected, m)?)?;

    // FSD submodule
    let fsd_mod = PyModule::new(m.py(), "fsd")?;
    fsd_mod.add_function(wrap_pyfunction!(fsd::calculate_fsd, &fsd_mod)?)?;
    m.add_submodule(&fsd_mod)?;

    // WPS submodule (also exposed as function above? Cleaned up duplication)
    let wps_mod = PyModule::new(m.py(), "wps")?;
    wps_mod.add_function(wrap_pyfunction!(wps::calculate_wps, &wps_mod)?)?;
    m.add_submodule(&wps_mod)?;

    // OCF submodule
    let ocf_mod = PyModule::new(m.py(), "ocf")?;
    ocf_mod.add_function(wrap_pyfunction!(ocf::calculate_ocf, &ocf_mod)?)?;
    m.add_submodule(&ocf_mod)?;
    
    // mFSD submodule
    let mfsd_mod = PyModule::new(m.py(), "mfsd")?;
    mfsd_mod.add_function(wrap_pyfunction!(mfsd::calculate_mfsd, &mfsd_mod)?)?;
    m.add_submodule(&mfsd_mod)?;
    
    // UXM submodule
    let uxm_mod = PyModule::new(m.py(), "uxm")?;
    uxm_mod.add_function(wrap_pyfunction!(uxm::calculate_uxm, &uxm_mod)?)?;
    m.add_submodule(&uxm_mod)?;
    
    // Unified Engine (Extract+Motif)
    let extract_motif_mod = PyModule::new(m.py(), "extract_motif")?;
    extract_motif_mod.add_function(wrap_pyfunction!(extract_motif::process_bam_parallel, &extract_motif_mod)?)?;
    m.add_submodule(&extract_motif_mod)?;
    
    // Unified Pipeline (FSC/FSD/WPS/OCF)
    m.add_function(wrap_pyfunction!(pipeline::run_unified_pipeline, m)?)?;
    

    


    
    // Version
    #[pyfn(m)]
    fn version() -> &'static str {
        env!("CARGO_PKG_VERSION")
    }
    
    Ok(())
}

/// Configure the number of threads for Rayon parallel processing.
/// Must be called before any parallel operations.
/// 
/// # Arguments
/// * `num_threads` - Number of threads to use (0 = use all available cores)
#[pyfunction]
#[pyo3(signature = (num_threads=0))]
fn configure_threads(num_threads: usize) -> PyResult<()> {
    let threads = if num_threads == 0 {
        num_cpus::get()
    } else {
        num_threads
    };
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
            format!("Failed to configure thread pool: {}. Note: can only be called once.", e)
        ))?;
    
    Ok(())
}
