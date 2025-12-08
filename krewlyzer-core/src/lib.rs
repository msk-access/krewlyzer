//! krewlyzer-core: High-performance Rust backend for cfDNA analysis
//! 
//! This crate provides fast implementations of fragment size analysis functions
//! for cell-free DNA analysis, exposed to Python via PyO3.

use pyo3::prelude::*;

pub mod filters;
pub mod bed;
pub mod fsc;

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

/// Python module initialization
#[pymodule]
fn krewlyzer_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReadFilters>()?;
    
    // FSC functions
    m.add_function(wrap_pyfunction!(fsc::count_fragments_by_bins, m)?)?;
    
    // Version
    #[pyfn(m)]
    fn version() -> &'static str {
        env!("CARGO_PKG_VERSION")
    }
    
    Ok(())
}
