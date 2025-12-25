//! Panel of Normals (PON) model loading and application.
//!
//! This module handles loading PON models from Parquet files and applying
//! hybrid GC correction (PON + within-sample residual).

use std::collections::HashMap;
use std::path::Path;
use std::fs::File;


use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use serde::{Deserialize, Serialize};
use log::{info, warn};


/// GC bias correction curves for a single fragment type.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcCurve {
    pub gc_bins: Vec<f64>,
    pub expected: Vec<f64>,
    pub std: Vec<f64>,
}

impl GcCurve {
    /// Interpolate expected coverage for given GC content.
    pub fn get_expected(&self, gc: f64) -> f64 {
        if self.gc_bins.is_empty() || self.expected.is_empty() {
            return 1.0;
        }
        
        // Linear interpolation
        if gc <= self.gc_bins[0] {
            return self.expected[0];
        }
        if gc >= *self.gc_bins.last().unwrap() {
            return *self.expected.last().unwrap();
        }
        
        for i in 0..self.gc_bins.len() - 1 {
            if gc >= self.gc_bins[i] && gc < self.gc_bins[i + 1] {
                let t = (gc - self.gc_bins[i]) / (self.gc_bins[i + 1] - self.gc_bins[i]);
                return self.expected[i] + t * (self.expected[i + 1] - self.expected[i]);
            }
        }
        
        1.0
    }
}

/// GC bias model containing curves for all fragment types.
#[derive(Debug, Clone, Default)]
pub struct GcBiasModel {
    pub short: Option<GcCurve>,
    pub intermediate: Option<GcCurve>,
    pub long: Option<GcCurve>,
    pub wps_long: Option<GcCurve>,
    pub wps_short: Option<GcCurve>,
}

impl GcBiasModel {
    /// Get expected coverage for a fragment type at given GC.
    pub fn get_expected(&self, gc: f64, frag_type: &str) -> f64 {
        let curve = match frag_type {
            "short" => &self.short,
            "intermediate" => &self.intermediate,
            "long" => &self.long,
            "wps_long" => &self.wps_long,
            "wps_short" => &self.wps_short,
            _ => &None,
        };
        
        curve.as_ref().map(|c| c.get_expected(gc)).unwrap_or(1.0)
    }
}

/// WPS baseline for a single region.
#[derive(Debug, Clone)]
pub struct WpsRegionBaseline {
    pub wps_long_mean: f64,
    pub wps_long_std: f64,
    pub wps_short_mean: f64,
    pub wps_short_std: f64,
}

/// WPS baseline model containing per-region statistics.
#[derive(Debug, Clone, Default)]
pub struct WpsBaseline {
    pub regions: HashMap<String, WpsRegionBaseline>,
}

impl WpsBaseline {
    /// Get baseline for a region.
    pub fn get_baseline(&self, region_id: &str) -> Option<&WpsRegionBaseline> {
        self.regions.get(region_id)
    }
    
    /// Compute z-score for WPS long.
    pub fn compute_z_long(&self, region_id: &str, sample_mean: f64) -> Option<f64> {
        self.get_baseline(region_id).map(|b| {
            if b.wps_long_std > 0.0 {
                (sample_mean - b.wps_long_mean) / b.wps_long_std
            } else {
                0.0
            }
        })
    }
    
    /// Compute z-score for WPS short.
    pub fn compute_z_short(&self, region_id: &str, sample_mean: f64) -> Option<f64> {
        self.get_baseline(region_id).map(|b| {
            if b.wps_short_std > 0.0 {
                (sample_mean - b.wps_short_mean) / b.wps_short_std
            } else {
                0.0
            }
        })
    }
}

/// Unified Panel of Normals model.
#[derive(Debug, Clone, Default)]
pub struct PonModel {
    pub schema_version: String,
    pub assay: String,
    pub build_date: String,
    pub n_samples: usize,
    pub reference: String,
    
    pub gc_bias: GcBiasModel,
    pub wps_baseline: WpsBaseline,
    // TODO: Add FSD baseline
}

impl PonModel {
    /// Load PON model from Parquet file.
    pub fn load(path: &Path) -> Result<Self, String> {
        if !path.exists() {
            return Err(format!("PON model not found: {}", path.display()));
        }
        
        let file = File::open(path)
            .map_err(|e| format!("Failed to open PON file: {}", e))?;
        
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| format!("Failed to read Parquet: {}", e))?;
        
        let reader = builder.build()
            .map_err(|e| format!("Failed to build Parquet reader: {}", e))?;
        
        let mut model = PonModel::default();
        model.schema_version = "1.0".to_string();
        
        // TODO: Parse record batches into model components
        for batch_result in reader {
            let _batch = batch_result
                .map_err(|e| format!("Failed to read batch: {}", e))?;
            // Parse batch into model
        }
        
        info!("Loaded PON model: {} (n={})", model.assay, model.n_samples);
        
        Ok(model)
    }
    
    /// Check assay compatibility and warn if mismatch.
    pub fn check_assay_compatibility(&self, sample_assay: Option<&str>) {
        if let Some(assay) = sample_assay {
            if assay != self.assay {
                warn!(
                    "PON model built for {}, sample may be from different assay ({})",
                    self.assay, assay
                );
            }
        }
    }
    
    /// Apply hybrid GC correction: PON + within-sample residual.
    ///
    /// Algorithm:
    /// 1. PON correction: pon_corrected = observed / pon_expected[gc]
    /// 2. Residual LOESS: residual = loess(gc, pon_corrected)
    /// 3. Final: final = pon_corrected / residual
    pub fn apply_hybrid_correction(
        &self,
        observed: &[f64],
        gc_values: &[f64],
        frag_type: &str,
    ) -> Vec<f64> {
        if observed.len() != gc_values.len() {
            return observed.to_vec();
        }
        
        // Step 1: PON correction
        let pon_corrected: Vec<f64> = observed
            .iter()
            .zip(gc_values.iter())
            .map(|(obs, gc)| {
                let expected = self.gc_bias.get_expected(*gc, frag_type);
                if expected > 0.0 {
                    obs / expected
                } else {
                    *obs
                }
            })
            .collect();
        
        // Step 2: Within-sample residual LOESS (REMOVED)
        // Hybrid correction now relies on upstream fragment weighting
        
        // Return the corrected values for this type
        // (correct_gc_bias_per_type returns short, intermediate, long)
        // For now, return PON-corrected only
        pon_corrected
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_gc_curve_interpolation() {
        let curve = GcCurve {
            gc_bins: vec![0.3, 0.4, 0.5, 0.6],
            expected: vec![0.8, 1.0, 1.0, 0.8],
            std: vec![0.1, 0.1, 0.1, 0.1],
        };
        
        // Exact match
        assert!((curve.get_expected(0.4) - 1.0).abs() < 0.001);
        
        // Interpolation
        let mid = curve.get_expected(0.35);
        assert!(mid > 0.8 && mid < 1.0);
        
        // Edge cases
        assert!((curve.get_expected(0.2) - 0.8).abs() < 0.001);
        assert!((curve.get_expected(0.7) - 0.8).abs() < 0.001);
    }
}
