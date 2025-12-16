//! Unified Single-Pass Engine for krewlyzer
//! 
//! This module coordinates the parallel processing of fragments through multiple "Consumers" (FSC, FSR, WPS, etc).
//! It uses Rayon for Map-Reduce parallelism.

use rayon::prelude::*;
use anyhow::{Result, Context};
use std::path::Path;
use std::fs::File;
use std::io::BufRead;
use noodles::bgzf;

use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;

use crate::bed::{Fragment, ChromosomeMap};

/// Trait for any consumer that wants to process fragments
/// 
/// Consumers must be thread-safe (Send + Sync + Clone).
/// Each thread gets a clone of the consumer to maintain thread-local state.
/// At the end, thread-local consumers are merged.
pub trait FragmentConsumer: Send + Sync + Clone {
    /// Name of the consumer (for logging)
    fn name(&self) -> &str;

    /// Process a single fragment.
    /// This is called in a tight loop, so it must be fast.
    /// The consumer should update its internal state (counters, histograms, etc).
    fn consume(&mut self, fragment: &Fragment);

    /// Merge results from another consumer (thread-local reduction).
    /// `other` is the partial result from another thread.
    fn merge(&mut self, other: Self);
}

/// The Main Engine
pub struct FragmentAnalyzer<C: FragmentConsumer> {
    consumer_template: C,
    chunk_size: usize,
}

impl<C: FragmentConsumer> FragmentAnalyzer<C> {
    pub fn new(consumer: C, chunk_size: usize) -> Self {
        Self {
            consumer_template: consumer,
            chunk_size,
        }
    }

    /// Process the BED file in parallel
    pub fn process_file(&self, bed_path: &Path, chrom_map: &mut ChromosomeMap) -> Result<C> {
        let file = File::open(bed_path)
            .with_context(|| format!("Failed to open BED file: {:?}", bed_path))?;
        let mut reader = bgzf::io::Reader::new(file);
        
        let mut final_consumer = self.consumer_template.clone();
        
        // Progress Bar (Spinner since total count unknown without full scan)
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"));
        pb.set_message("Processing fragments...");
        pb.enable_steady_tick(Duration::from_millis(100));

        // Buffer for current chunk
        let mut fragment_buffer: Vec<Fragment> = Vec::with_capacity(self.chunk_size);
        let mut line_buf = String::new();
        let mut total_processed: u64 = 0;
        
        loop {
            fragment_buffer.clear();
            
            // Fill chunk
            for _ in 0..self.chunk_size {
                line_buf.clear();
                match reader.read_line(&mut line_buf) {
                    Ok(0) => break, // EOF
                    Ok(_) => {
                        let line = line_buf.trim_end();
                        if line.is_empty() { continue; }
                        
                        // Parse fragment
                        if let Some(frag) = parse_line(line, chrom_map) {
                            fragment_buffer.push(frag);
                        }
                    },
                    Err(_) => break,
                }
            }
            
            if fragment_buffer.is_empty() {
                break;
            }
            
            let chunk_len = fragment_buffer.len() as u64;

            // Process chunk in parallel using Rayon's `fold` + `reduce` pattern
            let chunk_consumer = fragment_buffer.par_iter()
                .fold(
                    || self.consumer_template.clone(), // Init thread-local consumer from template
                    |mut c, frag| {
                        c.consume(frag);
                        c
                    }
                )
                .reduce(
                    || self.consumer_template.clone(), // Identity?
                    // Actually reduce expects (A, A) -> A.
                    // If we use reduce, we need a way to combine.
                    |mut a, b| {
                        a.merge(b);
                        a
                    }
                );
            
            // Merge chunk result into final result
            final_consumer.merge(chunk_consumer);
            
            // Update Progress
            total_processed += chunk_len;
            pb.set_message(format!("Processed {} fragments...", total_processed));
        }
        
        pb.finish_with_message(format!("Done! Processed {} fragments.", total_processed));
        
        Ok(final_consumer)
    }
}

fn parse_line(line: &str, chrom_map: &mut ChromosomeMap) -> Option<Fragment> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 4 { return None; }
    
    let chrom_str = fields[0];
    let start: u64 = fields[1].parse().ok()?;
    let end: u64 = fields[2].parse().ok()?;
    let gc: f64 = fields[3].parse().unwrap_or(0.0);
    
    let length = end.saturating_sub(start);
    
    // Normalize chromosome (remove "chr" if needed? Or let Consumer handle logic?)
    // Best to normalize globally here.
    let chrom_norm = chrom_str.trim_start_matches("chr");
    
    let chrom_id = chrom_map.get_id(chrom_norm);
    
    Some(Fragment {
        chrom_id,
        start,
        end,
        length,
        gc,
    })
}
