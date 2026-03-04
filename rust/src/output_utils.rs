//! Shared output utilities for Krewlyzer Rust modules.
//!
//! Provides two boilerplate-free write functions used by all output modules:
//! - [`write_parquet_batch`]: Write Arrow arrays to Parquet (Snappy compressed).
//! - [`write_tsv_batch`]: Write header + rows to TSV, with optional gzip.
//! - [`tsv_path`]: Resolve `.tsv` vs `.tsv.gz` output path.
//!
//! ## Design
//! Every Rust output module builds its schema + typed arrays and calls these
//! utilities — no repeated `ArrowWriter` or `BufWriter`/`GzEncoder` boilerplate.
//!
//! ## Usage
//! ```rust
//! use crate::output_utils::{write_parquet_batch, write_tsv_batch, tsv_path};
//!
//! // Parquet write
//! write_parquet_batch(&output_path, Arc::new(schema), arrays)?;
//!
//! // TSV write (with optional gzip)
//! write_tsv_batch(&tsv_path(&base_path, compress), &["col1", "col2"], rows, compress)?;
//! ```

use arrow::array::ArrayRef;
use arrow::datatypes::SchemaRef;
use arrow::record_batch::RecordBatch;
use flate2::{write::GzEncoder, Compression as GzCompression};
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};
use log::{debug, info};

/// Write typed Arrow arrays to a Parquet file (Snappy compressed).
///
/// Constructs a single `RecordBatch` from the provided schema and arrays,
/// then writes it using `ArrowWriter` with Snappy compression.
///
/// # Arguments
/// * `output_path` - Destination `.parquet` file path (created or overwritten).
/// * `schema` - Arrow schema matching the number and types of `arrays`.
/// * `arrays` - Column data; must align 1-to-1 with `schema` fields.
///
/// # Errors
/// Returns an error if the file cannot be created, the batch is invalid, or writing fails.
///
/// # Example
/// ```rust
/// use arrow::array::{ArrayRef, Float64Array};
/// use arrow::datatypes::{DataType, Field, Schema};
/// use std::sync::Arc;
/// use crate::output_utils::write_parquet_batch;
///
/// let schema = Schema::new(vec![Field::new("score", DataType::Float64, false)]);
/// let arrays: Vec<ArrayRef> = vec![Arc::new(Float64Array::from(vec![1.0, 2.0]))];
/// write_parquet_batch(Path::new("out.parquet"), Arc::new(schema), arrays).unwrap();
/// ```
pub fn write_parquet_batch(output_path: &Path, schema: SchemaRef, arrays: Vec<ArrayRef>) -> Result<()> {
    debug!("output_utils: writing Parquet → {:?}", output_path);

    let batch = RecordBatch::try_new(schema.clone(), arrays)
        .context("Building RecordBatch for Parquet output")?;

    let file = std::fs::File::create(output_path)
        .with_context(|| format!("Creating Parquet file: {:?}", output_path))?;

    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
        .context("Creating ArrowWriter")?;

    writer.write(&batch).context("Writing RecordBatch to Parquet")?;
    writer.close().context("Finalizing Parquet file")?;

    info!("output_utils: wrote {} rows → {:?}", batch.num_rows(), output_path);
    Ok(())
}

/// Write a header row and data rows to a TSV file, with optional gzip compression.
///
/// Uses `BufWriter` for efficiency. When `compress` is true, wraps in a `GzEncoder`;
/// the caller is responsible for passing a path with `.gz` extension (use [`tsv_path`]).
///
/// # Arguments
/// * `output_path` - Destination file path. Pass result of [`tsv_path`] for correct extension.
/// * `headers` - Column names written as the first line (tab-separated).
/// * `rows` - Data rows, each a `Vec<String>` aligned with `headers`.
/// * `compress` - If `true`, wraps output in gzip; path should end in `.gz`.
///
/// # Errors
/// Returns an error if the file cannot be created or if any write operation fails.
pub fn write_tsv_batch(
    output_path: &Path,
    headers: &[&str],
    rows: Vec<Vec<String>>,
    compress: bool,
) -> Result<()> {
    debug!(
        "output_utils: writing TSV ({} rows, compress={}) → {:?}",
        rows.len(), compress, output_path
    );

    let file = std::fs::File::create(output_path)
        .with_context(|| format!("Creating TSV file: {:?}", output_path))?;

    // Dynamic dispatch: gzip or plain buffered writer
    let mut writer: Box<dyn Write> = if compress {
        Box::new(GzEncoder::new(BufWriter::new(file), GzCompression::default()))
    } else {
        Box::new(BufWriter::new(file))
    };

    writeln!(writer, "{}", headers.join("\t"))
        .context("Writing TSV header")?;

    for row in &rows {
        writeln!(writer, "{}", row.join("\t"))
            .context("Writing TSV row")?;
    }

    // Explicitly flush so GzEncoder finalizes the gzip stream
    writer.flush().context("Flushing TSV writer")?;

    info!(
        "output_utils: wrote {} rows → {:?}",
        rows.len(), output_path
    );
    Ok(())
}

/// Resolve the output path for a TSV file, respecting the `--compress` flag.
///
/// Appends `.gz` to the path when `compress` is true. The base path should
/// already have the `.tsv` extension; this function appends `.gz` on top.
///
/// # Example
/// ```rust
/// use crate::output_utils::tsv_path;
/// use std::path::Path;
///
/// let p = tsv_path(Path::new("sample.FSD.tsv"), false); // → "sample.FSD.tsv"
/// let p = tsv_path(Path::new("sample.FSD.tsv"), true);  // → "sample.FSD.tsv.gz"
/// ```
pub fn tsv_path(base: &Path, compress: bool) -> PathBuf {
    if compress {
        // Append .gz to the full path (e.g. sample.FSD.tsv → sample.FSD.tsv.gz)
        PathBuf::from(format!("{}.gz", base.display()))
    } else {
        base.to_path_buf()
    }
}

/// Determine output_format from a string slice, with a safe default of "tsv".
///
/// Centralises validation so all modules handle unknown values consistently.
/// Unknown values fall back to "tsv" with a warning logged.
///
/// # Returns
/// One of `"tsv"`, `"parquet"`, or `"both"`.
pub fn validated_output_format(output_format: &str) -> &str {
    match output_format {
        "tsv" | "parquet" | "both" => output_format,
        other => {
            log::warn!(
                "output_utils: unknown output_format '{}'; defaulting to 'tsv'",
                other
            );
            "tsv"
        }
    }
}

/// Returns true if TSV output should be written for the given format.
#[inline]
pub fn should_write_tsv(output_format: &str) -> bool {
    output_format == "tsv" || output_format == "both"
}

/// Returns true if Parquet output should be written for the given format.
#[inline]
pub fn should_write_parquet(output_format: &str) -> bool {
    output_format == "parquet" || output_format == "both"
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::Float64Array;
    use arrow::datatypes::{DataType, Field, Schema};
    use tempfile::tempdir;

    #[test]
    fn test_tsv_path_no_compress() {
        let base = Path::new("sample.FSD.tsv");
        assert_eq!(tsv_path(base, false), PathBuf::from("sample.FSD.tsv"));
    }

    #[test]
    fn test_tsv_path_compress() {
        let base = Path::new("sample.FSD.tsv");
        assert_eq!(tsv_path(base, true), PathBuf::from("sample.FSD.tsv.gz"));
    }

    #[test]
    fn test_write_tsv_batch_plain() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.tsv");
        let headers = &["a", "b"];
        let rows = vec![vec!["1".into(), "2.5".into()], vec!["3".into(), "4.0".into()]];
        write_tsv_batch(&path, headers, rows, false).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with("a\tb\n"));
        assert!(content.contains("1\t2.5\n"));
    }

    #[test]
    fn test_write_tsv_batch_gzip() {
        use flate2::read::GzDecoder;
        use std::io::Read;
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.tsv.gz");
        let headers = &["x"];
        let rows = vec![vec!["42".into()]];
        write_tsv_batch(&path, headers, rows, true).unwrap();
        let f = std::fs::File::open(&path).unwrap();
        let mut gz = GzDecoder::new(f);
        let mut content = String::new();
        gz.read_to_string(&mut content).unwrap();
        assert!(content.contains("x\n"));
        assert!(content.contains("42\n"));
    }

    #[test]
    fn test_write_parquet_batch() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.parquet");
        let schema = Schema::new(vec![Field::new("score", DataType::Float64, false)]);
        let arrays: Vec<ArrayRef> = vec![Arc::new(Float64Array::from(vec![1.0, 2.5, 3.0]))];
        write_parquet_batch(&path, Arc::new(schema), arrays).unwrap();
        assert!(path.exists());
    }

    #[test]
    fn test_validated_output_format() {
        assert_eq!(validated_output_format("tsv"), "tsv");
        assert_eq!(validated_output_format("parquet"), "parquet");
        assert_eq!(validated_output_format("both"), "both");
        assert_eq!(validated_output_format("unknown"), "tsv");
    }

    #[test]
    fn test_should_write_flags() {
        assert!(should_write_tsv("tsv"));
        assert!(should_write_tsv("both"));
        assert!(!should_write_tsv("parquet"));
        assert!(should_write_parquet("parquet"));
        assert!(should_write_parquet("both"));
        assert!(!should_write_parquet("tsv"));
    }
}
