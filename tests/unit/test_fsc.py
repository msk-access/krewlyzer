"""
Unit tests for FSC (Fragment Size Coverage) calculations.

Tests the Rust-backed FSC counting functions directly.
"""
import pytest
import pysam
from pathlib import Path

from krewlyzer import _core


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_count_fragments_by_bins(tmp_path):
    """Test FSC fragment counting via Rust backend.
    
    FSC size categories (cumulative):
    - Ultra-short: ≤100bp
    - Short: ≤149bp (includes ultra-short)
    - Intermediate: 151-259bp
    - Long: ≥260bp
    """
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # 80bp - counts as ultra_short AND short
        f.write("chr1\t100\t180\t0.5\n")
        # 200bp - counts as intermediate
        f.write("chr1\t200\t400\t0.5\n")
        # 300bp - counts as long
        f.write("chr1\t500\t800\t0.5\n")
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
    
    ultra_shorts, shorts, ints, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Verify counts
    assert ultra_shorts[0] == 1  # 80bp
    assert shorts[0] == 1        # 80bp (cumulative but only 1 qualifying)
    assert ints[0] == 1          # 200bp
    assert longs[0] == 1         # 300bp
    assert totals[0] == 3        # Total


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_multiple_bins(tmp_path):
    """Test FSC counting across multiple bins."""
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # Bin 1 (0-1000): 2 fragments
        f.write("chr1\t100\t220\t0.5\n")  # 120bp - short
        f.write("chr1\t300\t500\t0.5\n")  # 200bp - intermediate
        # Bin 2 (1000-2000): 1 fragment
        f.write("chr1\t1100\t1400\t0.5\n")  # 300bp - long
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
        f.write("chr1\t1000\t2000\n")
    
    ultra_shorts, shorts, ints, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Bin 1
    assert totals[0] == 2
    assert shorts[0] == 1   # 120bp
    assert ints[0] == 1     # 200bp
    
    # Bin 2
    assert totals[1] == 1
    assert longs[1] == 1    # 300bp


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_empty_bin(tmp_path):
    """Test FSC handles empty bins correctly."""
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t100\t200\t0.5\n")
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t50\n")     # Empty bin
        f.write("chr1\t50\t500\n")   # Bin with fragment
    
    ultra_shorts, shorts, ints, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    assert totals[0] == 0  # Empty bin
    assert totals[1] == 1  # Bin with fragment
