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
    
    FSC 5 non-overlapping size channels:
    - Ultra-short: 65-100bp
    - Core-short: 101-149bp
    - Mono-nucleosomal: 150-220bp
    - Di-nucleosomal: 221-260bp
    - Long: 261-400bp
    """
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # 80bp - counts as ultra_short
        f.write("chr1\t100\t180\t0.5\n")
        # 170bp - counts as mono_nucl
        f.write("chr1\t200\t370\t0.5\n")
        # 300bp - counts as long
        f.write("chr1\t500\t800\t0.5\n")
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
    
    ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Verify counts
    assert ultra_shorts[0] == 1  # 80bp (65-100)
    assert core_shorts[0] == 0   # None (101-149)
    assert mono_nucls[0] == 1    # 170bp (150-220)
    assert di_nucls[0] == 0      # None (221-260)
    assert longs[0] == 1         # 300bp (261-400)
    assert totals[0] == 3        # Total


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_multiple_bins(tmp_path):
    """Test FSC counting across multiple bins with 5 channels."""
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # Bin 1 (0-1000): 2 fragments
        f.write("chr1\t100\t220\t0.5\n")  # 120bp - core_short (101-149)
        f.write("chr1\t300\t500\t0.5\n")  # 200bp - mono_nucl (150-220)
        # Bin 2 (1000-2000): 1 fragment
        f.write("chr1\t1100\t1400\t0.5\n")  # 300bp - long (261-400)
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
        f.write("chr1\t1000\t2000\n")
    
    ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Bin 1
    assert totals[0] == 2
    assert core_shorts[0] == 1   # 120bp (101-149)
    assert mono_nucls[0] == 1    # 200bp (150-220)
    
    # Bin 2
    assert totals[1] == 1
    assert longs[1] == 1         # 300bp (261-400)


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_empty_bin(tmp_path):
    """Test FSC handles empty bins correctly."""
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t100\t200\t0.5\n")  # 100bp - ultra_short
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t50\n")     # Empty bin
        f.write("chr1\t50\t500\n")   # Bin with fragment
    
    ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    assert totals[0] == 0  # Empty bin
    assert totals[1] == 1  # Bin with fragment
