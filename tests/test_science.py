"""
Science/calculation tests for krewlyzer.
These tests validate the core scientific calculations.
"""
import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import pysam
import os

from krewlyzer.helpers import reverse_complement, get_End_motif

# Rust backend is now required
import krewlyzer_core


def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("A") == "T"
    assert reverse_complement("") == ""
    assert reverse_complement("N") == "N"


def test_get_End_motif():
    # Initialize dict with keys because get_End_motif only updates existing keys
    d = {"AAA": 0, "CCC": 0, "GGG": 0}
    
    # Correct order: Emotif, end5, end3
    get_End_motif(d, "AAA", "CCC")
    assert d["AAA"] == 1
    assert d["CCC"] == 1
    
    get_End_motif(d, "AAA", "GGG")
    assert d["AAA"] == 2
    assert d["GGG"] == 1


def test_fsc_calculation(tmp_path):
    """Test FSC calculation using Rust backend."""
    # Create a dummy .bed.gz file with 2 fragments in different bins
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t100\t200\t0.5\n")  # Window 1: 1 short (len 100)
        f.write("chr1\t1100\t1300\t0.5\n")  # Window 2: 1 intermediate (len 200)
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    # Create bins file
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
        f.write("chr1\t1000\t2000\n")
    
    # Call Rust backend directly
    ultra_shorts, shorts, ints, longs, totals, gcs = krewlyzer_core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Verify counts
    # Window 1: 1 short fragment (len 100)
    # Window 2: 1 intermediate fragment (len 200)
    assert shorts[0] == 1  # Window 1 has 1 short
    assert shorts[1] == 0  # Window 2 has 0 shorts
    assert ints[0] == 0    # Window 1 has 0 intermediate
    assert ints[1] == 1    # Window 2 has 1 intermediate
    assert totals[0] == 1  # Window 1: 1 total
    assert totals[1] == 1  # Window 2: 1 total


def test_wps_calculation(tmp_path):
    """Test WPS calculation using Rust backend."""
    # Create BED file with reads
    bed_file = tmp_path / "wps.bed"
    with open(bed_file, "w") as f:
        # Long fragment (150bp) - should appear in wps_long
        f.write("chr1\t50\t200\t0.5\n")
        # Short fragment (50bp) - should appear in wps_short
        f.write("chr1\t100\t150\t0.5\n")
        
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    # Create TSV regions file
    tsv_file = tmp_path / "regions.tsv"
    with open(tsv_file, "w") as f:
        # id chrom start end strand
        f.write("region1\tchr1\t100\t200\t+\n")
    
    # Call Rust WPS calculation (new unified API)
    count = krewlyzer_core.wps.calculate_wps(
        bedgz,
        str(tsv_file),
        str(tmp_path),
        "wps",  # file_stem
        False,  # empty
        None    # total_fragments
    )
    
    assert count > 0  # At least one region processed
    
    # Check output file exists
    out_file = tmp_path / "wps.WPS.tsv.gz"
    assert out_file.exists()
    
    # Read and verify unified output structure
    df = pd.read_csv(out_file, sep="\t")
    assert len(df) > 0
    assert 'gene_id' in df.columns
    assert 'wps_long' in df.columns
    assert 'wps_short' in df.columns
    assert 'wps_ratio' in df.columns
    assert df['gene_id'].iloc[0] == 'region1'

