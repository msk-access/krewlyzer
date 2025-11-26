import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import pysam
import os
from krewlyzer.helpers import reverse_complement, get_End_motif
from krewlyzer.fsc import _calc_fsc
from krewlyzer.wps import _calc_wps

def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("A") == "T"
    assert reverse_complement("") == ""
    assert reverse_complement("N") == "N"

def test_get_End_motif():
    # Mocking motif counting
    # end5: 5' end sequence
    # end3: 3' end sequence (already reverse complemented if needed? No, get_End_motif expects raw sequences?)
    # Wait, get_End_motif in helpers.py takes (end5, end3, motif_dict).
    # It updates motif_dict.
    
    # Initialize dict with keys because get_End_motif only updates existing keys
    d = {"AAA": 0, "CCC": 0, "GGG": 0}
    
    from krewlyzer.helpers import get_End_motif
    # Correct order: Emotif, end5, end3
    get_End_motif(d, "AAA", "CCC")
    assert d["AAA"] == 1
    assert d["CCC"] == 1
    
    get_End_motif(d, "AAA", "GGG")
    assert d["AAA"] == 2
    assert d["GGG"] == 1

def test_fsc_calculation(tmp_path):
    # Create a dummy .bed.gz file
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # chrom start end gc
        # Lengths: 100, 200, 300
        f.write("chr1\t100\t200\t0.5\n") # len 100 (short)
        f.write("chr1\t300\t500\t0.5\n") # len 200 (intermediate)
        f.write("chr1\t600\t900\t0.5\n") # len 300 (long)
        
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    
    bedgz = str(bed_file) + ".gz"
    
    # Create bins file
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
        
    output_file = tmp_path / "fsc.txt"
    
    # Run _calc_fsc
    # windows=1000, continue_n=1
    _calc_fsc(bedgz, str(bins_file), 1000, 1, str(output_file))
    
    # Check output
    df = pd.read_csv(output_file, sep="\t")
    assert len(df) == 1
    assert df.iloc[0]['region'] == "chr1:0-999"
    # Z-scores might be nan if std is 0?
    # We have 1 window.
    # short_s = [1], inter_s = [1], long_s = [1].
    # mean = 1, std = 0.
    # (1-1)/0 = nan.
    # So we expect NaNs or 0s if handled.
    # The code does: (val - mean) / std.
    # If std is 0, it raises/logs error?
    # "Error calculating z-scores: invalid value encountered in true_divide" usually.
    # Or numpy returns nan/inf.
    
    # Let's add another window to make it calculable
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
        f.write("chr1\t1000\t2000\n")
        
    # Add reads for second window
    bed_file_2 = tmp_path / "test2.bed"
    with open(bed_file_2, "w") as f:
        f.write("chr1\t100\t200\t0.5\n") # Window 1: 1 short
        f.write("chr1\t1100\t1300\t0.5\n") # Window 2: 1 intermediate
    
    pysam.tabix_compress(str(bed_file_2), str(bed_file_2) + ".gz", force=True)
    pysam.tabix_index(str(bed_file_2) + ".gz", preset="bed", force=True)
    bedgz_2 = str(bed_file_2) + ".gz"
    
    _calc_fsc(bedgz_2, str(bins_file), 1000, 1, str(output_file))
    
    df = pd.read_csv(output_file, sep="\t")
    assert len(df) == 2
    # Window 1: short=1, inter=0
    # Window 2: short=0, inter=1
    # Short: [1, 0]. Mean=0.5. Std=0.7071 (ddof=1).
    # Z1 = (1-0.5)/0.7071 = 0.7071.
    # Z2 = (0-0.5)/0.7071 = -0.7071.
    
    print(df)
    assert np.isclose(df.iloc[0]['short-fragment-zscore'], 0.7071, atol=1e-4)
    assert np.isclose(df.iloc[1]['short-fragment-zscore'], -0.7071, atol=1e-4)

def test_wps_calculation(tmp_path):
    # Mock BED
    bed_file = tmp_path / "wps.bed"
    with open(bed_file, "w") as f:
        # Region: chr1:100-200 (length 101)
        # Read 1: 50-200 (len 150). Spans [110, 140]. Overlaps [-10, 260].
        # Read 2: 120-270 (len 150). Spans [180, 210]. Overlaps [60, 330].
        f.write("chr1\t50\t200\t0.5\n")
        f.write("chr1\t120\t270\t0.5\n")
        
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    # Mock TSV
    tsv_file = tmp_path / "regions.tsv"
    with open(tsv_file, "w") as f:
        # id chrom start end strand
        f.write("region1\tchr1\t100\t200\t+\n")
        
    output_pattern = str(tmp_path / "wps.%s.tsv.gz")
    
    _calc_wps(bedgz, str(tsv_file), output_pattern, protect_input=120)
    
    # Check output
    out_file = output_pattern % "region1"
    assert os.path.exists(out_file)
    
    df = pd.read_csv(out_file, sep="\t", header=None, names=["chrom", "pos", "cov", "starts", "wps"])
    
    # Check position 150
    # Read 1: Spanning? No (150 > 140). Overlapping? Yes.
    # Read 2: Spanning? No (150 < 180). Overlapping? Yes.
    # gcount=0. total=2. WPS = -2.
    row_150 = df[df['pos'] == 150].iloc[0]
    assert row_150['wps'] == -2
    
    # Check position 115
    # Read 1: Spanning? Yes (110 <= 115 <= 140). Overlapping? Yes.
    # Read 2: Spanning? No. Overlapping? Yes (60 <= 115 <= 330).
    # gcount=1. total=2. WPS = 0.
    row_115 = df[df['pos'] == 115].iloc[0]
    assert row_115['wps'] == 0
