import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import pysam
import os
from krewlyzer.motif import motif_process
from unittest.mock import patch, MagicMock

def test_motif_blacklist_filtering(tmp_path):
    # Create dummy BAM
    bam_file = tmp_path / "test.bam"
    header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}] }
    with pysam.AlignmentFile(bam_file, "wb", header=header) as outf:
        # Read 1: chr1:100-150 (overlaps blacklist)
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "A" * 50
        a.flag = 99 # paired, mapped, proper, first
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigar = ((0, 50),)
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 150
        outf.write(a)
        
        b = pysam.AlignedSegment()
        b.query_name = "read1"
        b.query_sequence = "T" * 50
        b.flag = 147 # paired, mapped, proper, second
        b.reference_id = 0
        b.reference_start = 200
        b.mapping_quality = 60
        b.cigar = ((0, 50),)
        b.next_reference_id = 0
        b.next_reference_start = 100
        b.template_length = -150
        outf.write(b)
        
        # Read 2: chr1:500-550 (no overlap)
        c = pysam.AlignedSegment()
        c.query_name = "read2"
        c.query_sequence = "C" * 50
        c.flag = 99
        c.reference_id = 0
        c.reference_start = 500
        c.mapping_quality = 60
        c.cigar = ((0, 50),)
        c.next_reference_id = 0
        c.next_reference_start = 600
        c.template_length = 150
        outf.write(c)
        
        d = pysam.AlignedSegment()
        d.query_name = "read2"
        d.query_sequence = "G" * 50
        d.flag = 147
        d.reference_id = 0
        d.reference_start = 600
        d.mapping_quality = 60
        d.cigar = ((0, 50),)
        d.next_reference_id = 0
        d.next_reference_start = 500
        d.template_length = -150
        outf.write(d)
        
    pysam.index(str(bam_file))
    
    # Create Blacklist
    # Overlaps Read 1 (100-250)
    # Blacklist: chr1:120-130
    blacklist_file = tmp_path / "blacklist.bed"
    with open(blacklist_file, "w") as f:
        f.write("chr1\t120\t130\n")
        
    # Create Genome FASTA
    genome_file = tmp_path / "genome.fa"
    with open(genome_file, "w") as f:
        f.write(">chr1\n")
        f.write("N" * 1000 + "\n")
    pysam.faidx(str(genome_file))
    
    output_bed = tmp_path / "output.bed"
    
    # Run motif_process
    motif_process(
        bamInput=str(bam_file),
        blacklistInput=str(blacklist_file),
        bedOutput=str(output_bed),
        genome_reference=str(genome_file),
        CHR=["chr1"],
        mapQuality=20,
        k_mer=3,
        fragFilter=False
    )
    
    # Check output
    # Read 1 should be filtered out.
    # Read 2 should remain.
    
    # Output is .bed.gz
    out_gz = str(output_bed) + ".gz"
    assert os.path.exists(out_gz)
    
    with pysam.TabixFile(out_gz) as tbx:
        rows = list(tbx.fetch("chr1", parser=pysam.asTuple()))
        assert len(rows) == 1
        # Read 2: start 500, end 650
        assert int(rows[0][1]) == 500
        assert int(rows[0][2]) == 650
