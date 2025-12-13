"""
Test extract blacklist/exclude regions filtering.
"""
import pytest
from pathlib import Path
import pysam
from typer.testing import CliRunner
from krewlyzer.cli import app


def test_extract_exclude_regions(tmp_path):
    # Create dummy BAM
    bam_file = tmp_path / "test.bam"
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
    with pysam.AlignmentFile(bam_file, "wb", header=header) as outf:
        # Read 1: chr1:100-250 (overlaps exclude at 120-130)
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "A" * 50
        a.flag = 99  # paired, mapped, proper, first
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
        b.flag = 147  # paired, mapped, proper, second
        b.reference_id = 0
        b.reference_start = 200
        b.mapping_quality = 60
        b.cigar = ((0, 50),)
        b.next_reference_id = 0
        b.next_reference_start = 100
        b.template_length = -150
        outf.write(b)
        
        # Read 2: chr1:500-650 (no overlap with exclude)
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
    
    # Create exclude regions (overlaps Read 1)
    exclude_file = tmp_path / "exclude.bed"
    with open(exclude_file, "w") as f:
        f.write("chr1\t120\t130\n")
        
    # Create Genome FASTA
    genome_file = tmp_path / "genome.fa"
    with open(genome_file, "w") as f:
        f.write(">chr1\n")
        f.write("N" * 1000 + "\n")
    pysam.faidx(str(genome_file))
    
    output_dir = tmp_path / "output"
    
    # Run extract via CLI
    runner = CliRunner()
    result = runner.invoke(app, [
        "extract", str(bam_file),
        "-g", str(genome_file),
        "-o", str(output_dir),
        "--exclude-regions", str(exclude_file),
        "--chromosomes", "chr1"
    ])
    
    assert result.exit_code == 0, f"CLI failed: {result.output}"
    
    # Check output
    out_gz = output_dir / "test.bed.gz"
    assert out_gz.exists()
    
    # Read 1 should be filtered (overlaps exclude)
    # Read 2 should remain
    with pysam.TabixFile(str(out_gz)) as tbx:
        rows = list(tbx.fetch("chr1", parser=pysam.asTuple()))
        # Only read2 should be present
        assert len(rows) == 1
        assert int(rows[0][1]) == 500  # start
