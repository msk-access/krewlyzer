"""
Integration tests for mFSD (Mutant Fragment Size Distribution).

Tests the new mFSD implementation with:
- All variant types (SNV, MNV, Insertion, Deletion, Complex)
- 4-way fragment classification (REF, ALT, NonREF, N)
- 39-column output format
- Optional distributions output
"""
import pytest
from pathlib import Path
import pysam
import pandas as pd
from krewlyzer.cli import app
from typer.testing import CliRunner


def create_mock_bam(path, reads_config=None):
    """Create BAM with configurable reads for testing different scenarios.
    
    reads_config: list of dicts with:
        - query_sequence: str
        - reference_start: int (0-based)
        - template_length: int
        - is_first_in_template: bool (default True for R1)
    """
    header = { 
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 10000, 'SN': 'chr1'}] 
    }
    
    if reads_config is None:
        # Default: 1 mutant (T) + 1 wildtype (A) at pos 1000
        reads_config = [
            {'query_sequence': 'TTTTT', 'reference_start': 1000, 'template_length': 150, 'name': 'read_mut'},
            {'query_sequence': 'AAAAA', 'reference_start': 1000, 'template_length': 200, 'name': 'read_wt'},
        ]
    
    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for i, cfg in enumerate(reads_config):
            a = pysam.AlignedSegment()
            a.query_name = cfg.get('name', f"read{i}")
            a.query_sequence = cfg['query_sequence']
            a.flag = 0x41 if cfg.get('is_first_in_template', True) else 0x81  # R1 or R2
            a.reference_id = 0
            a.reference_start = cfg['reference_start']
            a.mapping_quality = cfg.get('mapq', 60)
            a.cigar = ((0, len(cfg['query_sequence'])),)  # All match
            a.template_length = cfg['template_length']
            outf.write(a)

    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def create_mock_vcf(path, variants=None):
    """Create VCF with configurable variants.
    
    variants: list of tuples (chrom, pos, ref, alt)
    Default: SNV at chr1:1001 A>T
    """
    if variants is None:
        variants = [('chr1', 1001, 'A', 'T')]
    
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in variants:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")


# =============================================================================
# Basic Integration Tests
# =============================================================================

def test_mfsd_integration_snv(tmp_path):
    """Test basic SNV classification and output format."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    if result.exit_code != 0:
        print(result.stdout)
        if result.exception:
            import traceback
            traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)
            
    assert result.exit_code == 0
    
    # Check output exists
    output_file = output_dir / "test.mFSD.tsv"
    assert output_file.exists()
    
    # Validate output format (39 columns)
    df = pd.read_csv(output_file, sep='\t')
    assert len(df.columns) == 39, f"Expected 39 columns, got {len(df.columns)}"
    
    # Check expected columns exist
    expected_cols = [
        'Chrom', 'Pos', 'Ref', 'Alt', 'VarType',
        'REF_Count', 'ALT_Count', 'NonREF_Count', 'N_Count', 'Total_Count',
        'REF_MeanSize', 'ALT_MeanSize',
        'Delta_ALT_REF', 'KS_ALT_REF', 'KS_Pval_ALT_REF',
        'VAF_Proxy', 'Error_Rate', 'Quality_Score',
        'ALT_Confidence', 'KS_Valid'
    ]
    for col in expected_cols:
        assert col in df.columns, f"Missing column: {col}"
    
    # Verify variant type detection
    assert df.iloc[0]['VarType'] == 'SNV'


def test_mfsd_with_distributions(tmp_path):
    """Test --output-distributions flag generates distributions file."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test",
        "--output-distributions"
    ])
    
    assert result.exit_code == 0
    
    # Check main output
    assert (output_dir / "test.mFSD.tsv").exists()
    
    # Check distributions file (note: extension is .distributions.tsv)
    dist_files = list(output_dir.glob("*.distributions.tsv"))
    assert len(dist_files) == 1, f"Expected 1 distributions file, found {len(dist_files)}"
    
    # Validate distributions format
    df_dist = pd.read_csv(dist_files[0], sep='\t')
    expected_dist_cols = ['Chrom', 'Pos', 'Ref', 'Alt', 'Category', 'Size', 'Count']
    for col in expected_dist_cols:
        assert col in df_dist.columns, f"Missing distributions column: {col}"


# =============================================================================
# Variant Type Tests
# =============================================================================

def test_mfsd_detects_mnv(tmp_path):
    """Test MNV (multi-nucleotide variant) detection."""
    bam_file = tmp_path / "test.bam"
    # Create reads that span MNV at pos 1000-1001
    create_mock_bam(bam_file, [
        {'query_sequence': 'GCAAA', 'reference_start': 1000, 'template_length': 150, 'name': 'read_mut'},
        {'query_sequence': 'ATAAA', 'reference_start': 1000, 'template_length': 180, 'name': 'read_wt'},
    ])
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [('chr1', 1001, 'AT', 'GC')])  # MNV
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    assert result.exit_code == 0
    
    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep='\t')
    assert df.iloc[0]['VarType'] == 'MNV'


def test_mfsd_detects_insertion(tmp_path):
    """Test insertion detection."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)  # Default reads
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [('chr1', 1001, 'A', 'ATG')])  # Insertion
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    assert result.exit_code == 0
    
    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep='\t')
    assert df.iloc[0]['VarType'] == 'INS'


def test_mfsd_detects_deletion(tmp_path):
    """Test deletion detection."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [('chr1', 1001, 'ATG', 'A')])  # Deletion
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    assert result.exit_code == 0
    
    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep='\t')
    assert df.iloc[0]['VarType'] == 'DEL'


def test_mfsd_detects_complex(tmp_path):
    """Test complex variant detection."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [('chr1', 1001, 'ATG', 'CT')])  # Complex
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    assert result.exit_code == 0
    
    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep='\t')
    assert df.iloc[0]['VarType'] == 'COMPLEX'


# =============================================================================
# Edge Case Tests
# =============================================================================

def test_mfsd_low_count_handling(tmp_path):
    """Test low fragment count handling (MRD scenario)."""
    bam_file = tmp_path / "test.bam"
    # Only 1 ALT fragment
    create_mock_bam(bam_file, [
        {'query_sequence': 'TTTTT', 'reference_start': 1000, 'template_length': 150, 'name': 'read_mut'},
        {'query_sequence': 'AAAAA', 'reference_start': 1000, 'template_length': 180, 'name': 'read_wt1'},
        {'query_sequence': 'AAAAA', 'reference_start': 1000, 'template_length': 175, 'name': 'read_wt2'},
    ])
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test"
    ])
    
    assert result.exit_code == 0
    
    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep='\t')
    
    # Only 1 ALT, should have LOW confidence
    assert df.iloc[0]['ALT_Count'] == 1
    assert df.iloc[0]['ALT_Confidence'] == 'LOW'
    # KS_Valid should be FALSE (need >=2 in each group)
    assert df.iloc[0]['KS_Valid'] == False


def test_mfsd_verbose_logging(tmp_path):
    """Test verbose flag enables debug output."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)
    
    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)
    
    output_dir = tmp_path / "output"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "mfsd", str(bam_file), 
        "-i", str(vcf_file), 
        "-o", str(output_dir),
        "-s", "test",
        "--verbose"
    ])
    
    assert result.exit_code == 0
    # Verbose mode should work without errors
