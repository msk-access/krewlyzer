"""
Integration tests for mFSD (Mutant Fragment Size Distribution).

Tests the new mFSD implementation with:
- All variant types (SNV, MNV, Insertion, Deletion, Complex)
- 4-way fragment classification (REF, ALT, NonREF, N)
- 39-column output format
- Optional distributions output
"""

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
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10000, "SN": "chr1"}]}

    if reads_config is None:
        # Default: 1 mutant (T) + 1 wildtype (A) at pos 1000
        reads_config = [
            {
                "query_sequence": "TTTTT",
                "reference_start": 1000,
                "template_length": 150,
                "name": "read_mut",
            },
            {
                "query_sequence": "AAAAA",
                "reference_start": 1000,
                "template_length": 200,
                "name": "read_wt",
            },
        ]

    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for i, cfg in enumerate(reads_config):
            a = pysam.AlignedSegment()
            a.query_name = cfg.get("name", f"read{i}")
            a.query_sequence = cfg["query_sequence"]
            a.flag = 0x41 if cfg.get("is_first_in_template", True) else 0x81  # R1 or R2
            a.reference_id = 0
            a.reference_start = cfg["reference_start"]
            a.mapping_quality = cfg.get("mapq", 60)
            a.cigar = ((0, len(cfg["query_sequence"])),)  # All match
            a.template_length = cfg["template_length"]
            outf.write(a)

    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def create_mock_vcf(path, variants=None):
    """Create VCF with configurable variants.

    variants: list of tuples (chrom, pos, ref, alt)
    Default: SNV at chr1:1001 A>T
    """
    if variants is None:
        variants = [("chr1", 1001, "A", "T")]

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
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    if result.exit_code != 0:
        print(result.stdout)
        if result.exception:
            import traceback

            traceback.print_exception(
                type(result.exception), result.exception, result.exception.__traceback__
            )

    assert result.exit_code == 0

    # Check output exists
    output_file = output_dir / "test.mFSD.tsv"
    assert output_file.exists()

    # Validate output format (46 columns: 44 original + ALT_LLR + REF_LLR)
    df = pd.read_csv(output_file, sep="\t")
    assert (
        len(df.columns) == 46
    ), f"Expected 46 columns, got {len(df.columns)}: {list(df.columns)}"

    # Check expected columns exist
    expected_cols = [
        "Chrom",
        "Pos",
        "Ref",
        "Alt",
        "VarType",
        "REF_Count",
        "ALT_Count",
        "NonREF_Count",
        "N_Count",
        "Total_Count",
        "REF_MeanSize",
        "ALT_MeanSize",
        "Delta_ALT_REF",
        "KS_ALT_REF",
        "KS_Pval_ALT_REF",
        "VAF_Proxy",
        "Error_Rate",
        "Quality_Score",
        "ALT_Confidence",
        "KS_Valid",
        "ALT_LLR",
        "REF_LLR",  # New LLR columns for duplex/panel mode
    ]
    for col in expected_cols:
        assert col in df.columns, f"Missing column: {col}"

    # Verify variant type detection
    assert df.iloc[0]["VarType"] == "SNV"


def test_mfsd_with_distributions(tmp_path):
    """Test --output-distributions flag generates distributions file."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
            "--output-distributions",
        ],
    )

    assert result.exit_code == 0

    # Check main output
    assert (output_dir / "test.mFSD.tsv").exists()

    # Check distributions file (note: extension is .distributions.tsv)
    dist_files = list(output_dir.glob("*.distributions.tsv"))
    assert (
        len(dist_files) == 1
    ), f"Expected 1 distributions file, found {len(dist_files)}"

    # Validate distributions format
    df_dist = pd.read_csv(dist_files[0], sep="\t")
    expected_dist_cols = ["Chrom", "Pos", "Ref", "Alt", "Category", "Size", "Count"]
    for col in expected_dist_cols:
        assert col in df_dist.columns, f"Missing distributions column: {col}"


# =============================================================================
# Variant Type Tests
# =============================================================================


def test_mfsd_detects_mnv(tmp_path):
    """Test MNV (multi-nucleotide variant) detection."""
    bam_file = tmp_path / "test.bam"
    # Create reads that span MNV at pos 1000-1001
    create_mock_bam(
        bam_file,
        [
            {
                "query_sequence": "GCAAA",
                "reference_start": 1000,
                "template_length": 150,
                "name": "read_mut",
            },
            {
                "query_sequence": "ATAAA",
                "reference_start": 1000,
                "template_length": 180,
                "name": "read_wt",
            },
        ],
    )

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [("chr1", 1001, "AT", "GC")])  # MNV

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep="\t")
    assert df.iloc[0]["VarType"] == "MNV"


def test_mfsd_detects_insertion(tmp_path):
    """Test insertion detection."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)  # Default reads

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [("chr1", 1001, "A", "ATG")])  # Insertion

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep="\t")
    assert df.iloc[0]["VarType"] == "INS"


def test_mfsd_detects_deletion(tmp_path):
    """Test deletion detection."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [("chr1", 1001, "ATG", "A")])  # Deletion

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep="\t")
    assert df.iloc[0]["VarType"] == "DEL"


def test_mfsd_detects_complex(tmp_path):
    """Test complex variant detection."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file, [("chr1", 1001, "ATG", "CT")])  # Complex

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep="\t")
    assert df.iloc[0]["VarType"] == "COMPLEX"


# =============================================================================
# Edge Case Tests
# =============================================================================


def test_mfsd_low_count_handling(tmp_path):
    """Test low fragment count handling (MRD scenario)."""
    bam_file = tmp_path / "test.bam"
    # Only 1 ALT fragment
    create_mock_bam(
        bam_file,
        [
            {
                "query_sequence": "TTTTT",
                "reference_start": 1000,
                "template_length": 150,
                "name": "read_mut",
            },
            {
                "query_sequence": "AAAAA",
                "reference_start": 1000,
                "template_length": 180,
                "name": "read_wt1",
            },
            {
                "query_sequence": "AAAAA",
                "reference_start": 1000,
                "template_length": 175,
                "name": "read_wt2",
            },
        ],
    )

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep="\t")

    # Only 1 ALT, should have LOW confidence
    assert df.iloc[0]["ALT_Count"] == 1
    assert df.iloc[0]["ALT_Confidence"] == "LOW"
    # KS_Valid should be FALSE (need >=2 in each group)
    assert not df.iloc[0]["KS_Valid"]


def test_mfsd_verbose_logging(tmp_path):
    """Test verbose flag enables debug output."""
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
            "--verbose",
        ],
    )

    assert result.exit_code == 0
    # Verbose mode should work without errors


def test_mfsd_llr_scoring(tmp_path):
    """Test LLR scoring for low-N scenarios (duplex/panel mode).

    LLR = Log-Likelihood Ratio comparing tumor vs healthy fragment size models.
    Positive values indicate tumor-like distribution (shorter fragments ~145bp).
    Negative values indicate healthy-like distribution (nucleosomal ~167bp).
    """
    bam_file = tmp_path / "test.bam"
    # Create reads with typical tumor-like fragment sizes (~145bp)
    create_mock_bam(
        bam_file,
        [
            {
                "query_sequence": "T" * 50,
                "reference_start": 1000,
                "template_length": 140,
                "name": "tumor1",
            },
            {
                "query_sequence": "T" * 50,
                "reference_start": 1000,
                "template_length": 145,
                "name": "tumor2",
            },
            {
                "query_sequence": "T" * 50,
                "reference_start": 1000,
                "template_length": 150,
                "name": "tumor3",
            },
            # Ref reads with healthy-like sizes (~167bp)
            {
                "query_sequence": "A" * 50,
                "reference_start": 1000,
                "template_length": 165,
                "name": "healthy1",
            },
            {
                "query_sequence": "A" * 50,
                "reference_start": 1000,
                "template_length": 170,
                "name": "healthy2",
            },
        ],
    )

    vcf_file = tmp_path / "variants.vcf"
    create_mock_vcf(vcf_file)  # SNV at chr1:1001 A>T

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(vcf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    assert result.exit_code == 0

    df = pd.read_csv(output_dir / "test.mFSD.tsv", sep="\t")

    # LLR columns should exist
    assert "ALT_LLR" in df.columns
    assert "REF_LLR" in df.columns

    alt_llr = df.iloc[0]["ALT_LLR"]
    ref_llr = df.iloc[0]["REF_LLR"]

    # ALT fragments (tumor-like, ~145bp) should have positive LLR
    # REF fragments (healthy-like, ~167bp) should have negative LLR
    assert alt_llr != "NA", "ALT_LLR should not be NA when ALT fragments exist"
    assert ref_llr != "NA", "REF_LLR should not be NA when REF fragments exist"

    # Convert to float for comparison
    alt_llr_f = float(alt_llr)
    ref_llr_f = float(ref_llr)

    # Tumor-like ALT should have higher LLR than healthy-like REF
    assert (
        alt_llr_f > ref_llr_f
    ), f"Expected ALT_LLR ({alt_llr_f}) > REF_LLR ({ref_llr_f}) for tumor vs healthy"


# =============================================================================
# MAF Input Tests (Regression tests for column-index mismatch bug)
# =============================================================================


def create_mock_maf(path, variants=None, include_consequence=False):
    """Create MAF file with configurable column layout.

    Args:
        path: Path to write MAF file
        variants: list of tuples (chrom, pos, ref, alt, variant_type)
                  Default: SNV at chr1:1001 A>T
        include_consequence: If True, include extra Consequence column at index 8
                            (cBioPortal/VEP format that caused the original bug)
    """
    if variants is None:
        variants = [("chr1", 1001, "A", "T", "SNP")]

    with open(path, "w") as f:
        # Header
        if include_consequence:
            # cBioPortal MAF with extra Consequence column
            f.write(
                "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\t"
                "Start_Position\tEnd_Position\tStrand\tConsequence\t"
                "Variant_Classification\tVariant_Type\tReference_Allele\t"
                "Tumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n"
            )
            for chrom, pos, ref, alt, vtype in variants:
                f.write(
                    f"TEST_GENE\t0\tMSKCC\tGRCh37\t{chrom}\t{pos}\t{pos}\t+\t"
                    f"missense_variant\tMissense_Mutation\t{vtype}\t{ref}\t{ref}\t{alt}\tTEST_SAMPLE\n"
                )
        else:
            # Standard GDC MAF (no Consequence column)
            f.write(
                "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\t"
                "Start_Position\tEnd_Position\tStrand\t"
                "Variant_Classification\tVariant_Type\tReference_Allele\t"
                "Tumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n"
            )
            for chrom, pos, ref, alt, vtype in variants:
                f.write(
                    f"TEST_GENE\t0\tMSKCC\tGRCh37\t{chrom}\t{pos}\t{pos}\t+\t"
                    f"Missense_Mutation\t{vtype}\t{ref}\t{ref}\t{alt}\tTEST_SAMPLE\n"
                )


def test_mfsd_maf_standard_format(tmp_path):
    """Test mFSD with standard GDC MAF format (no Consequence column).

    Verifies that header-based column parsing correctly identifies
    Reference_Allele and Tumor_Seq_Allele2 in standard MAF layout.
    """
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    maf_file = tmp_path / "variants.maf"
    create_mock_maf(maf_file, [("chr1", 1001, "A", "T", "SNP")])

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(maf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    if result.exit_code != 0:
        print(result.stdout)
        if result.exception:
            import traceback

            traceback.print_exception(
                type(result.exception), result.exception, result.exception.__traceback__
            )

    assert result.exit_code == 0

    output_file = output_dir / "test.mFSD.tsv"
    assert output_file.exists()

    df = pd.read_csv(output_file, sep="\t")

    # Key assertions: parsed alleles should be actual nucleotides, not MAF metadata
    assert df.iloc[0]["Ref"] == "A", f"Expected Ref='A', got '{df.iloc[0]['Ref']}'"
    assert df.iloc[0]["Alt"] == "T", f"Expected Alt='T', got '{df.iloc[0]['Alt']}'"
    assert (
        df.iloc[0]["VarType"] == "SNV"
    ), f"Expected VarType='SNV', got '{df.iloc[0]['VarType']}'"


def test_mfsd_maf_cbio_format(tmp_path):
    """Regression test: cBioPortal MAF with extra Consequence column.

    This is the exact format that caused the original bug where:
    - fields[10] read 'SNP' (Variant_Type) instead of Reference_Allele
    - fields[12] read 'A' (Tumor_Seq_Allele1=ref) instead of Tumor_Seq_Allele2

    With header-based parsing, this should now work correctly.
    """
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    maf_file = tmp_path / "variants.maf"
    create_mock_maf(
        maf_file, [("chr1", 1001, "A", "T", "SNP")], include_consequence=True
    )

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(maf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    if result.exit_code != 0:
        print(result.stdout)
        if result.exception:
            import traceback

            traceback.print_exception(
                type(result.exception), result.exception, result.exception.__traceback__
            )

    assert result.exit_code == 0

    output_file = output_dir / "test.mFSD.tsv"
    assert output_file.exists()

    df = pd.read_csv(output_file, sep="\t")

    # The critical regression assertions:
    # Before the fix, Ref was 'SNP' and Alt was 'A' (both wrong)
    assert (
        df.iloc[0]["Ref"] == "A"
    ), f"Regression: Ref should be 'A', not '{df.iloc[0]['Ref']}' (was 'SNP' before fix)"
    assert (
        df.iloc[0]["Alt"] == "T"
    ), f"Regression: Alt should be 'T', not '{df.iloc[0]['Alt']}' (was ref allele before fix)"
    assert (
        df.iloc[0]["VarType"] == "SNV"
    ), f"Regression: VarType should be 'SNV', not '{df.iloc[0]['VarType']}' (was 'COMPLEX' before fix)"

    # With correct parsing, not everything should be NonREF
    total = df.iloc[0]["Total_Count"]
    nonref = df.iloc[0]["NonREF_Count"]
    if total > 0:
        error_rate = nonref / total
        assert (
            error_rate < 1.0
        ), f"Regression: Error_Rate should be < 1.0, got {error_rate} (was 1.0 before fix)"


def test_mfsd_maf_invalid_alleles_warns(tmp_path):
    """Test that MAF with missing required columns raises an error.

    Verifies that the parser correctly detects and reports missing
    required columns rather than silently using wrong indices.
    """
    bam_file = tmp_path / "test.bam"
    create_mock_bam(bam_file)

    # Create MAF with missing Tumor_Seq_Allele2 column
    maf_file = tmp_path / "bad_variants.maf"
    with open(maf_file, "w") as f:
        f.write("Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\n")
        f.write("KRAS\tchr1\t1001\tA\n")

    output_dir = tmp_path / "output"

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "mfsd",
            "-i",
            str(bam_file),
            "-V",
            str(maf_file),
            "-o",
            str(output_dir),
            "-s",
            "test",
        ],
    )

    # Should fail with informative error about missing column
    assert result.exit_code != 0, "Should fail when MAF is missing Tumor_Seq_Allele2"
