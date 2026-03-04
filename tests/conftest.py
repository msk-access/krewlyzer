"""
Shared pytest fixtures for krewlyzer test suite.

This file provides common fixtures used across unit, integration, and e2e tests.

DATA AVAILABILITY NOTE:
The entire src/krewlyzer/data/ folder is EXCLUDED from PyPI wheels to keep size <100MB.
Tests that require bundled data files use DATA_AVAILABLE to skip in PyPI installs.
For full test coverage, use: git clone + pip install -e .
"""

import pytest
from pathlib import Path
import pysam
import gzip

# =============================================================================
# Data Availability Detection
# =============================================================================


def _check_data_available():
    """
    Check if bundled data is available in the installed package.

    In CI: git checkout has data, but pip install . creates wheel without data.
    Tests run from checkout but imports come from installed package.
    So we must check the INSTALLED package path, not source.
    """
    try:
        import krewlyzer

        pkg_path = Path(krewlyzer.__file__).parent
        data_dir = pkg_path / "data"
        # Check if data directory and at least one key asset exist
        return data_dir.exists() and (data_dir / "genes").exists()
    except ImportError:
        return False


DATA_AVAILABLE = _check_data_available()

# Marker for tests that require bundled data
requires_data = pytest.mark.skipif(
    not DATA_AVAILABLE,
    reason="Bundled data not available (PyPI install). Use git clone + pip install -e .",
)


# =============================================================================
# Core File Fixtures
# =============================================================================


@pytest.fixture
def temp_dir(tmp_path):
    """Provide a temporary directory for test outputs."""
    yield tmp_path
    # Cleanup handled by pytest tmp_path


@pytest.fixture
def sample_name():
    """Standard sample name for tests."""
    return "test_sample"


@pytest.fixture
def temp_reference(tmp_path):
    """Create minimal FASTA reference with chr1."""
    ref = tmp_path / "genome.fa"
    seq = "ACGTACGTACGT" * 1000  # 12kb
    ref.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(ref))
    return ref


@pytest.fixture
def temp_bam(tmp_path):
    """Create minimal valid BAM with a single proper pair."""
    bam = tmp_path / "test.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10000, "SN": "chr1"}]}

    with pysam.AlignmentFile(str(bam), "wb", header=header) as outf:
        # Read 1 (R1)
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

        # Read 2 (R2)
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

    pysam.index(str(bam))
    return bam


@pytest.fixture
def temp_bedgz(tmp_path):
    """Create minimal BED.gz with a few fragments."""
    bed = tmp_path / "test.bed.gz"

    fragments = [
        ("chr1", 100, 250, 0.45),  # 150bp
        ("chr1", 500, 660, 0.50),  # 160bp
        ("chr1", 1000, 1170, 0.48),  # 170bp
    ]

    with gzip.open(bed, "wt") as f:
        for chrom, start, end, gc in fragments:
            f.write(f"{chrom}\t{start}\t{end}\t{gc:.4f}\n")

    pysam.tabix_index(str(bed), preset="bed", force=True)
    return bed


# =============================================================================
# Resource File Fixtures
# =============================================================================


@pytest.fixture
def temp_bins(tmp_path):
    """Create bins file for FSC testing."""
    bins = tmp_path / "bins.bed"
    bins.write_text("chr1\t0\t5000\tBin1\n")
    return bins


@pytest.fixture
def temp_arms(tmp_path):
    """Create arms file for FSD testing."""
    arms = tmp_path / "arms.bed"
    arms.write_text("chr1\t0\t10000\t1p\n")
    return arms


@pytest.fixture
def temp_transcripts(tmp_path):
    """Create transcript annotation for WPS testing."""
    tsv = tmp_path / "transcripts.tsv"
    tsv.write_text("Gene1\tchr1\t1000\t2000\t+\n")
    return tsv


@pytest.fixture
def temp_ocr(tmp_path):
    """Create OCR file for OCF testing."""
    ocr = tmp_path / "ocr.bed"
    ocr.write_text("chr1\t100\t200\tTissueA\n")
    return ocr


@pytest.fixture
def temp_markers(tmp_path):
    """Create markers file for UXM testing."""
    markers = tmp_path / "markers.bed"
    markers.write_text("chr1\t1000\t1005\n")
    return markers


@pytest.fixture
def temp_vcf(tmp_path):
    """Create VCF file for mFSD testing."""
    vcf = tmp_path / "variants.vcf"
    vcf.write_text("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1001\t.\tA\tT\t.\tPASS\t.
""")
    return vcf


# =============================================================================
# Complex Fixtures (Multiple Files)
# =============================================================================


@pytest.fixture
def full_test_data(
    tmp_path, temp_bam, temp_reference, temp_bins, temp_arms, temp_transcripts, temp_ocr
):
    """Bundle of all resource files needed for run-all testing."""
    return {
        "bam": temp_bam,
        "reference": temp_reference,
        "bins": temp_bins,
        "arms": temp_arms,
        "transcripts": temp_transcripts,
        "ocr": temp_ocr,
        "output_dir": tmp_path / "output",
    }


# =============================================================================
# Mock BAM Creators (for specific test scenarios)
# =============================================================================


def create_mock_bam_with_reads(path, reads_config):
    """
    Create BAM with configurable reads.

    Args:
        path: Path to write BAM
        reads_config: list of dicts with:
            - query_name: str
            - query_sequence: str
            - reference_start: int (0-based)
            - template_length: int
            - mapping_quality: int (default 60)
            - flag: int (default 99 for R1)
    """
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10000, "SN": "chr1"}]}

    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for read in reads_config:
            a = pysam.AlignedSegment()
            a.query_name = read.get("query_name", "read1")
            a.query_sequence = read.get("query_sequence", "A" * 50)
            a.flag = read.get("flag", 99)
            a.reference_id = 0
            a.reference_start = read.get("reference_start", 100)
            a.mapping_quality = read.get("mapping_quality", 60)
            a.cigar = ((0, len(a.query_sequence)),)
            a.next_reference_id = 0
            a.next_reference_start = read.get("mate_start", 200)
            a.template_length = read.get("template_length", 150)

            # Optional tags
            if "xm_tag" in read:
                a.set_tag("XM", read["xm_tag"])

            outf.write(a)

    pysam.index(str(path))
    return path


# =============================================================================
# Real Test Fixtures (from actual data)
# =============================================================================

FIXTURES_DIR = Path(__file__).parent / "data" / "fixtures"


@pytest.fixture
def real_bam():
    """Real BAM file subset (3144 reads from chr1:1-2000000)."""
    bam = FIXTURES_DIR / "test_sample.bam"
    if bam.exists():
        return bam
    pytest.skip("Real BAM fixture not available")


@pytest.fixture
def real_reference():
    """Real FASTA reference (2Mb chr1)."""
    ref = FIXTURES_DIR / "test_genome.fa"
    if ref.exists():
        return ref
    pytest.skip("Real reference fixture not available")


@pytest.fixture
def real_bins():
    """Real bins file (20 x 100kb bins on chr1)."""
    bins = FIXTURES_DIR / "test_bins.bed"
    if bins.exists():
        return bins
    pytest.skip("Real bins fixture not available")


@pytest.fixture
def real_arms():
    """Real arms file (chr1p, chr1q)."""
    arms = FIXTURES_DIR / "test_arms.bed"
    if arms.exists():
        return arms
    pytest.skip("Real arms fixture not available")


@pytest.fixture
def real_transcripts():
    """Real transcripts file (3 test genes)."""
    tsv = FIXTURES_DIR / "test_transcripts.tsv"
    if tsv.exists():
        return tsv
    pytest.skip("Real transcripts fixture not available")


@pytest.fixture
def real_pon():
    """Real PON model (msk-access-v1)."""
    pon = FIXTURES_DIR / "test.pon.parquet"
    if pon.exists():
        return pon
    pytest.skip("Real PON fixture not available")
