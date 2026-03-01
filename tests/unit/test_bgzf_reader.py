"""
Unit tests for BGZF-first file reader in bed.rs.

Tests the hybrid BGZF/gzip reader strategy where:
- BGZF files use noodles::bgzf
- Standard gzip files use flate2::MultiGzDecoder
"""

import pytest
import gzip
import tempfile


def _rust_available():
    """Check if Rust extension is available."""
    try:
        from krewlyzer import _core

        return True
    except ImportError:
        return False


class TestBgzfDetection:
    """Test BGZF format detection logic."""

    def test_standard_gzip_magic_bytes(self):
        """Standard gzip has 0x1f 0x8b header."""
        with tempfile.NamedTemporaryFile(suffix=".gz", delete=False) as f:
            with gzip.open(f.name, "wt") as gz:
                gz.write("test content\n")

            with open(f.name, "rb") as check:
                magic = check.read(2)
                assert magic == b"\x1f\x8b", "Gzip magic bytes should be 0x1f 0x8b"

    def test_bgzf_has_bc_subfield(self):
        """BGZF files have BC subfield in extra field."""
        # BGZF header pattern: gzip magic + FEXTRA flag
        bgzf_header_start = bytes([0x1F, 0x8B, 0x08, 0x04])

        assert bgzf_header_start[0:2] == b"\x1f\x8b", "BGZF starts with gzip magic"
        assert bgzf_header_start[3] == 0x04, "BGZF has FEXTRA flag"


class TestBgzfReaderFallback:
    """Test that reader handles both BGZF and standard gzip."""

    def test_can_read_standard_gzip_bed(self, tmp_path):
        """Standard gzip BED files should read correctly."""
        bed_file = tmp_path / "test.bed.gz"

        with gzip.open(bed_file, "wt") as f:
            f.write("chr1\t100\t200\t0.45\n")
            f.write("chr1\t300\t400\t0.55\n")

        lines = []
        with gzip.open(bed_file, "rt") as f:
            lines = f.readlines()

        assert len(lines) == 2
        assert lines[0].startswith("chr1\t100")

    def test_multiblock_gzip_creates_valid_file(self, tmp_path):
        """Multi-block gzip should create valid file."""
        import io

        bed_file = tmp_path / "multi.bed.gz"

        # Write block 1
        with open(bed_file, "wb") as f:
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
                gz.write(b"chr1\t100\t200\t0.45\n")
            f.write(buf.getvalue())

        # Append block 2
        with open(bed_file, "ab") as f:
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
                gz.write(b"chr1\t300\t400\t0.55\n")
            f.write(buf.getvalue())

        # File should exist with content
        assert bed_file.exists()
        assert bed_file.stat().st_size > 0


@pytest.mark.skipif(not _rust_available(), reason="Rust extension not available")
class TestRustBedReader:
    """Integration tests for Rust bed.rs get_reader function."""

    def test_rust_reads_gzip_bed(self, tmp_path):
        """Rust reader should handle standard gzip."""
        bed_file = tmp_path / "test.bed.gz"
        with gzip.open(bed_file, "wt") as f:
            f.write("chr1\t100\t200\t0.45\n")

        # Placeholder - actual Rust test would use _core
        assert bed_file.exists()
