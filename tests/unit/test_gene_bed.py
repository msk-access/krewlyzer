"""
Unit tests for gene_bed module.

Tests gene BED parsing, version detection, and bundled file loading
for MSK-ACCESS panel assays (xs1, xs2).
"""

import pytest
from pathlib import Path
import gzip

from conftest import requires_data

from krewlyzer.core.gene_bed import (
    detect_assay,
    detect_version_from_content,
    parse_gene_bed,
    load_gene_bed,
    get_bundled_gene_bed,
    _is_gene_region,
    _extract_gene_v1,
    _extract_gene_v2,
)

# =============================================================================
# detect_assay tests
# =============================================================================


class TestDetectAssay:
    """Tests for assay detection from filenames."""

    def test_detect_xs1_from_filename(self):
        assert detect_assay(Path("xs1.targets.bed")) == "xs1"
        assert detect_assay(Path("xs1.genes.bed.gz")) == "xs1"

    def test_detect_xs2_from_filename(self):
        assert detect_assay(Path("xs2.targets.bed")) == "xs2"
        assert detect_assay(Path("data/GRCh37/xs2.genes.bed.gz")) == "xs2"

    def test_detect_legacy_v1_patterns(self):
        assert detect_assay(Path("MSK-ACCESS-v1_0-probe-A.sorted.bed")) == "xs1"

    def test_detect_legacy_v2_patterns(self):
        assert detect_assay(Path("MSK-ACCESS-v2_targetsAllwFP.bed")) == "xs2"

    def test_unknown_assay_returns_none(self):
        assert detect_assay(Path("custom_panel.bed")) is None
        assert detect_assay(Path("random_file.bed")) is None


# =============================================================================
# Version detection tests
# =============================================================================


class TestDetectVersionFromContent:
    """Tests for version detection from BED file content."""

    def test_detect_v1_format(self, tmp_path):
        bed = tmp_path / "v1.bed"
        bed.write_text(
            "1\t100\t200\t-\texon_MTOR_48a.1_2\n1\t200\t300\t-\texon_NRAS_3.1_1\n"
        )
        assert detect_version_from_content(bed) == "V1"

    def test_detect_v2_format(self, tmp_path):
        bed = tmp_path / "v2.bed"
        bed.write_text(
            "1\t100\t200\tMTOR_target_02\t0\t.\n1\t200\t300\tNRAS_target_03\t0\t.\n"
        )
        assert detect_version_from_content(bed) == "V2"

    def test_detect_gene_bed_format(self, tmp_path):
        bed = tmp_path / "gene.bed"
        bed.write_text(
            "#chrom\tstart\tend\tgene\tname\n1\t100\t200\tMTOR\tMTOR_target_02\n"
        )
        assert detect_version_from_content(bed) == "GENE_BED"

    def test_detect_v1_gzipped(self, tmp_path):
        bed = tmp_path / "v1.bed.gz"
        with gzip.open(bed, "wt") as f:
            f.write("1\t100\t200\t-\texon_MTOR_48a.1_2\n")
        assert detect_version_from_content(bed) == "V1"


# =============================================================================
# Gene region filtering tests
# =============================================================================


class TestIsGeneRegion:
    """Tests for filtering non-gene regions."""

    def test_gene_regions_accepted(self):
        assert _is_gene_region("exon_MTOR_48a.1_2") is True
        assert _is_gene_region("BRCA2_target_01") is True

    def test_fingerprint_filtered(self):
        assert _is_gene_region("FP_rs2051068") is False
        assert _is_gene_region("snp_FP_rs123") is False

    def test_msi_filtered(self):
        assert _is_gene_region("msi_Target002") is False

    def test_bat_filtered(self):
        assert _is_gene_region("BAT-25") is False


# =============================================================================
# Gene extraction tests
# =============================================================================


class TestExtractGene:
    """Tests for gene name extraction from region names."""

    def test_extract_v1_standard(self):
        assert _extract_gene_v1("exon_MTOR_48a.1_2") == "MTOR"
        assert _extract_gene_v1("exon_ARID1A_1a.1_1") == "ARID1A"
        assert _extract_gene_v1("exon_H3F3A_2_1") == "H3F3A"

    def test_extract_v2_standard(self):
        assert _extract_gene_v2("MTOR_target_02") == "MTOR"
        assert _extract_gene_v2("BRCA2_target_01") == "BRCA2"

    def test_extract_returns_none_for_invalid(self):
        assert _extract_gene_v1("MTOR_target_02") is None  # V2 pattern
        assert _extract_gene_v2("exon_MTOR_48a.1_2") is None  # V1 pattern


# =============================================================================
# parse_gene_bed tests
# =============================================================================


class TestParseGeneBed:
    """Tests for parsing gene BED files."""

    def test_parse_v2_format(self, tmp_path):
        bed = tmp_path / "v2.bed"
        bed.write_text(
            "1\t100\t200\tMTOR_target_01\t0\t.\n"
            "1\t200\t300\tMTOR_target_02\t0\t.\n"
            "1\t400\t500\tBRCA2_target_01\t0\t.\n"
        )

        genes = parse_gene_bed(bed)

        assert len(genes) == 2
        assert "MTOR" in genes
        assert "BRCA2" in genes
        assert len(genes["MTOR"]) == 2
        assert len(genes["BRCA2"]) == 1

    def test_parse_filters_msi_regions(self, tmp_path):
        bed = tmp_path / "v2.bed"
        bed.write_text(
            "1\t100\t200\tMTOR_target_01\t0\t.\n" "1\t200\t300\tmsi_Target002\t0\t.\n"
        )

        genes = parse_gene_bed(bed)

        assert len(genes) == 1
        assert "MTOR" in genes

    def test_parse_gene_bed_format(self, tmp_path):
        bed = tmp_path / "gene.bed"
        bed.write_text(
            "#chrom\tstart\tend\tgene\tname\n"
            "1\t100\t200\tMTOR\tMTOR_target_01\n"
            "1\t200\t300\tMTOR\tMTOR_target_02\n"
        )

        genes = parse_gene_bed(bed)

        assert len(genes) == 1
        assert "MTOR" in genes
        assert len(genes["MTOR"]) == 2


# =============================================================================
# Bundled file tests
# =============================================================================


@requires_data
class TestBundledGeneBed:
    """Tests for bundled gene BED file loading."""

    def test_get_bundled_xs1(self):
        path = get_bundled_gene_bed("xs1", "GRCh37")
        assert path is not None
        assert path.exists()
        assert "xs1" in str(path)

    def test_get_bundled_xs2(self):
        path = get_bundled_gene_bed("xs2", "GRCh37")
        assert path is not None
        assert path.exists()
        assert "xs2" in str(path)

    def test_get_bundled_unknown_returns_none(self):
        assert get_bundled_gene_bed("unknown", "GRCh37") is None


@requires_data
class TestLoadGeneBed:
    """Tests for the unified load_gene_bed function."""

    def test_load_by_assay(self):
        genes = load_gene_bed(assay="xs2", genome="GRCh37")
        assert len(genes) > 100  # xs2 has 146 genes
        assert "MTOR" in genes
        assert "ATM" in genes

    def test_load_xs1_by_assay(self):
        genes = load_gene_bed(assay="xs1", genome="GRCh37")
        assert len(genes) > 100  # xs1 has 128 genes

    def test_load_raises_on_no_source(self):
        with pytest.raises(ValueError):
            load_gene_bed()  # No assay, no target_bed, no gene_bed
