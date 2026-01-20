"""
Unit tests for WPS anchor filtering module.

Tests filtering WPS anchors by gene names for panel-specific analysis.
"""

import pytest
from pathlib import Path
import tempfile
import gzip

from krewlyzer.core.wps_anchor_filter import (
    filter_anchors_by_genes,
    get_bundled_wps_anchors,
)


# =============================================================================
# filter_anchors_by_genes tests
# =============================================================================

class TestFilterAnchorsByGenes:
    """Tests for WPS anchor gene filtering."""
    
    def test_filter_tss_by_gene_name(self, tmp_path):
        """TSS anchors filtered by gene name."""
        # Create mock anchors file
        anchors = tmp_path / "anchors.bed.gz"
        with gzip.open(anchors, 'wt') as f:
            f.write("1\t100\t101\tTSS|MTOR|ENST001\t.\t+\n")
            f.write("1\t200\t201\tTSS|BRCA2|ENST002\t.\t+\n")
            f.write("1\t300\t301\tTSS|TP53|ENST003\t.\t+\n")
        
        output = tmp_path / "filtered.bed.gz"
        gene_names = {"MTOR", "BRCA2"}
        
        count = filter_anchors_by_genes(anchors, gene_names, output, include_ctcf_near_genes=False)
        
        assert count == 2
        
        # Verify content
        with gzip.open(output, 'rt') as f:
            lines = f.readlines()
        assert len(lines) == 2
        assert "MTOR" in lines[0]
        assert "BRCA2" in lines[1]
    
    def test_ctcf_near_tss_included(self, tmp_path):
        """CTCF anchors near target TSS are included."""
        anchors = tmp_path / "anchors.bed.gz"
        with gzip.open(anchors, 'wt') as f:
            f.write("1\t1000\t1001\tTSS|MTOR|ENST001\t.\t+\n")
            f.write("1\t1500\t1501\tCTCF|1:1501\t.\t.\n")  # 500bp away - close
            f.write("1\t500000\t500001\tCTCF|1:500001\t.\t.\n")  # 499kb away - far
        
        output = tmp_path / "filtered.bed.gz"
        gene_names = {"MTOR"}
        
        count = filter_anchors_by_genes(
            anchors, gene_names, output, 
            include_ctcf_near_genes=True, 
            ctcf_proximity_bp=10000
        )
        
        # Should get TSS + nearby CTCF, not far CTCF
        assert count == 2
    
    def test_ctcf_excluded_when_disabled(self, tmp_path):
        """CTCF anchors excluded when include_ctcf_near_genes=False."""
        anchors = tmp_path / "anchors.bed.gz"
        with gzip.open(anchors, 'wt') as f:
            f.write("1\t1000\t1001\tTSS|MTOR|ENST001\t.\t+\n")
            f.write("1\t1500\t1501\tCTCF|1:1501\t.\t.\n")
        
        output = tmp_path / "filtered.bed.gz"
        gene_names = {"MTOR"}
        
        count = filter_anchors_by_genes(
            anchors, gene_names, output,
            include_ctcf_near_genes=False
        )
        
        assert count == 1  # TSS only
    
    def test_case_insensitive_gene_matching(self, tmp_path):
        """Gene matching is case-insensitive."""
        anchors = tmp_path / "anchors.bed.gz"
        with gzip.open(anchors, 'wt') as f:
            f.write("1\t100\t101\tTSS|mtor|ENST001\t.\t+\n")
        
        output = tmp_path / "filtered.bed.gz"
        
        count = filter_anchors_by_genes(anchors, {"MTOR"}, output, include_ctcf_near_genes=False)
        assert count == 1


# =============================================================================
# Bundled anchor file tests
# =============================================================================

class TestBundledWpsAnchors:
    """Tests for bundled panel WPS anchors."""
    
    def test_get_bundled_xs1(self):
        """xs1 bundled anchors exist."""
        path = get_bundled_wps_anchors("xs1", "GRCh37")
        assert path is not None
        assert path.exists()
        assert "xs1" in str(path)
    
    def test_get_bundled_xs2(self):
        """xs2 bundled anchors exist."""
        path = get_bundled_wps_anchors("xs2", "GRCh37")
        assert path is not None
        assert path.exists()
        assert "xs2" in str(path)
    
    def test_get_bundled_unknown_returns_none(self):
        """Unknown assay returns None."""
        assert get_bundled_wps_anchors("unknown", "GRCh37") is None
