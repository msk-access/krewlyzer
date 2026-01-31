"""
Unit tests for FSC (Fragment Size Coverage) calculations.

Tests the Rust-backed FSC counting functions directly.
"""
import pytest
import pysam
from pathlib import Path

from krewlyzer import _core


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_count_fragments_by_bins(tmp_path):
    """Test FSC fragment counting via Rust backend.
    
    FSC 5 non-overlapping size channels:
    - Ultra-short: 65-100bp
    - Core-short: 101-149bp
    - Mono-nucleosomal: 150-220bp
    - Di-nucleosomal: 221-260bp
    - Long: 261-400bp
    """
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # 80bp - counts as ultra_short
        f.write("chr1\t100\t180\t0.5\n")
        # 170bp - counts as mono_nucl
        f.write("chr1\t200\t370\t0.5\n")
        # 300bp - counts as long
        f.write("chr1\t500\t800\t0.5\n")
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
    
    ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Verify counts
    assert ultra_shorts[0] == 1  # 80bp (65-100)
    assert core_shorts[0] == 0   # None (101-149)
    assert mono_nucls[0] == 1    # 170bp (150-220)
    assert di_nucls[0] == 0      # None (221-260)
    assert longs[0] == 1         # 300bp (261-400)
    assert totals[0] == 3        # Total


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_multiple_bins(tmp_path):
    """Test FSC counting across multiple bins with 5 channels."""
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        # Bin 1 (0-1000): 2 fragments
        f.write("chr1\t100\t220\t0.5\n")  # 120bp - core_short (101-149)
        f.write("chr1\t300\t500\t0.5\n")  # 200bp - mono_nucl (150-220)
        # Bin 2 (1000-2000): 1 fragment
        f.write("chr1\t1100\t1400\t0.5\n")  # 300bp - long (261-400)
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t1000\n")
        f.write("chr1\t1000\t2000\n")
    
    ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    # Bin 1
    assert totals[0] == 2
    assert core_shorts[0] == 1   # 120bp (101-149)
    assert mono_nucls[0] == 1    # 200bp (150-220)
    
    # Bin 2
    assert totals[1] == 1
    assert longs[1] == 1         # 300bp (261-400)


@pytest.mark.unit
@pytest.mark.rust
def test_fsc_empty_bin(tmp_path):
    """Test FSC handles empty bins correctly."""
    bed_file = tmp_path / "test.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t100\t200\t0.5\n")  # 100bp - ultra_short
    
    pysam.tabix_compress(str(bed_file), str(bed_file) + ".gz", force=True)
    pysam.tabix_index(str(bed_file) + ".gz", preset="bed", force=True)
    bedgz = str(bed_file) + ".gz"
    
    bins_file = tmp_path / "bins.bed"
    with open(bins_file, "w") as f:
        f.write("chr1\t0\t50\n")     # Empty bin
        f.write("chr1\t50\t500\n")   # Bin with fragment
    
    ultra_shorts, core_shorts, mono_nucls, di_nucls, longs, totals, gcs = _core.count_fragments_by_bins(
        bedgz, str(bins_file)
    )
    
    assert totals[0] == 0  # Empty bin
    assert totals[1] == 1  # Bin with fragment


@pytest.mark.unit
def test_aggregate_by_gene_gene_mode(tmp_path):
    """Test aggregate_by_gene with aggregate_by='gene' mode."""
    from krewlyzer.core.fsc_processor import aggregate_by_gene
    from types import SimpleNamespace
    
    # Create test fragments
    bed_file = tmp_path / "test.bed.gz"
    import gzip
    with gzip.open(bed_file, 'wt') as f:
        # 170bp fragments → mono_nucl (150-260bp)
        f.write("chr1\t100\t270\t0.5\n")  # Gene1 region
        f.write("chr1\t100\t270\t0.5\n")  # Gene1 region
        f.write("chr2\t500\t670\t0.5\n")  # Gene2 region
    
    # Create mock genes dict with regions
    genes = {
        'Gene1': [SimpleNamespace(chrom='chr1', start=0, end=500, name='Gene1_target_01')],
        'Gene2': [SimpleNamespace(chrom='chr2', start=400, end=800, name='Gene2_target_01')]
    }
    
    # Test gene mode
    output = tmp_path / "test.FSC.gene.tsv"
    aggregate_by_gene(bed_file, genes, output, aggregate_by='gene')
    
    import pandas as pd
    df = pd.read_csv(output, sep='\t')
    
    assert len(df) == 2  # 2 genes
    assert 'gene' in df.columns
    assert 'normalized_depth' in df.columns
    assert 'n_regions' in df.columns
    
    gene1 = df[df['gene'] == 'Gene1'].iloc[0]
    assert gene1['total'] == 2  # 2 fragments
    assert gene1['mono_nucl'] == 2


@pytest.mark.unit
def test_aggregate_by_gene_region_mode(tmp_path):
    """Test aggregate_by_gene with aggregate_by='region' mode."""
    from krewlyzer.core.fsc_processor import aggregate_by_gene
    from types import SimpleNamespace
    
    # Create test fragments
    bed_file = tmp_path / "test.bed.gz"
    import gzip
    with gzip.open(bed_file, 'wt') as f:
        f.write("chr1\t100\t270\t0.5\n")  # Gene1 region1
        f.write("chr1\t600\t770\t0.5\n")  # Gene1 region2
    
    # Create mock genes dict with multiple regions per gene
    genes = {
        'Gene1': [
            SimpleNamespace(chrom='chr1', start=0, end=500, name='Gene1_target_01'),
            SimpleNamespace(chrom='chr1', start=550, end=900, name='Gene1_target_02')
        ]
    }
    
    # Test region mode
    output = tmp_path / "test.FSC.regions.tsv"
    aggregate_by_gene(bed_file, genes, output, aggregate_by='region')
    
    import pandas as pd
    df = pd.read_csv(output, sep='\t')
    
    assert len(df) == 2  # 2 regions
    assert 'region_name' in df.columns
    assert 'chrom' in df.columns
    assert 'normalized_depth' in df.columns
    
    # Each region should have 1 fragment
    assert df['total'].tolist() == [1, 1]


@pytest.mark.unit
def test_filter_fsc_to_e1(tmp_path):
    """Test E1-only filtering (first exon per gene by position).
    
    E1 (first exon) serves as a proxy for promoter regions.
    Per Helzer et al. (2025), E1 has stronger cancer signal than
    whole-gene averages.
    """
    from krewlyzer.core.fsc_processor import filter_fsc_to_e1
    import pandas as pd
    
    # Create test FSC.regions.tsv with multiple exons per gene
    regions_file = tmp_path / "test.FSC.regions.tsv"
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
        'start': [1000, 500, 2000, 1500, 3000],  # Gene1 E1 at 500, Gene2 E1 at 1500
        'end': [1200, 700, 2200, 1700, 3200],
        'gene': ['Gene1', 'Gene1', 'Gene2', 'Gene2', 'Gene2'],
        'region_name': ['Gene1_ex2', 'Gene1_ex1', 'Gene2_ex2', 'Gene2_ex1', 'Gene2_ex3'],
        'normalized_depth': [100.0, 150.0, 200.0, 180.0, 120.0],
        'total': [50, 75, 100, 90, 60],
    })
    df.to_csv(regions_file, sep='\t', index=False)
    
    # Filter to E1 only
    e1_output = filter_fsc_to_e1(regions_file)
    
    # Verify output
    assert e1_output is not None
    assert e1_output.exists()
    assert 'e1only' in e1_output.name
    
    e1_df = pd.read_csv(e1_output, sep='\t')
    
    # Should have 2 genes
    assert len(e1_df) == 2
    
    # Should have E1 (first by position) for each gene
    gene1 = e1_df[e1_df['gene'] == 'Gene1'].iloc[0]
    gene2 = e1_df[e1_df['gene'] == 'Gene2'].iloc[0]
    
    assert gene1['start'] == 500   # Gene1 E1 at position 500
    assert gene2['start'] == 1500  # Gene2 E1 at position 1500


@pytest.mark.unit
def test_filter_fsc_to_e1_nonexistent_file(tmp_path):
    """Test filter_fsc_to_e1 handles missing input gracefully."""
    from krewlyzer.core.fsc_processor import filter_fsc_to_e1
    
    result = filter_fsc_to_e1(tmp_path / "nonexistent.tsv")
    assert result is None
