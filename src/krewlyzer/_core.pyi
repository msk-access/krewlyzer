# Type stubs for the Rust extension module (krewlyzer._core)
# This file tells mypy about the types in the native PyO3/Rust extension.
# Mirrors the public API exported in rust/src/lib.rs.
#
# Submodule namespaces correspond to sub-PyModules registered in lib.rs:
#   gc, mfsd, fsd, ocf, wps, region_entropy, region_mds, pon_builder,
#   extract_motif, uxm
#
# Keep in sync with Rust signatures in rust/src/*.rs when updating.
#
# NOTE: Path arguments are typed as `str` (not `pathlib.Path`) because PyO3's
# PathBuf extraction only guarantees `str` from Python's perspective when
# called via `str(path)` at all call sites.  Accept `str | None` for optional
# path args to match actual usage throughout the Python codebase.

from typing import Any

# ---------------------------------------------------------------------------
# Top-level functions
# ---------------------------------------------------------------------------

def configure_threads(num_threads: int = 0) -> None:
    """Configure the global Rayon thread pool used by all Rust compute kernels.

    Args:
        num_threads: Number of threads (0 = auto-detect CPU count).
    """
    ...

def run_unified_pipeline(
    bed_path: str,
    gc_ref_path: str | None,
    valid_regions_path: str | None,
    correction_out_path: str | None,
    correction_input_path: str | None,
    fsc_bins: str | None,
    fsc_output: str | None,
    wps_regions: str | None,
    wps_output: str | None,
    wps_background_regions: str | None,
    wps_background_output: str | None,
    wps_empty: bool,
    fsd_arms: str | None,
    fsd_output: str | None,
    ocf_regions: str | None,
    ocf_output: str | None,
    target_regions_path: str | None,
    bait_padding: int,
    output_format: str,
    compress: bool,
    silent: bool,
) -> None:
    """Single-pass Rust engine: Extract → GC-correct → FSC → WPS → FSD → OCF."""
    ...

def aggregate_by_gene(
    bed_path: str,
    gene_bed_path: str,
    output_path: str,
    gc_factors_path: str | None,
    aggregate_by: str,
) -> int:
    """Aggregate fragment counts by gene region from FSC BED output.

    Returns:
        Number of gene regions written.
    """
    ...

# ---------------------------------------------------------------------------
# krewlyzer._core.extract_motif
# ---------------------------------------------------------------------------

class extract_motif:
    @staticmethod
    def process_bam_parallel(
        bam_path: str,
        fasta_path: str,
        mapq: int,
        min_len: int,
        max_len: int,
        kmer: int,
        threads: int,
        output_bed_path: str | None,
        output_motif_prefix: str | None,
        exclude_path: str | None,
        target_regions_path: str | None,
        skip_duplicates: bool,
        require_proper_pair: bool,
        silent: bool,
    ) -> tuple[
        int,  # fragment count
        dict[str, int],  # end_motifs (off-target)
        dict[str, int],  # bp_motifs (off-target)
        dict[tuple[int, int], int],  # gc_observations (off-target)
        dict[str, int],  # end_motifs_on (on-target)
        dict[str, int],  # bp_motifs_on (on-target)
        dict[tuple[int, int], int],  # gc_observations_ontarget (on-target)
    ]:
        """Extract fragments and motifs from a BAM in parallel.

        Returns a 7-tuple of (fragment_count, end_motifs, bp_motifs,
        gc_obs, end_motifs_on, bp_motifs_on, gc_obs_on).
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.gc
# ---------------------------------------------------------------------------

class gc:
    @staticmethod
    def compute_and_write_gc_factors(
        observations: dict[tuple[int, int], int],
        gc_ref_path: str,
        valid_regions_path: str,
        output_path: str,
        output_format: str = "tsv",
        compress: bool = False,
    ) -> int:
        """Compute LOESS GC correction factors and write TSV/Parquet output.

        Returns:
            Number of correction factor bins written.
        """
        ...

    @staticmethod
    def generate_valid_regions(
        reference_path: str,
        exclude_regions_path: str,
        output_path: str,
        bin_size: int = 100000,
    ) -> int:
        """Generate valid region BED from reference FASTA.

        Args:
            reference_path: Path to indexed reference FASTA.
            exclude_regions_path: Path to exclude-regions BED (blacklist).
            output_path: Output path for the valid-regions BED.gz.
            bin_size: Bin size in bp (default: 100 000).

        Returns:
            Number of valid regions written.
        """
        ...

    @staticmethod
    def generate_valid_regions_ontarget(
        valid_regions_path: str,
        target_regions_path: str,
        output_path: str,
    ) -> int:
        """Filter valid regions to those overlapping target regions.

        Returns:
            Number of on-target valid regions written.
        """
        ...

    @staticmethod
    def generate_ref_genome_gc(
        reference_path: str,
        valid_regions_path: str,
        output_path: str,
    ) -> int:
        """Generate reference genome GC content table (Parquet).

        Returns:
            Total fragment count processed.
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.mfsd
# ---------------------------------------------------------------------------

class mfsd:
    @staticmethod
    def calculate_mfsd(
        bam_path: str,
        input_file: str,
        output_file: str,
        input_format: str,
        map_quality: int,
        min_frag_len: int,
        max_frag_len: int,
        output_distributions: bool,
        reference_path: str | None,
        correction_factors_path: str | None,
        require_proper_pair: bool,
        duplex_mode: bool,
        silent: bool = False,
        *,
        min_baseq: int = 20,
    ) -> None:
        """Compute mutant Fragment Size Distribution statistics from BAM + MAF/VCF."""
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.fsd
# ---------------------------------------------------------------------------

class fsd:
    @staticmethod
    def apply_pon_logratio(
        fsd_input_path: str,
        pon_parquet_path: str,
        output_path: str | None = None,
        baseline_table: str | None = None,
    ) -> int:
        """Apply PON log-ratio normalization to FSD TSV output.

        Returns:
            Number of arms processed.
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.ocf
# ---------------------------------------------------------------------------

class ocf:
    @staticmethod
    def apply_pon_zscore(
        ocf_path: str,
        pon_parquet_path: str,
        output_path: str | None = None,
    ) -> int:
        """Apply PON z-score normalization to OCF TSV output.

        Returns:
            Number of tissue types processed.
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.wps
# ---------------------------------------------------------------------------

class wps:
    @staticmethod
    def apply_pon_zscore(
        wps_parquet_path: str,
        pon_parquet_path: str,
        output_path: str | None = None,
    ) -> int:
        """Apply PON z-score normalization to WPS Parquet output.

        Returns:
            Number of regions processed.
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.region_entropy
# ---------------------------------------------------------------------------

class region_entropy:
    @staticmethod
    def run_region_entropy(
        bed_path: str,
        region_path: str,
        output_path: str,
        gc_correction_path: str | None,
        gc_correction_ontarget_path: str | None,
        target_regions_path: str | None = None,
        silent: bool = False,
    ) -> tuple[int, int]:
        """Compute per-region entropy from fragment BED file.

        Returns:
            (n_off_target_regions, n_on_target_regions)
        """
        ...

    @staticmethod
    def apply_pon_zscore(
        entropy_path: str,
        pon_parquet_path: str,
        output_path: str,
        baseline_table: str = "entropy_baseline",
    ) -> int:
        """Apply PON z-score normalization to region entropy output.

        Returns:
            Number of regions processed.
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.region_mds
# ---------------------------------------------------------------------------

class region_mds:
    @staticmethod
    def run_region_mds(
        bam_path: str,
        fasta_path: str,
        gene_bed_path: str,
        output_exon_path: str,
        output_gene_path: str,
        e1_only: bool,
        mapq: int,
        min_len: int,
        max_len: int,
        silent: bool,
    ) -> tuple[int, int]:
        """Compute per-region Motif Diversity Score.

        Returns:
            (n_exon_regions, n_gene_regions)
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.pon_builder
# ---------------------------------------------------------------------------

class pon_builder:
    @staticmethod
    def compute_gc_bias_model(
        all_gc_data: list[Any],
    ) -> Any:
        """Compute GC bias model from cohort GC observation data.

        Returns:
            Dict with GC model parameters.
        """
        ...

    @staticmethod
    def compute_fsd_baseline(
        fsd_paths: list[str],
    ) -> Any:
        """Compute FSD arm-length baseline from a cohort of FSD TSV files.

        Returns:
            Dict with per-arm mean/std baseline values.
        """
        ...

    @staticmethod
    def compute_wps_baseline(
        wps_paths: list[str],
    ) -> Any:
        """Compute WPS region baseline from a cohort of WPS Parquet files.

        Returns:
            Dict with per-region mean/std baseline values.
        """
        ...

    @staticmethod
    def compute_region_mds_baseline(
        mds_paths: list[str],
    ) -> Any:
        """Compute region MDS baseline from a cohort of MDS output files.

        Returns:
            Dict with per-gene mean/std baseline values.
        """
        ...

# ---------------------------------------------------------------------------
# krewlyzer._core.uxm
# ---------------------------------------------------------------------------

class uxm:
    @staticmethod
    def calculate_uxm(
        bam_path: str,
        marker_path: str,
        output_file: str,
        map_quality: int,
        min_cpg: int,
        methy_threshold: float,
        unmethy_threshold: float,
        pe_type: str,
        silent: bool = False,
    ) -> None:
        """Compute UXM (UnMethylated/eXtreme Methylation) scores from BAM."""
        ...
