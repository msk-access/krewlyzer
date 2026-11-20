"""
PON (Panel of Normals) integration module.

Provides shared functionality for loading PON models and computing
background/periodicity z-scores used by the WPS processor.

Note: Tool-specific z-score functions (FSC, FSD, WPS, OCF, TFBS, ATAC) have
been migrated to Rust for performance and are called via _core.<module>.
See the following Rust entry points:
    - _core.fsd.apply_pon_logratio      → fsd_processor.py
    - _core.wps.apply_pon_zscore        → wps_processor.py
    - _core.region_entropy.apply_pon_zscore → region_entropy_processor.py
    - _core.ocf.apply_pon_zscore        → ocf_processor.py
"""

from pathlib import Path
from typing import Optional
import logging

from krewlyzer.pon.model import GENOME_WIDE_GROUP

logger = logging.getLogger("core.pon_integration")


class UnusablePonError(ValueError):
    """A PON the run was asked to score against cannot be used: unreadable, or refused.

    An error, not a warning: carrying on gave an unscored run recorded as
    ``pon_applied=False`` -- indistinguishable from a deliberate ``--skip-pon`` -- and
    ``validate-output`` then stopped requiring every z-score column. A ``ValueError``, like
    every other bad ``--pon-model`` argument, so each tool reports it as an error line.
    """


def load_pon_model(pon_path: Path):
    """
    Load a PON model from disk.

    Refuses a model older than `MIN_PON_VERSION`. 0.9.0 changed what the
    features mean, so an earlier PON is not stale, it is measuring something
    else -- and every z-score computed against it is wrong in a way no schema
    check can see. Recording the version was only half of that guard; this is
    the half that acts on it.

    Refusing rather than warning, because the alternative is the failure this
    whole release exists to remove: a plausible number nobody questions. A
    warning in a log nobody reads is not a guard.

    Args:
        pon_path: Path to the PON parquet file

    Returns:
        The loaded PonModel.

    Raises:
        UnusablePonError: the file cannot be loaded, or the model is too old to score against.
    """
    from krewlyzer.pon.model import PonModel
    from krewlyzer.pon.provenance import check_pon_version

    try:
        pon = PonModel.load(pon_path)
    except Exception as e:
        raise UnusablePonError(
            f"PON {Path(pon_path).name} cannot be loaded ({type(e).__name__}: {e}); "
            "refusing to run unscored. Pass a readable PON, or --skip-pon to run without one."
        ) from e

    complaint = check_pon_version(pon.krewlyzer_version, Path(pon_path).name)
    if complaint:
        raise UnusablePonError(complaint)

    logger.info(
        f"Loaded PON model: {pon.assay} (n={pon.n_samples}, "
        f"built for {pon.krewlyzer_version or 'unknown'})"
    )
    return pon


def load_requested_pon(
    pon_path: Optional[Path],
    record_dir: Optional[Path] = None,
    sample: Optional[str] = None,
    output_format: str = "tsv",
    compress: bool = False,
):
    """The PON a run resolved, loaded **before any work**; None when no PON was resolved.

    Every caller loads here, right after `resolve_pon_model`, so an unusable PON stops the run
    before a single table is written. Loading after the pass left unscored tables under their
    final names, beside whatever completion marker the directory already had.

    A single tool passes ``record_dir`` and ``sample``: if that sample's record says it was
    scored against another PON, scoring it against this one is refused (run-all, which rewrites
    the whole record, passes neither).
    """
    if not pon_path:
        return None
    pon = load_pon_model(pon_path)
    if record_dir is not None and sample:
        from .output_provenance import refuse_a_second_pon

        refuse_a_second_pon(
            Path(record_dir), sample, pon, Path(pon_path), output_format, compress
        )
    return pon


def compute_nrl_zscore(
    observed_nrl: float, pon, group_id: str = GENOME_WIDE_GROUP
) -> Optional[float]:
    """
    Compute NRL (Nucleosome Repeat Length) z-score.

    Compares sample NRL (from Rust FFT) against healthy reference
    from WPS background baseline.

    NRL ~180-200bp in healthy plasma; deviations indicate
    altered nucleosome spacing (cancer signature).

    Args:
        observed_nrl: NRL in bp from FFT analysis (nrl_bp column)
        pon: Loaded PonModel with wps_background_baseline
        group_id: Group identifier (default: the genome-wide group)

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.wps_background_baseline is None:
        return None

    return pon.wps_background_baseline.compute_nrl_zscore(observed_nrl, group_id)


def compute_periodicity_zscore(
    observed_periodicity: float, pon, group_id: str = GENOME_WIDE_GROUP
) -> Optional[float]:
    """
    Compute periodicity z-score.

    Periodicity measures the strength of nucleosome signal from FFT.
    Higher values indicate more regular nucleosome positioning.

    Args:
        observed_periodicity: Periodicity score from FFT
        pon: Loaded PonModel with wps_background_baseline
        group_id: Group identifier (default: the genome-wide group)

    Returns:
        Z-score or None if not computable
    """
    if pon is None or pon.wps_background_baseline is None:
        return None

    stats = pon.wps_background_baseline.get_periodicity_stats(group_id)
    if stats is None:
        return None

    mean, std = stats
    if std > 0:
        return (observed_periodicity - mean) / std
    return 0.0
