"""
BAM utility functions shared across krewlyzer commands.

Provides filter compatibility checking and BAM statistics.
"""

from pathlib import Path
from typing import Any, Dict, List
import logging
import pysam

# Module-level logger so callers can set verbosity independently of this file.
logger = logging.getLogger("core.bam_utils")

# Fraction thresholds for issuing filter compatibility warnings.
_PROPER_PAIR_WARN_THRESHOLD = 0.1  # < 10% proper pairs → flag for duplex BAMs
_DUPLICATE_WARN_THRESHOLD = 0.5  # > 50% duplicates → flag for investigation
_MAPQ_FAIL_WARN_THRESHOLD = 0.5  # > 50% low MAPQ → suggest threshold change


def check_bam_compatibility(
    bam_path: Path,
    require_proper_pair: bool,
    skip_duplicates: bool,
    mapq_threshold: int,
    sample_size: int = 10000,
) -> Dict[str, Any]:
    """
    Sample reads from BAM to check filter compatibility.

    This function samples reads from a BAM file to determine if the
    current filter settings will work well. It's particularly useful
    for detecting duplex/consensus BAMs that have 0% proper pairs.

    Args:
        bam_path: Path to the input BAM file.
        require_proper_pair: Whether proper-pair filter is enabled.
        skip_duplicates: Whether duplicate filter is enabled.
        mapq_threshold: Minimum MAPQ threshold.
        sample_size: Number of reads to sample (default: 10000).

    Returns:
        Dict with keys:
        - total_sampled (int): number of reads sampled
        - would_pass (int): number that would pass all filters
        - pass_rate (float): fraction that would pass (0.0–1.0)
        - issues (List[str]): human-readable issue descriptions
        - suggested_flags (List[str]): CLI flag suggestions to fix issues
    """
    # Use typed local lists so mypy can verify .append() calls.
    # They are assembled into the return dict at the end of the function.
    issues: List[str] = []
    suggested_flags: List[str] = []

    logger.debug(
        f"Checking BAM compatibility: {bam_path.name} "
        f"(sample_size={sample_size}, mapq>={mapq_threshold}, "
        f"proper_pair={require_proper_pair}, skip_dups={skip_duplicates})"
    )

    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as e:
        msg = f"Could not open BAM: {e}"
        logger.warning(msg)
        issues.append(msg)
        return {
            "total_sampled": 0,
            "would_pass": 0,
            "pass_rate": 0.0,
            "issues": issues,
            "suggested_flags": suggested_flags,
        }

    # Sample from first chromosome that has reads
    sampled = 0
    passed = 0
    proper_pair_count = 0
    duplicate_count = 0
    low_mapq_count = 0

    for read in bam.fetch(until_eof=True):
        if sampled >= sample_size:
            break

        # Skip completely unusable reads (unmapped, secondary, supplementary)
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        sampled += 1

        # Check individual filters
        is_proper = read.is_proper_pair
        is_dup = read.is_duplicate
        has_mapq = read.mapping_quality >= mapq_threshold

        if is_proper:
            proper_pair_count += 1
        if is_dup:
            duplicate_count += 1
        if not has_mapq:
            low_mapq_count += 1

        # Would this read pass all active filters?
        passes = True
        if require_proper_pair and not is_proper:
            passes = False
        if skip_duplicates and is_dup:
            passes = False
        if not has_mapq:
            passes = False

        if passes:
            passed += 1

    bam.close()

    pass_rate = passed / sampled if sampled > 0 else 0.0
    logger.debug(
        f"  BAM sample: {sampled} reads, {passed} pass ({pass_rate:.1%}), "
        f"proper={proper_pair_count}, dups={duplicate_count}, lo_mapq={low_mapq_count}"
    )

    # Analyze sampled reads and generate actionable issue descriptions
    if sampled > 0:
        proper_rate = proper_pair_count / sampled
        dup_rate = duplicate_count / sampled
        mapq_fail_rate = low_mapq_count / sampled

        # Very few proper pairs → likely duplex/consensus BAM; flag for user
        if require_proper_pair and proper_rate < _PROPER_PAIR_WARN_THRESHOLD:
            msg = (
                f"Only {proper_rate:.1%} of reads are marked as proper pairs. "
                "This is common for duplex/consensus BAMs."
            )
            logger.warning(f"  BAM compatibility: {msg}")
            issues.append(msg)
            suggested_flags.append("--no-require-proper-pair")

        # High duplicate rate → investigate before running (possible QC issue)
        if skip_duplicates and dup_rate > _DUPLICATE_WARN_THRESHOLD:
            msg = (
                f"{dup_rate:.1%} of reads are marked as duplicates. "
                "Consider using --no-skip-duplicates if these are valid reads."
            )
            logger.warning(f"  BAM compatibility: {msg}")
            issues.append(msg)
            suggested_flags.append("--no-skip-duplicates")

        # High MAPQ failure rate → threshold may be too strict for this BAM
        if mapq_fail_rate > _MAPQ_FAIL_WARN_THRESHOLD:
            msg = (
                f"{mapq_fail_rate:.1%} of reads have MAPQ < {mapq_threshold}. "
                "Consider lowering --mapq threshold."
            )
            logger.warning(f"  BAM compatibility: {msg}")
            issues.append(msg)
            suggested_flags.append("--mapq 0")

    return {
        "total_sampled": sampled,
        "would_pass": passed,
        "pass_rate": pass_rate,
        "issues": issues,
        "suggested_flags": suggested_flags,
    }
