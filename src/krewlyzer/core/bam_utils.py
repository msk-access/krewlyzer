"""
BAM utility functions shared across krewlyzer commands.

Provides filter compatibility checking and BAM statistics.
"""

from pathlib import Path
import pysam


def check_bam_compatibility(
    bam_path: Path,
    require_proper_pair: bool,
    skip_duplicates: bool,
    mapq_threshold: int,
    sample_size: int = 10000
) -> dict:
    """
    Sample reads from BAM to check filter compatibility.
    
    This function samples reads from a BAM file to determine if the
    current filter settings will work well. It's particularly useful
    for detecting duplex/consensus BAMs that have 0% proper pairs.
    
    Args:
        bam_path: Path to the input BAM file
        require_proper_pair: Whether proper-pair filter is enabled
        skip_duplicates: Whether duplicate filter is enabled
        mapq_threshold: Minimum MAPQ threshold
        sample_size: Number of reads to sample (default: 10000)
    
    Returns:
        dict with:
        - total_sampled: number of reads sampled
        - would_pass: number that would pass filters
        - pass_rate: fraction that would pass
        - issues: list of issue descriptions
        - suggested_flags: list of suggested flag changes
    """
    result = {
        "total_sampled": 0,
        "would_pass": 0,
        "pass_rate": 0.0,
        "issues": [],
        "suggested_flags": []
    }
    
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as e:
        result["issues"].append(f"Could not open BAM: {e}")
        return result
    
    # Sample from first chromosome that has reads
    sampled = 0
    passed = 0
    proper_pair_count = 0
    duplicate_count = 0
    low_mapq_count = 0
    
    for read in bam.fetch(until_eof=True):
        if sampled >= sample_size:
            break
            
        # Skip completely unusable reads
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
        
        # Would this read pass all filters?
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
    
    result["total_sampled"] = sampled
    result["would_pass"] = passed
    result["pass_rate"] = passed / sampled if sampled > 0 else 0.0
    
    # Analyze issues
    if sampled > 0:
        proper_rate = proper_pair_count / sampled
        dup_rate = duplicate_count / sampled
        mapq_fail_rate = low_mapq_count / sampled
        
        # If very few proper pairs, suggest disabling that filter
        if require_proper_pair and proper_rate < 0.1:
            result["issues"].append(
                f"Only {proper_rate:.1%} of reads are marked as proper pairs. "
                "This is common for duplex/consensus BAMs."
            )
            result["suggested_flags"].append("--no-require-proper-pair")
        
        # If high duplicate rate, flag it
        if skip_duplicates and dup_rate > 0.5:
            result["issues"].append(
                f"{dup_rate:.1%} of reads are marked as duplicates. "
                "Consider using --no-skip-duplicates if these are valid reads."
            )
            result["suggested_flags"].append("--no-skip-duplicates")
        
        # If MAPQ fails too many reads
        if mapq_fail_rate > 0.5:
            result["issues"].append(
                f"{mapq_fail_rate:.1%} of reads have MAPQ < {mapq_threshold}. "
                "Consider lowering --mapq threshold."
            )
            result["suggested_flags"].append(f"--mapq 0")
    
    return result
