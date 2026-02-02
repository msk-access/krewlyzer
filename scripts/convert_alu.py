#!/usr/bin/env python3
"""
Filter RepeatMasker BED files to extract full-length Alu elements for WPS background.

Alu elements are used for "stacking" to generate global chromatin quality metrics.
We filter for full-length Alus (280-320bp) to ensure clean signal stacking.

Output includes subfamily classification for hierarchical stacking:
- Global_All: All Alus stacked together
- Family_AluY/AluS/AluJ: By evolutionary age (active vs silent chromatin)
- Chr{N}_All: Per-chromosome stacking for arm-level analysis

Input: UCSC RepeatMasker BED (chrom, start, end, repName, score, strand)
Output: BED7 with subfamily column: chrom, start, end, name, score, strand, subfamily
"""

import argparse
import gzip
import logging
from pathlib import Path
from collections import Counter

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Full-length Alu is ~300bp (280-320bp range)
MIN_LENGTH = 280
MAX_LENGTH = 320

# Valid chromosomes
VALID_CHROMS = set([str(i) for i in range(1, 23)] + ['X', 'Y'])


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome - strip chr prefix."""
    return chrom.lstrip('chr')


def classify_subfamily(rep_name: str) -> str:
    """
    Classify Alu element into major subfamily.
    
    AluY: Young, active chromatin (CpG-rich, often methylated)
    AluS: Middle age, most abundant
    AluJ: Old/ancient, often heterochromatin
    """
    # Strip Alu_ prefix if present
    name = rep_name.replace('Alu_', '').replace('Alu', '')
    
    # Classify by first letter after "Alu"
    if name.startswith('Y'):
        return 'AluY'
    elif name.startswith('S'):
        return 'AluS'
    elif name.startswith('J'):
        return 'AluJ'
    else:
        # FAM, FLAM, FRAM etc go to "Other"
        return 'AluOther'


def filter_alu(input_bed: Path, output_bed: Path) -> tuple[int, int, Counter]:
    """
    Filter RepeatMasker BED for full-length Alu elements.
    
    Returns: (total_alu_count, filtered_count, subfamily_counts)
    """
    total_alu = 0
    filtered = 0
    subfamily_counts = Counter()
    
    opener = gzip.open if str(input_bed).endswith('.gz') else open
    
    with opener(input_bed, 'rt') as f_in, open(output_bed, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            chrom, start, end, rep_name, score, strand = fields[:6]
            
            # Check if Alu family (starts with Alu)
            if not rep_name.startswith('Alu'):
                continue
            
            total_alu += 1
            
            # Check chromosome
            chrom_norm = normalize_chrom(chrom)
            if chrom_norm not in VALID_CHROMS:
                continue
            
            # Check length (full-length Alu only)
            length = int(end) - int(start)
            if length < MIN_LENGTH or length > MAX_LENGTH:
                continue
            
            filtered += 1
            
            # Classify subfamily
            subfamily = classify_subfamily(rep_name)
            subfamily_counts[subfamily] += 1
            
            # Output BED7: chrom, start, end, name, score, strand, subfamily
            f_out.write(f"{chrom_norm}\t{start}\t{end}\t{rep_name}\t{score}\t{strand}\t{subfamily}\n")
    
    return total_alu, filtered, subfamily_counts


def main():
    parser = argparse.ArgumentParser(description='Filter RepeatMasker for full-length Alu elements')
    parser.add_argument('--input', type=Path, required=True, help='Input RepeatMasker BED(.gz)')
    parser.add_argument('--output', type=Path, required=True, help='Output BED7 file with subfamily')
    
    args = parser.parse_args()
    
    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Filtering Alu elements from: {args.input}")
    total, filtered, subfamilies = filter_alu(args.input, args.output)
    
    logger.info(f"Total Alu elements: {total:,}")
    logger.info(f"Full-length Alu (280-320bp): {filtered:,} ({100*filtered/total:.1f}%)")
    logger.info("Subfamily breakdown:")
    for sf, count in sorted(subfamilies.items()):
        logger.info(f"  {sf}: {count:,} ({100*count/filtered:.1f}%)")
    logger.info(f"Output: {args.output}")


if __name__ == '__main__':
    main()
