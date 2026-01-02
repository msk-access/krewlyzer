#!/usr/bin/env python3
"""
Convert Ensembl BioMart TSS/CTCF files to ML-optimized WPS anchor BED6 format.

Key Optimizations:
1. Canonical TSS: One per gene (longest transcript)
2. BED6 format: chrom, start, end, name, score, strand
3. TSS naming: TSS|GeneName|TranscriptID
4. CTCF: Midpoint centered, strand="."
5. Strand-aware for ML vector flipping

Output compatible with strand-aware WPS processing.
"""

import argparse
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome - strip chr prefix if present."""
    return str(chrom).lstrip('chr')


def convert_tss_canonical(input_tsv: Path) -> pd.DataFrame:
    """
    Convert TSS TSV to BED6 format with canonical selection.
    
    Strategy:
    1. Group by Gene stable ID
    2. Select longest transcript per gene
    3. Output 1bp TSS point (for windowing by Rust)
    4. Include strand for ML vector flipping
    
    Output: BED6 (chrom, start, end, name, score, strand)
    """
    logger.info(f"Reading TSS file: {input_tsv}")
    df = pd.read_csv(input_tsv, sep='\t', compression='gzip', low_memory=False)
    
    logger.info(f"Total transcripts: {len(df)}")
    
    # Required columns
    gene_id_col = 'Gene stable ID'
    transcript_id_col = 'Transcript stable ID'
    gene_name_col = 'Gene name'
    chrom_col = 'Chromosome/scaffold name'
    tss_col = 'Transcription start site (TSS)'
    strand_col = 'Strand'
    length_col = 'Transcript length (including UTRs and CDS)'
    
    # Sort by gene and transcript length (descending) to get longest first
    df = df.sort_values([gene_id_col, length_col], ascending=[True, False])
    
    # Keep only the canonical (longest) transcript per gene
    df_canonical = df.drop_duplicates(subset=gene_id_col, keep='first').copy()
    logger.info(f"Canonical transcripts (one per gene): {len(df_canonical)}")
    
    # Map strand: Ensembl uses 1/-1, BED uses +/-
    strand_map = {1: '+', -1: '-', '1': '+', '-1': '-'}
    df_canonical['bed_strand'] = df_canonical[strand_col].map(strand_map)
    
    # Filter valid chromosomes (1-22, X, Y)
    df_canonical['chrom'] = df_canonical[chrom_col].astype(str).apply(normalize_chrom)
    valid = df_canonical['chrom'].str.match(r'^([1-9]|1[0-9]|2[0-2]|X|Y)$')
    df_canonical = df_canonical[valid].copy()
    
    # BED format is 0-based, half-open
    # TSS point: start = TSS-1, end = TSS (1bp window)
    df_canonical['bed_start'] = df_canonical[tss_col].astype(int) - 1
    df_canonical['bed_end'] = df_canonical[tss_col].astype(int)
    
    # Name: TSS|GeneName|TranscriptID
    df_canonical['name'] = 'TSS|' + df_canonical[gene_name_col].fillna(df_canonical[gene_id_col]).astype(str) + '|' + df_canonical[transcript_id_col].astype(str)
    
    # Score: use "." for now
    df_canonical['score'] = '.'
    
    # Build BED6
    bed = df_canonical[['chrom', 'bed_start', 'bed_end', 'name', 'score', 'bed_strand']].copy()
    bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    
    # Ensure start >= 0
    bed['start'] = bed['start'].clip(lower=0)
    
    logger.info(f"TSS: {len(bed)} canonical anchors (one per gene)")
    return bed


def convert_ctcf_midpoint(input_tsv: Path) -> pd.DataFrame:
    """
    Extract CTCF sites from regulatory features, centered on midpoint.
    
    Strategy:
    1. Filter for CTCF binding sites
    2. Calculate midpoint
    3. Output 1bp window at midpoint
    4. Strand = "." (unstranded for symmetric analysis)
    
    Output: BED6 (chrom, start, end, name, score, strand)
    """
    logger.info(f"Reading CTCF/regulatory file: {input_tsv}")
    df = pd.read_csv(input_tsv, sep='\t', compression='gzip', low_memory=False)
    
    # Column names
    chrom_col = 'Chromosome/scaffold name'
    start_col = 'Start (bp)'
    end_col = 'End (bp)'
    type_col = 'Feature type'
    
    # Filter for CTCF only
    ctcf_mask = df[type_col].str.contains('CTCF', case=False, na=False)
    ctcf = df[ctcf_mask].copy()
    logger.info(f"CTCF entries: {len(ctcf)} of {len(df)}")
    
    # Filter valid chromosomes
    ctcf['chrom'] = ctcf[chrom_col].astype(str).apply(normalize_chrom)
    valid = ctcf['chrom'].str.match(r'^([1-9]|1[0-9]|2[0-2]|X|Y)$')
    ctcf = ctcf[valid].copy()
    
    # Calculate midpoint (centered 1bp window)
    ctcf['midpoint'] = ((ctcf[start_col] + ctcf[end_col]) // 2).astype(int)
    ctcf['start'] = ctcf['midpoint'] - 1
    ctcf['end'] = ctcf['midpoint']
    
    # Name: CTCF|chrom:position
    ctcf['name'] = 'CTCF|' + ctcf['chrom'] + ':' + ctcf['midpoint'].astype(str)
    
    # Score and strand (unstranded)
    ctcf['score'] = '.'
    ctcf['strand'] = '.'
    
    # Build BED6
    bed = ctcf[['chrom', 'start', 'end', 'name', 'score', 'strand']].copy()
    
    # Ensure start >= 0
    bed['start'] = bed['start'].clip(lower=0)
    
    logger.info(f"CTCF: {len(bed)} midpoint anchors")
    return bed


def merge_and_save(tss_bed: pd.DataFrame, ctcf_bed: pd.DataFrame, output_bed: Path):
    """Merge TSS and CTCF into single sorted BED6 file."""
    merged = pd.concat([tss_bed, ctcf_bed], ignore_index=True)
    
    # Sort by chromosome then position
    def chrom_sort_key(c):
        if c == 'X':
            return 23
        elif c == 'Y':
            return 24
        else:
            try:
                return int(c)
            except:
                return 99
    
    merged['_sort'] = merged['chrom'].apply(chrom_sort_key)
    merged = merged.sort_values(['_sort', 'start']).drop(columns=['_sort'])
    
    # Write BED6 (no header)
    merged.to_csv(output_bed, sep='\t', index=False, header=False)
    
    tss_count = len(tss_bed)
    ctcf_count = len(ctcf_bed)
    logger.info(f"Merged anchors: {len(merged)} total ({tss_count} TSS + {ctcf_count} CTCF)")
    logger.info(f"Output: {output_bed}")


def main():
    parser = argparse.ArgumentParser(description='Convert Ensembl BioMart files to ML-optimized WPS anchors (BED6)')
    parser.add_argument('--tss', type=Path, required=True, help='TSS TSV file from BioMart')
    parser.add_argument('--ctcf', type=Path, required=True, help='Regulatory features TSV from BioMart')
    parser.add_argument('--output', type=Path, required=True, help='Output BED6 file')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.tss.exists():
        raise FileNotFoundError(f"TSS file not found: {args.tss}")
    if not args.ctcf.exists():
        raise FileNotFoundError(f"CTCF file not found: {args.ctcf}")
    
    # Convert with canonical selection
    tss_bed = convert_tss_canonical(args.tss)
    ctcf_bed = convert_ctcf_midpoint(args.ctcf)
    
    # Merge and save
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merge_and_save(tss_bed, ctcf_bed, args.output)


if __name__ == '__main__':
    main()
