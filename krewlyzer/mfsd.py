import typer
from pathlib import Path
from typing import Optional
import logging
import pysam
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import track
import os

console = Console()
logging.basicConfig(level="INFO", handlers=[RichHandler(console=console)], format="%(message)s")
logger = logging.getLogger("mfsd")

def classify_read(read: pysam.AlignedSegment, pos: int, ref: str, alt: str) -> str:
    """
    Classify a read as Mutant or Wild-Type at a specific genomic position.
    Currently supports SNVs.
    pos: 0-based genomic position.
    """
    try:
        # Check if read covers the position
        if read.reference_start > pos or read.reference_end <= pos:
            return "Unknown"

        # Get aligned pairs to map reference position to query position
        # matches_only=True ensures we only get aligned bases (no deletions/insertions in cigar at this pos)
        # But for SNV, we want to see the base.
        # aligned_pairs returns (query_pos, ref_pos).
        
        # Optimization: use get_aligned_pairs(matches_only=True) might skip if it's a mismatch? 
        # No, matches_only=True means "not None", i.e., aligned columns. It includes mismatches.
        
        pairs = read.get_aligned_pairs(matches_only=True)
        for q_pos, r_pos in pairs:
            if r_pos == pos:
                base = read.query_sequence[q_pos].upper()
                if base == alt:
                    return "Mutant"
                if base == ref:
                    return "WildType"
                return "Other" # Different base
    except Exception:
        return "Unknown"
    return "Unknown"

def parse_input_file(input_file: Path, input_format: str) -> pd.DataFrame:
    """
    Parse VCF or MAF file into a standardized DataFrame.
    Returns DataFrame with columns: [chrom, pos, ref, alt] (pos is 0-based)
    """
    if input_format == "auto":
        if input_file.suffix.lower() in ['.vcf', '.gz']: # .vcf.gz
            input_format = "vcf"
        elif input_file.suffix.lower() in ['.maf', '.txt', '.tsv']:
            input_format = "maf"
        else:
            raise ValueError(f"Could not determine format for {input_file}. Please specify --format.")

    variants = []
    
    if input_format == "vcf":
        try:
            vcf = pysam.VariantFile(str(input_file))
            for record in vcf:
                # VCF is 1-based, pysam.VariantFile.pos is 1-based.
                # record.start is 0-based.
                # We handle only the first ALT allele for now if multiple exist.
                if len(record.alts) > 0:
                    variants.append({
                        'chrom': record.chrom,
                        'pos': record.start, # 0-based
                        'ref': record.ref,
                        'alt': record.alts[0]
                    })
        except Exception as e:
            logger.error(f"Error parsing VCF: {e}")
            raise typer.Exit(1)
            
    elif input_format == "maf":
        try:
            # MAF is tab-delimited.
            # Columns: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2
            df = pd.read_csv(input_file, sep='\t', comment='#')
            # Check required columns
            required = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
            if not all(col in df.columns for col in required):
                 # Try alternative column names if standard ones fail? 
                 # For now assume standard MAF.
                 raise ValueError(f"MAF file missing required columns: {required}")
            
            for _, row in df.iterrows():
                variants.append({
                    'chrom': str(row['Chromosome']),
                    'pos': int(row['Start_Position']) - 1, # MAF is 1-based
                    'ref': row['Reference_Allele'],
                    'alt': row['Tumor_Seq_Allele2']
                })
        except Exception as e:
            logger.error(f"Error parsing MAF: {e}")
            raise typer.Exit(1)
            
    return pd.DataFrame(variants)

def calc_mfsd(
    bam_file: Path,
    input_file: Path,
    output_file: Path,
    input_format: str = "auto",
    map_quality: int = 20
) -> None:
    """
    Calculate Mutant Fragment Size Distribution metrics.
    """
    try:
        logger.info(f"Parsing variants from {input_file}...")
        variants_df = parse_input_file(input_file, input_format)
        logger.info(f"Found {len(variants_df)} variants.")
        
        bam = pysam.AlignmentFile(str(bam_file), "rb")
        
        results = []
        
        for _, var in track(variants_df.iterrows(), total=len(variants_df), description="Processing variants..."):
            chrom = var['chrom']
            pos = var['pos']
            ref = var['ref']
            alt = var['alt']
            
            # Skip if ref/alt are not single bases (Indels) for now?
            # Let's try to handle SNVs primarily.
            if len(ref) > 1 or len(alt) > 1:
                # Simple Indel logic or skip?
                # classify_read logic above assumes SNV (base comparison).
                # Let's skip Indels for this version to ensure correctness of SNV first.
                # Or we can try.
                # For now, let's log warning and skip complex indels to avoid noise.
                # logger.warning(f"Skipping Indel at {chrom}:{pos} ({ref}->{alt}) - only SNVs supported currently.")
                continue

            mutant_lengths = []
            wt_lengths = []
            
            try:
                # Fetch reads
                # pos is 0-based.
                for read in bam.fetch(chrom, pos, pos + 1):
                    if read.mapping_quality < map_quality:
                        continue
                    if read.is_duplicate or read.is_unmapped or read.is_secondary:
                        continue
                        
                    cls = classify_read(read, pos, ref, alt)
                    
                    length = abs(read.template_length)
                    # template_length (TLEN) is the insert size.
                    # 0 means single ended or info not available.
                    if length == 0:
                        length = read.query_length # Fallback to read length?
                        
                    if cls == "Mutant":
                        mutant_lengths.append(length)
                    elif cls == "WildType":
                        wt_lengths.append(length)
                        
            except Exception as e:
                logger.warning(f"Error fetching reads at {chrom}:{pos}: {e}")
                continue
                
            # Calculate metrics
            n_mut = len(mutant_lengths)
            n_wt = len(wt_lengths)
            
            if n_mut > 0 and n_wt > 0:
                mut_mean = np.mean(mutant_lengths)
                wt_mean = np.mean(wt_lengths)
                delta_size = wt_mean - mut_mean
                
                # KS Test
                ks_stat, ks_pval = ks_2samp(mutant_lengths, wt_lengths)
            else:
                mut_mean = np.nan
                wt_mean = np.nan
                delta_size = np.nan
                ks_stat = np.nan
                ks_pval = np.nan
                
            results.append({
                'Chrom': chrom,
                'Pos': pos + 1, # Output 1-based for user convenience? Or 0-based?
                                # VCF/MAF users expect 1-based usually. Let's stick to 1-based for output.
                'Ref': ref,
                'Alt': alt,
                'Mut_Count': n_mut,
                'WT_Count': n_wt,
                'Mut_MeanSize': mut_mean,
                'WT_MeanSize': wt_mean,
                'Delta_Size': delta_size,
                'KS_Stat': ks_stat,
                'KS_Pval': ks_pval
            })
            
        # Write output
        out_df = pd.DataFrame(results)
        out_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"mFSD analysis complete. Results written to {output_file}")
        
    except Exception as e:
        logger.error(f"Fatal error in calc_mfsd: {e}")
        raise typer.Exit(1)

def mfsd(
    bam_path: Path = typer.Argument(..., help="Input BAM file"),
    input_file: Path = typer.Option(..., "--input", "-i", help="Input VCF or MAF file containing variants"),
    output: Path = typer.Option(..., "--output", "-o", help="Output file path (TSV)"),
    format: str = typer.Option("auto", "--format", "-f", help="Input format: 'auto', 'vcf', or 'maf'"),
    map_quality: int = typer.Option(20, "--map-quality", "-q", help="Minimum mapping quality")
) -> None:
    """
    Calculate Mutant Fragment Size Distribution (mFSD) features.
    Compares fragment sizes of mutant vs. wild-type reads at variant sites.
    """
    if not bam_path.exists():
        logger.error(f"BAM file not found: {bam_path}")
        raise typer.Exit(1)
    if not input_file.exists():
        logger.error(f"Input variant file not found: {input_file}")
        raise typer.Exit(1)
        
    # Create parent dir for output if needed
    output.parent.mkdir(parents=True, exist_ok=True)
    
    calc_mfsd(bam_path, input_file, output, format, map_quality)
