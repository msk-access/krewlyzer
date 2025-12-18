import pysam
import itertools
import os
import numpy as np
import pandas as pd
import math
from collections import defaultdict
from rich.logging import RichHandler
import logging
from pathlib import Path

logging.basicConfig(level="INFO", handlers=[RichHandler()], format="%(message)s")
logger = logging.getLogger("krewlyzer-helpers")

# Note: gc_correct function removed in v0.3.1 - GC correction now handled by Rust LOESS
# See count_fragments_gc_corrected() in _core module


class commonError(Exception):
    def __init__(self, message):
        logger.error(f"commonError: {message}")
        self.message = message

def maxCore(nCore: int | None = None) -> int | None:
    if nCore and nCore > 16:
        logger.warning("Requested nCore > 16; capping to 16.")
        return 16
    else:
        return nCore

# Alias for CLI import consistency
max_core = maxCore

def rmEndString(x: str, y: list[str]) -> str:
    for item in y:
        if x.endswith(item):
            x = x.replace(item, "")
    return x

def isSoftClipped(cigar: list[tuple[int, int]]) -> bool:
    """
    cigar information:
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    """
    for (op, count) in cigar:
        if op in [4, 5, 6]:
            return True
    return False

def GCcontent(seq: str) -> float:
    try:
        nA = seq.count("a") + seq.count("A")
        nT = seq.count("t") + seq.count("T")
        nG = seq.count("g") + seq.count("G")
        nC = seq.count("c") + seq.count("C")
        percent_GC = (nG + nC) / (nA + nT + nG + nC) if (nA + nT + nG + nC) > 0 else 0
        return percent_GC
    except Exception as e:
        logger.error(f"GCcontent calculation failed: {e}")
        return 0

def read_pair_generator(bam: pysam.AlignmentFile, region_string: str | None = None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    Reference: https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    try:
        for read in bam.fetch(region=region_string):
            if read.is_unmapped or read.is_qcfail or read.is_duplicate:
                continue
            if not read.is_paired or not read.is_proper_pair:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.mate_is_unmapped:
                continue
            if read.rnext != read.tid:
                continue
            if read.template_length == 0:
                continue
            if isSoftClipped(read.cigar):
                continue
            qname = read.query_name
            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                else:
                    yield read_dict[qname][0], read
                del read_dict[qname]
    except Exception as e:
        logger.error(f"Error during BAM read pair generation: {e}")
        return

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    trans_table = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return seq.translate(trans_table)[::-1]

def get_End_motif(Emotif: dict[str, int], end5: str, end3: str) -> dict[str, int]:
    """
    Update End Motif frequency dictionary.
    end5: 5' end of the fragment (from Read 1)
    end3: 3' end of the fragment (from Read 2, already reverse complemented to be on forward strand relative to fragment)
    """
    if "N" in end5 or "n" in end5 or "N" in end3 or "n" in end3:
        return Emotif
    
    # In cfDNAFE, they used:
    # seq2 = reverse_seq(seq2) -> which was just complement, not reverse.
    # And they passed forward_end3 twice.
    #
    # Correct logic:
    # We want the 4-mer at the 5' end and the 4-mer at the 3' end.
    # The 5' end sequence is just the sequence.
    # The 3' end sequence (on the + strand) is what we want.
    #
    # However, standard motif analysis often looks at the 5' end of the fragment on both strands.
    # If we treat the fragment as double stranded:
    # Strand 1: 5' [Seq] 3'
    # Strand 2: 3' [Seq_RC] 5'
    #
    # The 5' end of Strand 1 is `end5`.
    # The 5' end of Strand 2 corresponds to the 3' end of Strand 1, reverse complemented.
    #
    # If `end3` passed here is the sequence of the 3' end of the fragment on the forward strand:
    # Then the 5' end of the reverse strand is `reverse_complement(end3)`.
    #
    # Let's assume the caller passes the raw sequences from the reads.
    # Read 1 (Forward): 5' -> 3'. 5' end is the start of fragment.
    # Read 2 (Reverse): 5' -> 3'. 5' end is the other start of fragment (on reverse strand).
    #
    # So we should just take the 5' end of Read 1 and the 5' end of Read 2.
    #
    # But let's look at how `motif.py` calls this.
    # It will be updated to pass:
    # forward_end5 (Read 1 5' end)
    # reverse_end5 (Read 2 5' end)
    #
    # So we just count both of them.
    
    if end5 in Emotif:
        Emotif[end5] += 1
    if end3 in Emotif:
        Emotif[end3] += 1
    return Emotif

def calc_MDS(inputEndMotifFile: str | Path, outputfile: str | Path) -> None:
    inputfile = pd.read_table(inputEndMotifFile, header=None, names=['bases', 'frequency'])
    k_mer = math.log(len(inputfile), 4)
    frequency = inputfile['frequency'].to_numpy()
    MDS = np.sum(-frequency * np.log2(frequency) / np.log2(4 ** k_mer))
    with open(outputfile, 'a') as f:
        f.write(str(inputEndMotifFile) + '\t' + str(MDS) + '\n')

def get_Breakpoint_motif(Bpmotif: dict[str, int], seq1: str, seq2: str) -> dict[str, int]:
    """
    Update Breakpoint Motif frequency dictionary.
    seq1: Sequence around the 5' end of the fragment.
    seq2: Sequence around the 3' end of the fragment.
    """
    if "N" in seq1 or "n" in seq1 or "N" in seq2 or "n" in seq2:
        return Bpmotif
    
    # Similar to End Motif, we just count the motifs at both breakpoints.
    # The caller should ensure seq1 and seq2 are the correct sequences surrounding the breakpoints.
    
    if seq1 in Bpmotif:
        Bpmotif[seq1] += 1
    if seq2 in Bpmotif:
        Bpmotif[seq2] += 1
    return Bpmotif
