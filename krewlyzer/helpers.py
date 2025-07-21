import pysam
import itertools
import os
import numpy as np
import pandas as pd
import math
from collections import defaultdict
from rich.logging import RichHandler
import logging
from skmisc.loess import loess
import numpy as np

logging.basicConfig(level="INFO", handlers=[RichHandler()], format="%(message)s")
logger = logging.getLogger("krewlyzer-helpers")

def gc_correct(coverage, bias):
    """
    Perform GC bias correction on coverage values using LOESS regression.
    Logs errors and raises commonError if fitting fails.
    """
    covl = len(coverage)
    valid = [True for _ in range(covl)]
    temp_cov = []
    temp_bias = []
    for i in range(covl):
        if np.isnan(bias[i]):
            valid[i] = False
        else:
            temp_cov.append(coverage[i])
            temp_bias.append(bias[i])
    if not temp_cov or not temp_bias:
        logger.error("No valid coverage/bias values for GC correction.")
        raise commonError("No valid coverage/bias values for GC correction.")
    med = np.median(temp_cov)
    correct_cov = []
    try:
        i = np.arange(np.min(temp_bias), np.max(temp_bias), 0.001)
        coverage_trend = loess(temp_bias, temp_cov, span=0.75)
        coverage_trend.fit()
        coverage_model = loess(i, coverage_trend.predict(i, stderror=True).values)
        coverage_model.fit()
        coverage_pred = coverage_model.predict(temp_bias, stderror=True)
        pred = np.array(coverage_pred.values)
        coverage_corrected = temp_cov - pred + med
    except Exception as e:
        logger.error(f"GC correction failed: {e}")
        raise commonError(f"GC correction failed: {e}")
    i, j = 0, 0
    while i < covl:
        if valid[i]:
            if coverage_corrected[j] < 0:
                correct_cov.append(0)
            else:
                correct_cov.append(coverage_corrected[j])
            j += 1
        else:
            correct_cov.append(0)
        i += 1
    return correct_cov

class commonError(Exception):
    def __init__(self, message):
        logger.error(f"commonError: {message}")
        self.message = message

def maxCore(nCore=None):
    if nCore and nCore > 16:
        logger.warning("Requested nCore > 16; capping to 16.")
        return 16
    else:
        return nCore

# Alias for CLI import consistency
max_core = maxCore

def rmEndString(x, y):
    for item in y:
        if x.endswith(item):
            x = x.replace(item, "")
    return x

def isSoftClipped(cigar):
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

def GCcontent(seq):
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

def read_pair_generator(bam, region_string=None):
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

def reverse_seq(seq):
    r_seq = ''
    for i in seq:
        if i == 'A':
            r_seq += 'T'
        elif i == 'T':
            r_seq += 'A'
        elif i == 'C':
            r_seq += 'G'
        elif i == 'G':
            r_seq += 'C'
        else:
            r_seq += i
    return r_seq

def get_End_motif(Emotif, seq1, seq2):
    if seq1.count('N') + seq1.count('n') + seq2.count('N') + seq2.count('n') != 0:
        return Emotif
    seq2 = reverse_seq(seq2)
    if seq1 in Emotif.keys():
        Emotif[seq1] += 1
    if seq2 in Emotif.keys():
        Emotif[seq2] += 1
    return Emotif

def calc_MDS(inputEndMotifFile, outputfile):
    inputfile = pd.read_table(inputEndMotifFile, header=None, names=['bases', 'frequency'])
    k_mer = math.log(len(inputfile), 4)
    frequency = inputfile['frequency'].to_numpy()
    MDS = np.sum(-frequency * np.log2(frequency) / np.log2(4 ** k_mer))
    with open(outputfile, 'a') as f:
        f.write(inputEndMotifFile + '\t' + str(MDS) + '\n')

def get_Breakpoint_motif(Bpmotif, seq1, seq2):
    # seq1 and seq2 do not include N
    if seq1.count('N') + seq1.count('n') + seq2.count('N') + seq2.count('n') != 0:
        return Bpmotif
    seq2 = reverse_seq(seq2)
    if seq1 in Bpmotif.keys():
        Bpmotif[seq1] += 1
    if seq2 in Bpmotif.keys():
        Bpmotif[seq2] += 1
    return Bpmotif
