#!/usr/bin/env python3
"""
Purpose: Evolutionary bioinformatics module file
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-02-05

Details: Library for evolutionary bioinformatics

Available tuples:
- NUCS_DEFINED
- GAPS
- IUPAC
- IUPAC_AMBIG
- STOP_CODONS
- AA_STOP
- AA_AMBIG

Available functions:
- hamming(): returns the number of differences between two aligned sequences
- summary_string(): returns summary statistics on an input list of numbers
"""

import argparse
import dendropy
import gffutils
import gzip
import math
import numpy as np
import operator as op
import os
import re
import sys
import vcf
from Bio import Align, AlignIO, Seq, SeqIO
from collections import Counter, defaultdict
from ete3 import Tree
from ete3.parser import newick
from itertools import starmap, zip_longest
from numpy import nan as NA
from pprint import pprint
from scipy.stats import binom
from typing import Any, Dict, List, NamedTuple, Set, TextIO


# -----------------------------------------------------------------------------
# TUPLES
# -----------------------------------------------------------------------------

NUCS_DEFINED = ('A', 'C', 'G', 'T', 'U')
NUCS_INDETERMINATE = ('N')
GAPS = ('-')
IUPAC = ('A', 'C', 'G', 'T', 'U', 'N', 'B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'Y')
IUPAC_AMBIG = ('N', 'B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'Y')
# IUPAC_AMBIG_NONN = ('B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'Y')
# CODONS_START = ('ATG', 'AUG')
# CODONS_START_ALT = ('CTG', 'CUG')  # TODO
STOP_CODONS = ('TAA', 'UAA', 'TAG', 'UAG', 'TGA', 'UGA')
AA_STOP = ('*')
AA_AMBIG = ('X')


"""
IUPAC KEY:
'A'=A
'B'=C,G,T; all bases except A (B follows A)
'C'=C
'D'=A,G,T; all bases except C (D follows C)
'G'=G
'H'=A,C,T; all bases except G (H follows G)
'K'=G,T (Keto)
'M'=A,C (aMino)
'N'=A,C,G,T (aNy)
'R'=A,G (puRine)
'S'=C,G (Strongly bonded)
'T'=T
'U'=U
'V'=A,C,G; all bases except T or U (V follows U)
'W'=A,T (Weakly bonded)
'Y'=C,T (pYrimidine)
"""


# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def hamming(seq1: str, seq2: str) -> int:
    """ Calculate Hamming distance - Ken Youens-Clark """
    """
Purpose: <FILL IN>

Details: <FILL IN>

Args:
    <ARG1: DESCRIPTION>
    <ARG2: DESCRIPTION>
    """
    return sum(starmap(op.ne, zip_longest(seq1, seq2)))


# -----------------------------------------------------------------------------
def test_hamming() -> None:
    """ Test hamming() """
    assert hamming('', '') == 0
    assert hamming('AC', 'ACGT') == 2
    assert hamming('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT') == 7


# --------------------------------------------------
def find_kmers(seq: str, k: int) -> List[str]:
    """ Find k-mers in string - Ken Youens-Clark """
    n = len(seq) - k + 1
    return [] if n < 1 else [seq[i:i + k] for i in range(n)]


# --------------------------------------------------
def test_find_kmers() -> None:
    """ Test find_kmers() """

    assert find_kmers('', 1) == []
    assert find_kmers('ACTG', 1) == ['A', 'C', 'T', 'G']
    assert find_kmers('ACTG', 2) == ['AC', 'CT', 'TG']
    assert find_kmers('ACTG', 3) == ['ACT', 'CTG']
    assert find_kmers('ACTG', 4) == ['ACTG']
    assert find_kmers('ACTG', 5) == []


# -----------------------------------------------------------------------------
def summary_string(num_list: list, decimals: int = 3, suffix: str = "") -> str:
    """ Calculate summary statistics on an input list of numbers """
    n_num = len(num_list)

    # All numbers?
    all_numeric = all(str(x).replace('-', '').replace('.', '').isnumeric() for x in num_list)

    if n_num == 0 or not all_numeric:
        min_num = 'NA'
        Q1_num = 'NA'
        mean_num = 'NA'
        std_num = 'NA'
        median_num = 'NA'
        Q3_num = 'NA'
        max_num = 'NA'
    else:  # all work if only 1 element
        min_num = np.round(np.min(num_list), decimals)
        Q1_num = np.round(np.quantile(num_list, 0.25), decimals)
        mean_num = np.round(np.mean(num_list), decimals)
        std_num = np.round(np.std(num_list), decimals)
        median_num = np.round(np.median(num_list), decimals)
        Q3_num = np.round(np.quantile(num_list, 0.75), decimals)
        max_num = np.round(np.max(num_list), decimals)

    return f'n={n_num};min={min_num}{suffix};Q1={Q1_num}{suffix};mean={mean_num}{suffix};std={std_num}{suffix}' \
           f';median={median_num}{suffix};Q3={Q3_num}{suffix};max={max_num}{suffix}'


# -----------------------------------------------------------------------------
def test_summary_string() -> None:
    """ Test summary_string() """
    assert summary_string([]) == 'n=0;min=NA;Q1=NA;mean=NA;std=NA;median=NA;Q3=NA;max=NA'
    assert summary_string(['']) == 'n=1;min=NA;Q1=NA;mean=NA;std=NA;median=NA;Q3=NA;max=NA'
    assert summary_string(['Hello']) == 'n=1;min=NA;Q1=NA;mean=NA;std=NA;median=NA;Q3=NA;max=NA'
    assert summary_string([1]) == 'n=1;min=1;Q1=1.0;mean=1.0;std=0.0;median=1.0;Q3=1.0;max=1'
    assert summary_string([1, 5.6, 298.3, 55, 9]) == 'n=5;min=1.0;Q1=5.6;mean=73.78;std=113.933;median=9.0;Q3=55.0;max=298.3'