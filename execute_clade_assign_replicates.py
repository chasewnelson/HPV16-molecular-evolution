#!/usr/bin/env python3
"""
Purpose: Execute replicates using clade_assigner.py in working directory
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-12-14
"""

import argparse
# import dendropy
# import gffutils
# import gzip
# import math
# import numpy as np
# import operator as op
import os
# import pandas as pd
# import re
import sys
import time
# import vcf
# from Bio import Align, AlignIO, Seq, SeqIO
from collections import Counter, defaultdict
# from ete3 import Tree
# from ete3.parser import newick
from itertools import starmap, zip_longest
# from numpy import nan as NA
# from pprint import pprint
# from scipy.stats import binom
from subprocess import getstatusoutput
from typing import Any, Dict, List, NamedTuple, Set, TextIO


usage = """# -----------------------------------------------------------------------------
execute_clade_assign_replicates.py - Execute replicates using clade_assigner.py in working directory
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ execute_clade_assign_replicates.py --help
    $ pydoc ./execute_clade_assign_replicates.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ execute_clade_assign_replicates.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    replicates: int
    permutations: int
    tree: TextIO
    representatives: str
    start: int


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Execute replicates using clade_assigner.py. HELP: execute_clade_assign_replicates.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-r',
                        '--replicates',
                        metavar='int',
                        help='Number of replicates to perform, each using the given number of permutations [REQUIRED]',
                        required=True,
                        type=int)

    parser.add_argument('-p',
                        '--permutations',
                        metavar='int',
                        help='Number of permutations (tree root / pruning order) for clade_assigner.py [REQUIRED]',
                        required=True,
                        type=int)

    parser.add_argument('-t',
                        '--tree',
                        metavar='FILE',
                        help='NEWICK file containing the tree(s) [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-R',
                        '--representatives',
                        metavar='FILE|str',
                        help='Sequence ID(s) (comma-separated) of clade representatives [REQUIRED]',
                        required=True,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-s',
                        '--start',
                        metavar='int',
                        help='Starting replicate number (DEFAULT: 1) [OPTIONAL]',
                        required=False,
                        type=int,
                        default=1)

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    if not args.replicates > 0:
        parser.error(f'\n### ERROR: replicates="{args.replicates}" must be > 0')

    if not args.permutations > 0:
        parser.error(f'\n### ERROR: permutations="{args.permutations}" must be > 0')

    if not os.path.isfile(args.tree):
        parser.error(f'\n### ERROR: file tree="{args.tree}" does not exist')

    if not args.start >= 1:
        parser.error(f'\n### ERROR: start="{args.start}" must be >= 1')

    return Args(replicates=args.replicates,
                permutations=args.permutations,
                tree=args.tree,
                representatives=args.representatives,
                start=args.start)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    replicates = args.replicates
    permutations = args.permutations
    tree = args.tree
    representatives = args.representatives
    start = args.start

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # PASS

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:replicates="{replicates}"')
    print(f'LOG:permutations="{permutations}"')
    print(f'LOG:tree="{tree}"')
    print(f'LOG:representatives="{representatives}"')
    print(f'LOG:start="{start}"')

    # -------------------------------------------------------------------------
    # PERFORM REPLICATES IN WORKING DIRECTORY
    for replicate in range(start, start + replicates):
        COMMAND = f'clade_assigner.py -t {tree} -r {representatives} -p {permutations} ' \
                  f'-o p{permutations}_r{replicate} -S > p{permutations}_r{replicate}.out'
        rv, out = getstatusoutput(COMMAND)
        assert rv == 0

    # -------------------------------------------------------------------------
    # DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time / 60, ndigits=1)} minutes; {round(elapsed_time / 60 / 60, ndigits=1)} hours)')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')

# -----------------------------------------------------------------------------
# CALL MAIN
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
