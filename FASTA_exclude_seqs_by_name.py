#!/usr/bin/env python3
"""
Purpose: Exclude a given set of sequences from a FASTA file
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-01-20

Details: Exclude a given set of sequences from a FASTA file

Returns:
    - STDOUT: documentation of retained and excluded sequence names
    - FILE: --out_file, a FASTA file with excluded sequences removed
"""

import argparse
import os
import re
import sys
from Bio import SeqIO
from typing import List, NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
FASTA_exclude_seqs_by_name.py - Exclude a given set of sequences from a FASTA file
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_exclude_seqs_by_name.py --help
    $ pydoc ./FASTA_exclude_seqs_by_name.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_exclude_seqs_by_name.py --seq_file=seqs.fasta --seq_names=seq_names.txt --out_file=seqs_exclude.fasta
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    seq_file: TextIO
    seq_names: str
    out_file: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Exclude a given set of sequences from a FASTA file. HELP: FASTA_exclude_seqs_by_name.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--seq_file',
                        metavar='FILE',
                        help='FASTA file(s) [REQUIRED]',
                        required=True,
                        type=argparse.FileType('rt'))

    parser.add_argument('-I',
                        '--seq_names',
                        metavar='str/FILE',
                        help='A comma-separated list of sequence IDs (names; beginning of FASTA headers) [REQUIRED]',
                        required=True,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Output file prefix for printing sequences; otherwise STDOUT [OPTIONAL]',
                        required=False,
                        type=str,
                        default='')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # if seq_names is file, convert to string
    if os.path.isfile(args.seq_names):
        args.seq_names = open(args.seq_names).read().rstrip()

    # die if seq_names contains anything besides word, dash, whitespace, and commas
    re_NOT_seqID = re.compile(r'[^\w\-\s\n,]')
    if re_NOT_seqID.search(args.seq_names):
        parser.error(f'\n### ERROR: custom_sites="{args.seq_names[:20]}..." contains invalid characters')

    # die if out_file exists
    if args.out_file != '' and os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(seq_file=args.seq_file,
                seq_names=args.seq_names,
                out_file=args.out_file)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    seq_fh = args.seq_file
    seq_names = args.seq_names
    out_file = args.out_file

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:seq_file="{seq_fh.name}"')
    print(f'LOG:seq_names="{seq_names}"')
    print(f'LOG:out_file="{out_file}"')

    # -------------------------------------------------------------------------
    # INITIALIZE set of seq_names
    seq_name_list = list(seq_names.split(','))
    print(f'LOG:seq_names_count={len(seq_name_list)}')
    seq_name_set = set(seq_name_list)
    unique_seq_names_count = len(seq_name_set)
    print(f'LOG:unique_seq_names_count={unique_seq_names_count}')

    # -------------------------------------------------------------------------
    # PRINT the wanted recs

    # OPEN out_file for writing, to which SeqIO.write() will append
    out_fh = sys.stdout
    if out_file != '':
        out_fh = open(out_file, 'wt')

    # initialize counters
    nseqs = 0
    nseqs_included = 0
    seq_names_included: List[str] = []
    nseqs_excluded = 0
    seq_names_excluded: List[str] = []
    for rec in SeqIO.parse(seq_fh, 'fasta'):
        nseqs += 1
        if rec.id in seq_name_set:  # EXCLUDE
            nseqs_excluded += 1
            seq_names_excluded.append(rec.id)
        else:  # INCLUDE
            nseqs_included += 1
            seq_names_included.append(rec.id)
            SeqIO.write(rec, out_fh, 'fasta')  # APPENDS to an open filehandle

    print(f'LOG:nseqs={nseqs}')
    print(f'LOG:nseqs_included={nseqs_included}')
    print(f'LOG:seq_names_included={",".join(seq_names_included)}')
    print(f'LOG:nseqs_excluded={nseqs_excluded}')
    print(f'LOG:seq_names_excluded={",".join(seq_names_excluded)}')

    # -------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
# CALL MAIN
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
