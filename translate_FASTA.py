#!/usr/bin/env python3
"""
Purpose: Translate a nucleotide FASTA file into amino acid sequences
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2023-05-23
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os
import sys
# from Bio import Align, AlignIO, SeqIO
# from collections import Counter, defaultdict
from evobioinfo import GAPS, hamming, IUPAC, IUPAC_AMBIG, NUCS_DEFINED, NUCS_INDETERMINATE, summary_string
from numpy import nan as NA
from typing import NamedTuple


# Usage
usage = """# -----------------------------------------------------------------------------
translate_FASTA.py - Translate a nucleotide FASTA file into amino acid sequences
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ translate_FASTA.py --help
    $ pydoc ./translate_FASTA.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ translate_FASTA.py --fasta_file=seqs.fasta -p
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    fasta_file: str
    out_file: str
    # p_dist: bool

# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Translate a nucleotide FASTA file into amino acid sequences. HELP: translate_FASTA.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--fasta_file',
                        metavar='FILE',
                        help='FASTA file containing sequence(s) [REQUIRED]',
                        required=True,
                        # nargs='+',
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Output file name [OPTIONAL]',
                        required=False,
                        type=str,
                        default='<INPUT>_aa.fasta')  # detect this later

    # parser.add_argument('-p',
    #                     '--p_dist',
    #                     help='Activate with alignments to calculate p-distance between sequences with identical IDs '
    #                          '[OPTIONAL]',
    #                     action='store_true')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments
    if not os.path.isfile(args.fasta_file):
        parser.error(f'\n### ERROR: fasta_file="{args.fasta_file}" does not exist')

    # FORM the output file name if not given
    if args.out_file == '<INPUT>_aa.fasta':
        args.out_file = os.path.splitext(args.fasta_file)[0] + '_aa.fasta'

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(fasta_file=args.fasta_file,
                out_file=args.out_file)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    fasta_file = args.fasta_file
    out_file = args.out_file
    # p_dist = args.p_dist

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print('LOG:')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')

    # arguments received
    print(f'LOG:fasta_file="{fasta_file}"')
    print(f'LOG:out_file="{out_file}"')
    # print(f'LOG:p_dist="{p_dist}"')

    # -------------------------------------------------------------------------
    # REGEX & TUPLES
    # re_date = re.compile(r'(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+).(\d+)')

    # fasta_file = sys.argv[1]
    translate_fasta(fasta_file, out_file)

    # -------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
def translate_fasta(input_file, output_file):
    translated_records = []

    # Parse the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        # Translate the sequence
        translated_seq = record.seq.translate()

        # Create a new record for the translated sequence
        new_record = SeqRecord(
            translated_seq,
            id=record.id,
            description=""
        )

        # Add the new record to our list
        translated_records.append(new_record)

    # Write the translated sequences to the output FASTA file
    SeqIO.write(translated_records, output_file, "fasta")


# -----------------------------------------------------------------------------
# CALL MAIN
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

