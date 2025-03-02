#!/usr/bin/env python3
"""
Purpose: Mask specified range(s) of sites in a FASTA file
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2023-03-09
"""

import argparse
import os
import sys
import time
from Bio import AlignIO
from typing import NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
FASTA_mask_sites.py - Mask specified range(s) of sites in a FASTA file
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_mask_sites.py --help
    $ pydoc ./FASTA_mask_sites.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_mask_sites.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    in_file: TextIO
    out_file: str
    start: str  # string because may be multiple, comma-separated
    end: str
    mask_char: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Mask specified range(s) of sites in a FASTA file. HELP: FASTA_mask_sites.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--in_file',
                        metavar='FILE',
                        help='FASTA multiple sequence alignment (MSA) [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of FASTA output file to contain masked sequences [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-s',
                        '--start',
                        metavar='int(s)',
                        help='Start site(s) for the range(s) to be masked; 1-based; comma-separated if multiple '
                             '(e.g., for sites 9-11 and 27-30, use "9,27") [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-e',
                        '--end',
                        metavar='int(s)',
                        help='End site(s) for the range(s) to be masked; 1-based; comma-separated if multiple '
                             '(e.g., for sites 9-11 and 27-30, use (e.g., "11,30") [REQUIRED]',
                        required=True,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-k',
                        '--mask_char',
                        metavar='str',
                        help='Character used for masking [OPTIONAL]',
                        required=False,
                        type=str,
                        default='N')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # CONVERT start(s) and end(s) to lists
    start_list = [int(start) for start in args.start.split(',')]
    end_list = [int(end) for end in args.end.split(',')]

    # All starts <= ends
    for start, end in zip(start_list, end_list):
        if start < 1:
            parser.error(f'\n### ERROR: start={start} smaller than position 1')
        if end < start:
            parser.error(f'\n### ERROR: start={start} is greater than end={end}')

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(in_file=args.in_file[0],
                out_file=args.out_file,
                start=args.start,
                end=args.end,
                mask_char=args.mask_char)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    in_fh = args.in_file
    out_file = args.out_file
    start = args.start
    end = args.end
    mask_char = args.mask_char

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # CONVERT start(s) and end(s) to lists
    start_list = [int(start) for start in start.split(',')]
    end_list = [int(end) for end in end.split(',')]

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    #    # group ID(s)
    #    regex_group_list = []
    #    for this_group_key in group_key_list:  # group_key.split(','):
    #        regex_group = re.compile(f'{this_group_key}="([^"]+)"')
    #        regex_group_list.append(regex_group)
    #
    #    # Nucleotides
    #    defined_nucs = ('A', 'C', 'G', 'T', 'U')

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:in_file="{in_fh.name}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:start={start}')
    print(f'LOG:end={end}')
    print(f'LOG:mask_char="{mask_char}"')

    # -------------------------------------------------------------------------
    # TODO: CHECK FOR DUPLICATE SEQUENCE IDs

    # -------------------------------------------------------------------------
    # OPEN FASTA FILE; LOOP SEQS; for each, FIND ITS BED AND MASK GIVEN COVERAGE

    # get recs - using AlignIO for now because it IS an MSA
    recs = AlignIO.read(in_fh, 'fasta')  # unit is ALIGNMENT, and .read() because we have only one ALN per file

    # Get alignment (sequence) length and make sure no ranges exceed
    seq_length = recs.get_alignment_length()  # None

    for iv_end in end_list:
        if iv_end > seq_length:
            sys.exit(f'\n### ERROR: end={iv_end} falls outside the sequence, length={seq_length}')

    # track num sequences examined
    nseqs = 0

    # FASTA: loop recs and MASK
    for rec in recs:
        nseqs += 1

        for iv_start, iv_end in zip(start_list, end_list):
            # convert to indices
            start_idx = iv_start - 1
            end_idx = iv_end

            # compute segment length
            seg_length = end_idx - start_idx

            # mask and overwrite
            rec.seq = rec.seq[:start_idx] + (seg_length * mask_char) + rec.seq[end_idx:]

        # PRINT MODIFIED rec to new file
        # SeqIO.write(rec, out_file, 'fasta')

    # # PRINT MODIFIED recs to new file
    AlignIO.write(recs, out_file, 'fasta')

    # -------------------------------------------------------------------------
    # LOG sequence counts and groups
    # Number of sequences
    print(f'LOG:nseqs={nseqs}')

    # -------------------------------------------------------------------------
    # TIME & DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_time_per_seq = elapsed_time / nseqs

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time / 60, ndigits=1)} minutes; {round(elapsed_time / 60 / 60, ndigits=1)} hours)')

    print(f'TIME ELAPSED PER SEQ: {round(elapsed_time_per_seq, ndigits=9)} seconds')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
