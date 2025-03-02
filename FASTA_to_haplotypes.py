#!/usr/bin/env python3
"""
Purpose: Tally the unique haplotypes in a FASTA file (nucleotide or amino acid)
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2023-05-23
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import argparse
import os
import sys
# from Bio import Align, AlignIO, SeqIO
# from collections import Counter, defaultdict
# from evobioinfo import GAPS, hamming, IUPAC, IUPAC_AMBIG, NUCS_DEFINED, NUCS_INDETERMINATE, summary_string
# from numpy import nan as NA
from typing import NamedTuple  # FASTA_to_haplotypes.py -i HPV16_A_n4013_E6_aa.fasta


# Usage
usage = """# -----------------------------------------------------------------------------
FASTA_to_haplotypes.py - Tally the unique haplotypes in a FASTA file (nucleotide or amino acid)
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_to_haplotypes.py --help
    $ pydoc ./FASTA_to_haplotypes.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_to_haplotypes.py -i HPV16_A_n4013_E6_aa.fasta
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    fasta_file: str
    out_file: str
    prefix: str
    header: str
    # p_dist: bool


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Tally the unique haplotypes in a FASTA file (nucleotide or amino acid). HELP: FASTA_to_haplotypes.py --help',
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
                        default='<INPUT>_haplotypes.tsv')  # detect this later

    parser.add_argument('-p',
                        '--prefix',
                        metavar='str',
                        help='Prefix of FASTA name for most common haplotype [OPTIONAL]',
                        required=False,
                        type=str,
                        default='<INPUT>_haplotype1')  # detect this later

    parser.add_argument('-H',
                        '--header',
                        metavar='str',
                        help='Header for FASTA of most common haplotype [OPTIONAL]',
                        required=False,
                        type=str,
                        default='<INPUT>_haplotype1')  # detect this later

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
    if args.out_file == '<INPUT>_haplotypes.tsv':
        args.out_file = os.path.splitext(args.fasta_file)[0] + '_haplotypes.tsv'

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    # FORM the prefix FASTA name if not given
    if args.prefix == '<INPUT>_haplotype1':
        args.prefix = os.path.splitext(args.fasta_file)[0] + '_haplotype1'

    # FORM the HEADER of FASTA
    if args.header == '<INPUT>_haplotype1':
        args.header = os.path.splitext(args.fasta_file)[0] + '_haplotype1'

    return Args(fasta_file=args.fasta_file,
                out_file=args.out_file,
                prefix=args.prefix,
                header=args.header)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    fasta_file = args.fasta_file
    out_file = args.out_file
    prefix = args.prefix
    header = args.header
    # p_dist = args.p_dist

    # form output fasta - HACK
    out_fasta = out_file.replace('.tsv', '.fasta')

    # ensure out_fasta not present
    if os.path.isfile(out_fasta):
        sys.exit(f'\n### ERROR: out_file="{out_fasta}" already exists')

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
    print(f'LOG:out_fasta="{out_fasta}"')
    print(f'LOG:prefix="{prefix}"')
    print(f'LOG:header="{header}"')
    # print(f'LOG:p_dist="{p_dist}"')

    # -------------------------------------------------------------------------
    # REGEX & TUPLES
    # re_date = re.compile(r'(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+).(\d+)')

    # -------------------------------------------------------------------------
    # BUSINESS
    sequences = read_fasta_file(fasta_file)
    sequence_counts = count_unique_sequences(sequences)

    # write the tsv
    write_sequence_counts(sequence_counts, out_file)

    # write the fasta
    write_seq_counts_fasta(sequence_counts, out_fasta)

    # get and print min_X_count
    min_X_count = get_min_x_count(sequence_counts)
    print(f"Minimum X_count: {min_X_count}")

    # get most common seq and its count
    most_common_tuple = sequence_counts.most_common(1)[0]
    most_common_sequence, most_common_count = most_common_tuple[0], most_common_tuple[1]

    # write the most common sequence (haplotype)
    if header == prefix:
        header = f'{prefix}_n{most_common_count}'

    out_fasta_name = f'{prefix}_n{most_common_count}.fasta'
    most_common_Seq = Seq(most_common_sequence)
    record = SeqRecord(most_common_Seq, id=header, description="")
    SeqIO.write(record, out_fasta_name, "fasta")

    # WARNING if the most commmon sequence didn't have the minimum nujber of Xs
    if count_x_chars(most_common_Seq) > min_X_count:
        print(f'\n### WARNING: min_X_count={min_X_count} but most common seq has X={count_x_chars(most_common_Seq)}')
    else:
        print(f'CONFIRMED: most common seq has min_X_count={min_X_count}')

    # -------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
def read_fasta_file(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences


# -----------------------------------------------------------------------------
def count_unique_sequences(sequences):
    sequence_counts = Counter(sequences)
    return sequence_counts


# -----------------------------------------------------------------------------
def count_x_chars(sequence):
    return sequence.count('X')


# -----------------------------------------------------------------------------
def write_sequence_counts(sequence_counts, output_file):
    with open(output_file, 'w') as f:
        f.write("Seq_ID\tSequence\tCount\tX_Count\n")
        counter = 0
        for sequence, count in sequence_counts.most_common():
            counter += 1
            x_count = count_x_chars(sequence)
            f.write(f"h{counter}\t{sequence}\t{count}\t{x_count}\n")


# -----------------------------------------------------------------------------
def write_seq_counts_fasta(sequence_counts, output_file):
    with open(output_file, 'w') as f:
        counter = 0
        for sequence, count in sequence_counts.most_common():
            counter += 1
            x_count = count_x_chars(sequence)

            f.write(f">h{counter}_n{count}_x{x_count}\n")
            f.write(sequence + "\n")


# -----------------------------------------------------------------------------
def get_min_x_count(sequence_counts):
    min_X_count = min(count_x_chars(seq) for seq in sequence_counts.keys())
    return(min_X_count)


# -----------------------------------------------------------------------------
# CALL MAIN
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

