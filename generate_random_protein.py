#! /usr/bin/env python3

"""
Purpose: Generate randomized peptides from an input amino acid sequence
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2023-06/28
"""

import argparse
import os
import random
import sys
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
generate_random_protein.py - Generate randomized peptides from an input amino acid sequence
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ generate_random_protein.py --help
    $ pydoc ./generate_random_protein.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ generate_random_protein.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    in_file: TextIO
    out_file: str
    length: int
    number: int
    prefix: str  # optional


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Generate randomized peptides from an input amino acid sequence. HELP: generate_random_protein.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--in_file',
                        metavar='FILE',
                        help='Peptide multiple sequence alignment (FASTA) [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of output file to contain randomized peptides (e.g., "input_random_l9_n1000.fasta" '
                             '[REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-l',
                        '--length',
                        metavar='int',
                        help='Length of randomized peptides to generate'
                             '[REQUIRED]',
                        required=True,
                        type=int)

    parser.add_argument('-n',
                        '--number',
                        metavar='int',
                        help='Number of randomized peptides to generate [REQUIRED]',
                        required=True,
                        type=int)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-p',
                        '--prefix',
                        metavar='str',
                        help='Prefix for FASTA ids [OPTIONAL]',
                        required=False,
                        type=str,
                        default='')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    # Both integers are >0
    if args.length < 1:
        parser.error(f'\n### ERROR: length={args.length} is <1')
    if args.number < 1:
        parser.error(f'\n### ERROR: number={args.number} is <1')

    return Args(in_file=args.in_file[0],
                out_file=args.out_file,
                length=args.length,
                number=args.number,
                prefix=args.prefix)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    proteome_file = args.in_file
    out_file = args.out_file
    peptide_length = args.length
    num_peptides = args.number
    prefix = args.prefix

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

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
    print(f'LOG:in_file="{proteome_file.name}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:length={peptide_length}')
    print(f'LOG:number={num_peptides}')
    print(f'LOG:prefix="{prefix}"')

    # -------------------------------------------------------------------------
    # OPEN PROTEIN FASTA
    recs = SeqIO.parse(proteome_file, "fasta")
    #print(recs)

    proteome = ""
    for rec in recs:
        #print(rec.seq)
        this_seq = str(rec.seq)
        if str(rec.seq).endswith('*'):
            print(f'Removing trailing termination symbol * from: {this_seq}')
            this_seq = this_seq.rstrip('*')

        # add seq to proteome
        proteome = proteome + str(this_seq)

    #proteome[:10]
    #'MESLVPGFNE'

    #proteome = list(proteome)
    #proteome[:10]
    #['M', 'E', 'S', 'L', 'V', 'P', 'G', 'F', 'N', 'E']

    #random.shuffle(proteome) # acts in situ!
    #proteome[:10]
    #['A', 'F', 'A', 'K', 'K', 'L', 'Y', 'I', 'T', 'P']

    print(f'LOG:proteome="{proteome}"')
    print(f'LOG:proteome_length="{len(proteome)}"')

    # Produce a FASTA of random peptides using the proteome
    # Following: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc284
    with open(out_file, "w") as outfile_hdl:
        for i in range(num_peptides):
            # generate this randomized peptide
            #random.shuffle(proteome)
            this_random_peptide = random.choices(proteome, k=peptide_length)

            # create as record
            shuffled_rec = SeqRecord(
                Seq("".join(this_random_peptide)),  # rec.seq.alphabet deprecated as of Biopython version 1.78
                id=f'{prefix}s{i + 1}',  # id="s%i" % (i + 1),
                description=""  # "random peptide drawn from %s" % proteome_file.name
            )

            # write record to file (appends)
            outfile_hdl.write(shuffled_rec.format("fasta"))

    # -------------------------------------------------------------------------
    # LOG sequence counts, unnecessary here
    # print(f'LOG:nseqs={nseqs}')

    # -------------------------------------------------------------------------
    # TIME & DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time / 60, ndigits=1)} minutes; {round(elapsed_time / 60 / 60, ndigits=1)} hours)')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

