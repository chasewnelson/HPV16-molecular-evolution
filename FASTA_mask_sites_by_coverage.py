#!/usr/bin/env python3
"""
Purpose: Mask sites in a FASTA file based on coverage values from a BED file
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2023-03-06

Details: the COVERAGE VALUES must be present in the NAME attribute, i.e., directly following the interval
"""

import argparse
import HTSeq
import os
import sys
import time
from Bio import AlignIO
from typing import NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
FASTA_mask_sites_by_coverage.py - Mask sites in a FASTA file based on coverage values from a BED file
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_mask_sites_by_coverage.py --help
    $ pydoc ./FASTA_mask_sites_by_coverage.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_mask_sites_by_coverage.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    seq_file: TextIO
    # bed_file: str  # TextIO
    bed_dir: str  # TextIO
    bed_suffix: str
    out_file: str
    min_cov: int
    chrom: str
    mask_char: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Mask sites in a FASTA file based on coverage values from a BED file. HELP: FASTA_mask_sites_by_coverage.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-s',
                        '--seq_file',
                        metavar='FILE',
                        help='FASTA file containing one or more sequences with IDs matching BED files [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    # parser.add_argument('-b',
    #                     '--bed_file',
    #                     metavar='FILE',
    #                     help='BED file(s) containing coverage data for the sequence(s) in the FASTA [REQUIRED]',
    #                     required=True,
    #                     nargs='+',
    #                     type=str)  # argparse.FileType('rt'))

    parser.add_argument('-b',
                        '--bed_dir',
                        metavar='DIR',
                        help='Directory containing BED file(s) named for the sequence(s) in the FASTA [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)  # argparse.FileType('rt'))

    parser.add_argument('-e',
                        '--bed_suffix',
                        metavar='str',
                        help='EXT/suffix of bed files to be excluded from file names to yield sequence IDs [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of FASTA output file to contain masked sequences [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-m',
                        '--min_cov',
                        metavar='int',
                        help='Minimum coverage (≥1) to require for a site, else it will be masked [REQUIRED]',
                        required=True,
                        type=int)

    parser.add_argument('-c',
                        '--chrom',
                        metavar='str',
                        help='Name of chromosome to use for coverage data [REQUIRED]',
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

    # Minimum coverage acceptable
    if args.min_cov < 1:
        parser.error(f'\n### ERROR: min_cov="{args.min_cov}" must be ≥1')

    # # All bed files exist
    # for this_file in args.bed_file:
    #     if not os.path.isfile(this_file):
    #         parser.error(f'\n### ERROR: bed_file="{this_file}" does not exist')

    # BED dir exists
    if not os.path.isdir(args.bed_dir[0]):
        parser.error(f'\n### ERROR: bed_dir="{args.bed_dir[0]}" does not exist')

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(seq_file=args.seq_file[0],
                # bed_file=args.bed_file,
                bed_dir=args.bed_dir[0],
                bed_suffix=args.bed_suffix,
                out_file=args.out_file,
                min_cov=args.min_cov,
                chrom=args.chrom,
                mask_char=args.mask_char)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    seq_fh = args.seq_file
    # bed_file_list = args.bed_file
    bed_dir = args.bed_dir
    bed_suffix = args.bed_suffix
    out_file = args.out_file
    min_cov = args.min_cov
    chrom = args.chrom
    mask_char = args.mask_char

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    #    # -------------------------------------------------------------------------
    #    # INITIALIZE lists from comma-separated input
    #    group_key_list = group_key.split(',')
    #    exclude_aln_list = exclude_seq.split(',')
    #    exclude_group_list = exclude_group.split(',')
    #
    #    # SITE list for overlap checking
    #    # num_custom_sites = 0
    #    # if custom_sites != '':
    #    re_s = re.compile(r'\s+')  # whitespace
    #    custom_sites = custom_sites.replace(',', ' ')
    #    # custom_sites = set(sorted(map(int, re_s.split(custom_sites))))  # actually a set
    #    custom_site_list = re_s.split(custom_sites)
    #    len_custom_site_list = len(custom_site_list)
    #
    #    if len_custom_site_list == 1 and custom_site_list[0] == '':
    #        custom_site_list = []
    #        len_custom_site_list = 0
    #    else:
    #        custom_site_list = list(map(int, custom_site_list))
    #    custom_sites = set(custom_site_list)  # actually a set - pointless to sort
    #    num_custom_sites = len(custom_sites)
    #
    #    if len_custom_site_list != num_custom_sites:
    #        sys.exit(f'\n### ERROR: Duplicates provided in --custom_sites: custom_sites={len_custom_site_list} but ' + \
    #                 f'only {num_custom_sites} are unique.\n')

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
    # # custom bed file string if there are a lot
    # bed_file_count = len(bed_file_list)
    # bed_file_str = ''
    # if bed_file_count > 10:
    #     bed_file_str = ",".join([bed_file for bed_file in bed_file_list[:10]]) + f' and {bed_file_count - 10} others'
    # else:
    #     bed_file_str = ",".join([bed_file for bed_file in bed_file_list])

    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:seq_file="{seq_fh.name}"')
    print(f'LOG:bed_dir="{bed_dir}"')
    # print(f'LOG:bed_file="{bed_file_str}"')
    # print(f'LOG:bed_file="{",".join([bed_fh.name for bed_fh in bed_fh_list])}"')
    print(f'LOG:bed_suffix="{bed_suffix}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:min_cov={min_cov}')
    print(f'LOG:chrom="{chrom}"')
    print(f'LOG:mask_char="{mask_char}"')
    # print(f'RESULT:bed_file_count={bed_file_count}')

    # -------------------------------------------------------------------------
    # Make a dictionary of BED file IDs to file names - UNNECESSARY given the current method of file naming

    # -------------------------------------------------------------------------
    # TODO: CHECK FOR DUPLICATE SEQUENCE IDs

    # -------------------------------------------------------------------------
    # OPEN FASTA FILE; LOOP SEQS; for each, FIND ITS BED AND MASK GIVEN COVERAGE

    # get recs - using AlignIO for now because it IS an MSA
    recs = AlignIO.read(seq_fh, 'fasta')  # unit is ALIGNMENT, and .read() because we have only one ALN per file

    # Keep track of alignment length
    seq_length = recs.get_alignment_length()  # None
    nseqs = 0
    # nseqs_wMeta = 0
    # nseqs_woMeta = 0

    # loop recs and make sure each has an associated BED file
    for rec in recs:

        # DIE if a BED file does not exist
        bed_file_name = os.path.join(bed_dir, rec.id + bed_suffix)
        # print(f'bed_file_name="{bed_file_name}"')
        if not os.path.isfile(bed_file_name):
            sys.exit(f'\n### ERROR: bed_file="{bed_file_name}" does not exist')

#        for site, nuc in enumerate(rec.seq, start=1):
#            recs_woMeta.append(rec.id)
#           # print(f'### WARNING: no group metadata in header of {rec.id}')

    # FASTA: loop recs and MASK
    # rec_counter = 0
    for rec in recs:
        nseqs += 1

        # rec_counter += 1
        # if rec_counter == 10:
        #     break

        # form BED file
        bed_file_name = os.path.join(bed_dir, rec.id + bed_suffix)
        # print(f'seq={}:bed="{bed_file_name}"')

        # BED: open file
        bed_recs = HTSeq.BED_Reader(bed_file_name)

        for bed_rec in bed_recs:
            # bed_rec_gff = bed_rec.get_gff_line()
            this_chrom = bed_rec.iv.chrom

            if this_chrom == chrom:
                # cov = int(bed_rec.name)

                # if cov too low, MASK THIS STRETCH
                if int(bed_rec.name) < min_cov:
                    iv_start = int(bed_rec.iv.start)  # 0-based
                    iv_end = int(bed_rec.iv.end)
                    # iv_length = int(bed_rec.iv.length)

                    # prepare if overlapping amplicons
                    iv_start_list = []
                    iv_end_list = []

                    if iv_end > 2 * seq_length:
                        sys.exit(f'\n### ERROR: can accomodate wrap-around, but iv_end={iv_end} is > (2 * seq_length)')
                    elif iv_start >= seq_length:  # TODO double check not > | end index 7906) for pos 7906
                        # EXAMPLE
                        # POS 7907..7969, IDX 7906..7968, BED [7906,7968)
                        # POS 1..62, IDX 0..61, BED [0,62)
                        # start - length is 7906-7906 = 0, CORRECT FOR BED
                        # end - length is 7968-7906 = 62, CORRECT FOR BED
                        iv_start_list = [iv_start - seq_length]
                        iv_end_list = [iv_end - seq_length]
                    elif iv_end > seq_length:  # only end but not start is over length; must be split into TWO segments
                        iv_start_list = [iv_start, 0]  # inclusive index
                        iv_end_list = [seq_length, iv_end - seq_length]  # non-inclusive index
                    else:
                        # TAKE 1
                        iv_start_list = [iv_start]  # inclusive index
                        iv_end_list = [iv_end]  # non-inclusive index

                    # MASK - starts are inclusive, ends are not
                    for start, end in zip(iv_start_list, iv_end_list):
                        seg_length = end - start

                        # EXAMPLE START
                        # POS 3..5, BED [2,5)
                        # Before this range for masking, POS 1..2, the last exclusive POS==3, IDX==2, range [:2)

                        # EXAMPLE END
                        # POS 1..62, BED [0,62)
                        # After this range for masking, the next inclusive POS==63, IDX==62, same range [62:)
                        rec.seq = rec.seq[:start] + (seg_length * mask_char) + rec.seq[end:]

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
