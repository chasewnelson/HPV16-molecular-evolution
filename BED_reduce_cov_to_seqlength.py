#!/usr/bin/env python3
"""
Purpose: For BED files with ranges wrapping the circular sequence end, reduce coverage values to sequence length
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
from collections import defaultdict
from typing import List, NamedTuple

usage = """# -----------------------------------------------------------------------------
BED_reduce_cov_to_seqlength.py - For BED files with ranges wrapping the circular sequence end, reduce coverage values to sequence length
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ BED_reduce_cov_to_seqlength.py --help
    $ pydoc ./BED_reduce_cov_to_seqlength.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ BED_reduce_cov_to_seqlength.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    bed_dir: str  # TextIO
    seq_length: int
    chrom: str
    out_suffix: str
    bed_ext: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='For BED files with ranges wrapping the circular sequence end, reduce coverage values to '
                    'sequence length. HELP: BED_reduce_cov_to_seqlength.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-b',
                        '--bed_dir',
                        metavar='DIR',
                        help='Directory containing BED file(s) named for the sequence(s) in the FASTA [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)  # argparse.FileType('rt'))

    parser.add_argument('-s',
                        '--seq_length',
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
    # OPTIONAL - NONE
    parser.add_argument('-o',
                        '--out_suffix',
                        metavar='str',
                        help='Name of FASTA output file to contain masked sequences [REQUIRED]',
                        required=False,
                        type=str,
                        default='reduced')

    parser.add_argument('-e',
                        '--bed_ext',
                        metavar='str',
                        help='Extension used by BED files in the bed_dir directory [OPTIONAL]',
                        required=False,
                        type=str,
                        default='.bed')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # BED dir exists
    if not os.path.isdir(args.bed_dir[0]):
        parser.error(f'\n### ERROR: bed_dir="{args.bed_dir[0]}" does not exist')

    # for each file in bed_dir, ensure out_file implied by out_suffix not present
    bed_file_list = [file for file in os.listdir(args.bed_dir[0]) if file.endswith(args.bed_ext)]
    for bed_file in bed_file_list:
        # form output
        out_file = os.path.splitext(bed_file)[0] + f'_{args.out_suffix}{args.bed_ext}'
        # print(f'file={bed_file} out_file={out_file}')

        if os.path.isfile(out_file):
            parser.error(f'\n### ERROR: out_file="{out_file}" already exists')

    return Args(bed_dir=args.bed_dir,
                seq_length=args.seq_length,
                chrom=args.chrom,
                out_suffix=args.out_suffix,
                bed_ext=args.bed_ext)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    bed_dir = args.bed_dir[0]
    seq_length = args.seq_length
    chrom = args.chrom
    out_suffix = args.out_suffix
    bed_ext = args.bed_ext

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
    print(f'LOG:bed_dir="{bed_dir}"')
    print(f'LOG:seq_length="{seq_length}"')
    print(f'LOG:chrom="{chrom}"')
    print(f'LOG:out_suffix="{out_suffix}"')
    print(f'LOG:bed_ext={bed_ext}')

    # -------------------------------------------------------------------------
    # Make a dictionary of BED file IDs to file names - UNNECESSARY given the current method of file naming

    # -------------------------------------------------------------------------
    # LOOP BED FILES and COMPUTE WRAP AROUND COVERAGES
    bed_file_list = [file for file in os.listdir(bed_dir) if file.endswith(bed_ext)]
    nfiles = 0

    # FILES
    for bed_file in bed_file_list:
        nfiles += 1

        # form OUTPUT BED file
        out_file = os.path.splitext(bed_file)[0] + f'_{args.out_suffix}{args.bed_ext}'
        # bed_file_name = os.path.join(bed_dir, rec.id + bed_suffix)

        # INITIALIZE DICT
        # bed_file_dict = dict()

        # OPEN BED
        bed_recs = HTSeq.BED_Reader(bed_file)
        rec_num = 0
        max_iv_end = 0

        # prepare start / end lists
        iv_start_list: List[int] = []
        iv_end_list: List[int] = []
        iv_coverage_list: List[int] = []

        # LOOP RECS
        for bed_rec in bed_recs:
            # bed_rec_gff = bed_rec.get_gff_line()
            this_chrom = str(bed_rec.iv.chrom)
            this_coverage = int(bed_rec.name)

            if this_chrom == chrom:
                rec_num += 1  # only those pertinent to our chosen chrom
                # cov = int(bed_rec.name)

                # process range, splitting in two if spans the circular end
                iv_start = int(bed_rec.iv.start)  # 0-based
                iv_end = int(bed_rec.iv.end)
                # iv_length = int(bed_rec.iv.length)

                # split into segments corresponding to the reference sites
                if iv_end > 2 * seq_length:  # error if double wrap
                    sys.exit(f'\n### ERROR: can accomodate one wrap-around, but iv_end={iv_end} is > (2 * seq_length)')

                elif iv_start >= seq_length:  # end index 7906) for pos 7906
                    # EXAMPLE
                    # POS 7907..7969, IDX 7906..7968, BED [7906,7968)
                    # POS 1..62, IDX 0..61, BED [0,62)
                    # start - length is 7906-7906 = 0, CORRECT FOR BED
                    # end - length is 7968-7906 = 62, CORRECT FOR BED
                    iv_start_list.append(iv_start - seq_length)
                    iv_end_list.append(iv_end - seq_length)
                    iv_coverage_list.append(this_coverage)

                    # replace max if current is bigger
                    if iv_end > max_iv_end:
                        max_iv_end = iv_end - seq_length

                elif iv_end > seq_length:  # only end but not start is over length; must be split into TWO segments
                    iv_start_list.extend([iv_start, 0])  # inclusive index
                    iv_end_list.extend([seq_length, iv_end - seq_length])  # non-inclusive index
                    iv_coverage_list.extend([this_coverage, this_coverage])

                    # replace max if current is bigger
                    if iv_end > max_iv_end:
                        max_iv_end = iv_end - seq_length

                else:  # no overlap end, works as-is
                    # TAKE 1
                    iv_start_list.append(iv_start)  # inclusive index
                    iv_end_list.append(iv_end)  # non-inclusive index
                    iv_coverage_list.append(this_coverage)

        # LOOP STORED INTERVALS, split again into before / after the max_end
        iv_start_list_mod = []
        iv_end_list_mod = []
        iv_coverage_list_mod = []

        # initialize dict for pre-max values
        pre_max_cov_dict = defaultdict(int)

        # if before the max_end, store each site separately for ease
        for iv_start, iv_end, iv_coverage in zip(iv_start_list, iv_end_list, iv_coverage_list):
            # print(f'file={bed_file} iv_start={iv_start} iv_end={iv_end} iv_coverage={iv_coverage}')

            if iv_start >= max_iv_end:  # whole interval is after max_iv_end, store as-is
                iv_start_list_mod.append(iv_start)
                iv_end_list_mod.append(iv_end)
                iv_coverage_list_mod.append(iv_coverage)

            elif iv_end > max_iv_end:  # only end but not start is at or over length; must be split into TWO segments
                # pre-max segment: store each site separately
                for i in range(iv_start, max_iv_end):
                    # print(f'\ti={i}')
                    pre_max_cov_dict[i] += iv_coverage

                # post-max segment: store as-is
                iv_start_list_mod.append(max_iv_end)  # inclusive index
                iv_end_list_mod.append(iv_end)  # non-inclusive index
                iv_coverage_list_mod.append(iv_coverage)  # non-inclusive index

            else:  # whole interval before max_iv_end, save all site-by-site
                for i in range(iv_start, iv_end):
                    # print(f'\ti={i}')
                    pre_max_cov_dict[i] += iv_coverage

        # REDUCE runs of identical values? Should be sorted, but can use as-is ordering for now - BUG, FORGET IT
        pre_max_cov_dict_indices = sorted(pre_max_cov_dict.keys())
        pre_max_cov_start_list = []
        pre_max_cov_end_list = []
        pre_max_cov_coverage_list = []
        for idx in pre_max_cov_dict_indices:
            idx_coverage = pre_max_cov_dict[idx]

            end_idx = idx + 1

            while end_idx in pre_max_cov_dict_indices:
                if pre_max_cov_dict[end_idx] == idx_coverage:
                    del pre_max_cov_dict[end_idx]
                    pre_max_cov_dict_indices.remove(end_idx)
                    end_idx += 1
                else:
                    break

            # append values
            pre_max_cov_start_list.append(idx)
            pre_max_cov_end_list.append(end_idx)
            pre_max_cov_coverage_list.append(idx_coverage)

        # # ADD site-by-site from before max_iv_end to lists - JUST USING DICT, NO REDUCING RUNS
        # for idx in sorted(pre_max_cov_dict.keys()):
        #     iv_start_list_mod.append(idx)
        #     iv_end_list_mod.append(idx + 1)
        #     iv_coverage_list_mod.append(pre_max_cov_dict[idx])

        # OPEN OUT_FILE amd PRINT
        out_fh = open(out_file, 'wt')

        # FIRST PRINT BEFORE or REDUCED data from before max_iv_end to lists
        for iv_start, iv_end, iv_coverage in zip(pre_max_cov_start_list, pre_max_cov_end_list, pre_max_cov_coverage_list):
            # # if first only
            # iv_start_list_mod.append(iv_start)
            # iv_end_list_mod.append(iv_end)
            # iv_coverage_list_mod.append(iv_coverage)

            # if REDUCED
            out_fh.write(f'{chrom}\t{iv_start}\t{iv_end}\t{iv_coverage}\n')

        # NEXT PRINT data from AFTER max_iv_end
        for iv_start, iv_end, iv_coverage in zip(iv_start_list_mod, iv_end_list_mod, iv_coverage_list_mod):
            out_fh.write(f'{chrom}\t{iv_start}\t{iv_end}\t{iv_coverage}\n')

        out_fh.close()

    # -------------------------------------------------------------------------
    # LOG sequence counts and groups
    # Number of sequences
    print(f'LOG:nfiles={nfiles}')

    # -------------------------------------------------------------------------
    # TIME & DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_time_per_seq = elapsed_time / nfiles

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time / 60, ndigits=1)} minutes; {round(elapsed_time / 60 / 60, ndigits=1)} hours)')

    print(f'TIME ELAPSED PER SEQ: {round(elapsed_time_per_seq, ndigits=9)} seconds')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
