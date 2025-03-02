#!/usr/bin/env python3
"""
Purpose: Split records having multiple ALT values onto their own lines; LATER: convert MNV to SNV when appropriate
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-03-10
"""

import argparse
import os
import numpy as np
import re
import sys
from typing import List, NamedTuple, Optional, TextIO


usage = """# -----------------------------------------------------------------------------
VCF_splitter.py - Split records having multiple ALT values onto their own lines; LATER: convert MNV to SNV when appropriate
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ VCF_splitter.py --help
    $ pydoc ./VCF_splitter.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ VCF_splitter.py --seq_file=big_genomes.fa
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    in_file: TextIO
    out_file: str
    ALT_name: str
    INFO_name: str
    AF_key: str
    # DP_key: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Split records having multiple ALT values onto their own lines; LATER: convert MNV to SNV when ' 
                    'appropriate. HELP: VCF_splitter.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--in_file',
                        metavar='FILE',
                        help='A VCF or VCF_like TSV (TAB-separated values) tables [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of output file [REQUIRED]',
                        required=True,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-A',
                        '--ALT_name',
                        metavar='str',
                        help='Name of ALT column [OPTIONAL]',
                        required=False,
                        type=str,
                        default='ALT')

    parser.add_argument('-I',
                        '--INFO_name',
                        metavar='str',
                        help='Name of INFO column [OPTIONAL]',
                        required=False,
                        type=str,
                        default='INFO')

    parser.add_argument('-F',
                        '--AF_key',
                        metavar='str',
                        help='FORMAT key to use for ALT allele frequency [OPTIONAL]',
                        required=False,
                        type=str,
                        default='AF')

    # parser.add_argument('-d',
    #                     '--DP_key',
    #                     metavar='str',
    #                     help='FORMAT key to use for read depth (coverage) [OPTIONAL]',
    #                     required=False,
    #                     type=str,
    #                     default='DP')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # ONE input file
    if len(args.in_file) > 1:
        parser.error(f'\n### ERROR: must provide exactly 1 input file')

    # DIE if out_file already exists
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(in_file=args.in_file[0],
                out_file=args.out_file,
                ALT_name=args.ALT_name,
                INFO_name=args.INFO_name,
                AF_key=args.AF_key)  # DP_key=args.DP_key


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    in_fh = args.in_file
    out_file = args.out_file
    ALT_name = args.ALT_name
    INFO_name = args.INFO_name
    AF_key = args.AF_key
    # DP_key = args.DP_key

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')

    # arguments received
    print(f'LOG:in_file="{in_fh.name}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:ALT_name="{ALT_name}"')
    print(f'LOG:INFO_name={INFO_name}')
    print(f'LOG:AF_key={AF_key}')
    # print(f'LOG:DP_key={DP_key}')

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # match user-supplied rules
    # re_rule = re.compile(r'(\w+)([!=<>]+)([.\d]+)')

    # VCF file regex
    # NOTE: DP and AF are sufficient because the allele count, so inconsistently coded across VCFs, can simply be inferred

    # VCF_line_pattern = r"([\w\_\.]+)\t(\d+)\t([^\s]+)\t([acgtuACGTU]+)\t([acgtuACGTU]+)\t(\d+)\t(\w+)\t([\w\=\;\.\,]+)"
    # re_VCF_INFO = re.compile(VCF_INFO_pattern)

    re_VCF_INFO = re.compile(r'(\w+)=([\w\d/,.\-<>]+)')

    # -------------------------------------------------------------------------
    # OPEN FILES; SPLIT and PRINT AS YOU GO

    # open out_file
    out_fh = open(out_file, "wt")

    # prepare variables
    header_list: List[str] = []
    ALT_idx = None
    INFO_idx = None

    # initialize counters
    line_num = 0
    rec_num = 0

    # LINES
    for line in in_fh:
        line_num += 1

        # chomp newline
        line = line.rstrip()

        if line.startswith("##"):  # this is metadata
            # WRITE metadata line
            out_fh.write(line + "\n")

        elif line.startswith("#"):  # this is the header
            header_list = line[1:].split("\t")  # remove #, split

            if ALT_name not in header_list or INFO_name not in header_list:
                sys.exit(f'\n### ERROR: "{ALT_name}" or "{INFO_name} provided but are not present in the header')

            ALT_idx = header_list.index(ALT_name)
            INFO_idx = header_list.index(INFO_name)

            # WRITE header line
            out_fh.write(f'{line}\tALT_num\tVAF\n')  # out_fh.write(f'{line}\t{AF_key}\n')

        elif line_num == 1:  # first line but doesn't begin with #, so much be VCF-like table without metadata
            header_list = line.split("\t")  # remove #, split

            if ALT_name not in header_list or INFO_name not in header_list:
                sys.exit(f'\n### ERROR: "{ALT_name}" or "{INFO_name} provided but are not present in the header')

            ALT_idx = header_list.index(ALT_name)
            INFO_idx = header_list.index(INFO_name)

            # WRITE header line
            out_fh.write(f'{line}\tALT_num\tVAF\n')  # out_fh.write(f'{line}\t{AF_key}\n')
        else:  # data line
            rec_num += 1

            line_list = line.split("\t")

            # -------------------------------------------------------------
            # EXTRACT DATA for this SNP
            # this_CHROM = line_list[0]
            # this_POS = int(line_list[1])
            # this_ID = line_list[2]
            # this_REF = str(line_list[3])
            this_ALT_rec = str(line_list[ALT_idx])
            # this_ALT_rec = str(line_list[4])
            # this_QUAL = line_list[5]
            # this_FILTER = line_list[6]  # formatted as a semicolon-separated list, so not possible to CSV by variant
            this_INFO_rec = str(line_list[INFO_idx])
            # this_INFO_rec = line_list[7]
            # this_FORMAT_rec = line_list[8]
            # this_sample_rec = line_list[9]

            # -------------------------------------------------------------
            # PROCESS DATA
            # ==> ALT <==
            # this_ALT_rec_comma_n = this_ALT_rec.count(',')
            this_ALT_list = this_ALT_rec.split(',')
            # this_ALT_n = len(this_ALT_list)  # ALREADY USE this_ALT_n for another purpose

            # ==> INFO <==
            # re_this_INFOgroups = re_VCF_INFO.match(this_INFO)
            this_INFO_list = this_INFO_rec.split(';')
            this_INFO_list_keys = list(this_INFO_list)
            this_INFO_list_keys = [re_VCF_INFO.sub(r'\1', x) for x in this_INFO_list_keys]

            this_INFO_list_values: List[Optional[str]] = [re_VCF_INFO.sub(r'\2', x) if '=' in x else None for x in
                                                          list(this_INFO_list)]
            # this_INFO_list_values = list(this_INFO_list)
            # this_INFO_list_values = [re_VCF_INFO.sub(r'\2', x) if '=' in x else None for x in this_INFO_list_values]
            this_INFO_data = dict(zip(this_INFO_list_keys, this_INFO_list_values))
            # this_DP = re_VCF_DP.search(this_INFO_rec)
            # this_DP = this_INFO_rec[this_DP.start():this_DP.end()]  # AF=0.367003
            # this_DP = int(this_DP.replace("DP=", ""))
            # this_INFO_DP = this_INFO_data[DP_key]  # TODO: option later?
            # this_AF_rec = re_VCF_AF.search(this_INFO_rec)
            # this_AF_rec = this_INFO_rec[this_AF_rec.start():this_AF_rec.end()]  # AF=0.367003; AF=0,1
            # this_AF_rec = this_AF_rec.replace('AF=', '')
            this_INFO_AF = this_INFO_data[AF_key]  # TODO: option later?

            # print(f'this_INFO_rec={this_INFO_rec}')
            # print(f'this_INFO_list_keys={this_INFO_list_keys}')
            # print(f'this_INFO_list_values={this_INFO_list_values}')
            # print(f'this_INFO_data={this_INFO_data}')
            # print(f'this_INFO_DP={this_INFO_DP}')
            # print(f'this_INFO_AF={this_INFO_AF}')

            this_AF_list: List[float] = list(map(float, this_INFO_AF.split(',')))
            # this_AF_list = this_AF_rec.split(',')
            # this_AF_list = list(map(float, this_AF_list))

            # print(f'this_DP={this_DP}')
            # # print(f'type(this_DP)={type(this_DP)}')
            # print(f'this_ALT_list={this_ALT_list}')
            # print(f'this_AF_rec={this_AF_rec}')
            # print(f'this_AF_list={this_AF_list}')
            # print(f'this_AC_list={this_AC_list}')

            if len(this_ALT_list) != len(this_AF_list):
                sys.exit('\n### ERROR: number of ALT alleles and AF values do not match for the following record:\n'
                         f'{line}\n')

            # -----------------------------------------------------------------
            # print one line for each ALT/AF
            ALT_num = 0
            for temp_ALT, temp_AF in zip(this_ALT_list, this_AF_list):
                ALT_num += 1
                temp_line_list = line_list
                temp_line_list[ALT_idx] = temp_ALT
                temp_line = '\t'.join(temp_line_list) + f'\t{ALT_num}\t{temp_AF}'

                # WRITE info line for this one ALT/AF
                out_fh.write(f'{temp_line}\n')

    # close output file
    out_fh.close()

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
