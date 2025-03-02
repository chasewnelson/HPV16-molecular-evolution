#!/usr/bin/env python3
"""
Purpose: Process raxml-ng log files for ML tree information
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-11-27

Details: Process raxml-ng log files for ML tree information
"""

import argparse
import os
import re
import sys
import time
from collections import defaultdict
from typing import NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
raxml-ng-log-process.py - Process raxml-ng log files for ML tree informations
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ raxml-ng-log-process.py --help
    $ pydoc ./raxml-ng-log-process.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ raxml-ng-log-process.py --in_file
# -----
"""


# -----------------------------------------------------------------------------
class Args(NamedTuple):
    """ Command-line arguments """
    in_file: TextIO
    out_file: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Process raxml-ng log files for ML tree information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--in_file',
                        metavar='FILE',
                        help='raxml-ng .log file(s) [REQUIRED]',
                        required=True,
                        nargs='+',
                        type=argparse.FileType('rt'))

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Output file name [OPTIONAL]',
                        required=False,
                        type=str,
                        default='<in_file>.TABLE.tsv')  # detect this later

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # FORM the output file name if not given
    if args.out_file == '<in_file>.TABLE.tsv':
        # concatenate in_file base names
        in_file_names = [os.path.splitext(os.path.basename(this_fh.name))[0] for this_fh in args.in_file]

        # get file name extension
        in_file_extension = os.path.splitext(os.path.basename(args.in_file[0].name))[1]

        # form output file name
        out_file_name = '_'.join(in_file_names) + in_file_extension + '.TABLE.tsv'
        args.out_file = out_file_name

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(in_file=args.in_file,
                out_file=args.out_file)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    in_fh_list = args.in_file
    out_file = args.out_file

    # form in_file name list
    in_file_name_list = [os.path.basename(this_fh.name) for this_fh in args.in_file]
    in_file_name_list_str = ','.join(in_file_name_list)

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # REGEX & TUPLES
    # re_rule = re.compile(r'(\w+)([!=<>]+)([.\d]+)')
    re_tree_log = re.compile(r'\[(\d\d):(\d\d):(\d\d)\] \[worker #(\d+)\] ML tree search #(\d+), logLikelihood: (-\d*.\d*)')
    # '[00:02:11] [worker #13] ML tree search #14, logLikelihood: -23385.588027'

    # re_start_trees_both = re.compile(r'start tree\(s\): random \((\d+)\) \+ parsimony \((\d+)\)')
    re_start_trees_rand = re.compile(r'start tree\(s\): random \((\d+)\)')
    # re_start_trees_pars = re.compile(r'start tree\(s\): parsimony \((\d+)\)')

    re_random_seed = re.compile(r'random seed: (\d+)')

    re_model = re.compile(r'^Model: (.+)$')

    # -------------------------------------------------------------------------
    # Record the current working directory
    curr_dir = str(os.getcwd())

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{curr_dir}"')
    print(f'LOG:in_file="{in_file_name_list_str}"')
    print(f'LOG:out_file="{out_file}"')

    # -------------------------------------------------------------------------
    # PREPARE OUTPUT FILE
    out_file_hdl = open(out_file, "wt")

    # header
    header_list = ['file_name', 'search', 'seed', 'model', 'start_tree_type',
                   'time', 'hour', 'min', 'sec', 'worker', 'logLikelihood']
    header = '\t'.join(header_list)
    out_file_hdl.write(header + '\n')

    # -------------------------------------------------------------------------
    # LOOP INPUT FILE(S)
    for log_fh in in_fh_list:
        log_fh_name = log_fh.name

        # ---------------------------------------------------------------------
        # LOOP FILE LINES, searching for regex and saving

        # EXAMPLE: # '[00:02:11] [worker #13] ML tree search #14, logLikelihood: -23385.588027'
        search_to_info_dl = defaultdict(list)
        rand_trees = 0
        # pars_trees = 0
        random_seed = 'NA'
        model = 'NA'

        for line in log_fh:
            # chomp newline
            line = line.rstrip()
            # print(f'line:"{line}"')

            if re_tree_log_MATCH := re_tree_log.match(line):
                TREE_LOG_MATCH_GROUPS = re_tree_log_MATCH.groups()
                log_time = f'{TREE_LOG_MATCH_GROUPS[0]}:{TREE_LOG_MATCH_GROUPS[1]}:{TREE_LOG_MATCH_GROUPS[2]}'
                log_hour = TREE_LOG_MATCH_GROUPS[0]
                log_min = TREE_LOG_MATCH_GROUPS[1]
                log_sec = TREE_LOG_MATCH_GROUPS[2]
                log_worker = TREE_LOG_MATCH_GROUPS[3]
                log_search = int(TREE_LOG_MATCH_GROUPS[4])
                log_LL = TREE_LOG_MATCH_GROUPS[5]
                search_to_info_dl[log_search] = [log_time, log_hour, log_min, log_sec, log_worker, log_LL]
            # elif re_start_trees_both_SEARCH := re_start_trees_both.search(line):
            #     START_TREES_BOTH_SEARCH_GROUPS = re_start_trees_both_SEARCH.groups()
            #     rand_trees = int(START_TREES_BOTH_SEARCH_GROUPS[0])
            #     pars_trees = int(START_TREES_BOTH_SEARCH_GROUPS[1])
            elif re_start_trees_rand_SEARCH := re_start_trees_rand.search(line):
                START_TREES_RAND_SEARCH_GROUPS = re_start_trees_rand_SEARCH.groups()
                rand_trees = int(START_TREES_RAND_SEARCH_GROUPS[0])
            # elif re_start_trees_pars_SEARCH := re_start_trees_pars.search(line):
            #     START_TREES_PARS_SEARCH_GROUPS = re_start_trees_pars_SEARCH.groups()
            #     pars_trees = int(START_TREES_PARS_SEARCH_GROUPS[0])
            elif re_random_seed_SEARCH := re_random_seed.search(line):
                RANDOM_SEED_SEARCH_GROUPS = re_random_seed_SEARCH.groups()
                random_seed = RANDOM_SEED_SEARCH_GROUPS[0]
            elif re_model_MATCH := re_model.match(line):
                MODEL_MATCH_GROUPS = re_model_MATCH.groups()
                model = MODEL_MATCH_GROUPS[0]

        # search_type_list: List[str] = []
        # for i in range(rand_trees):
        #     num = i + 1
        #     search_type_list[num] = 'rand'

        # for i in range(pars_trees):
        #     num = i + 1
        #     search_type_list[num] = 'pars'

        # -------------------------------------------------------------------------
        # PRINT TO OUT_FILE
        for search_num in sorted(search_to_info_dl.keys()):

            # Determine starting tree type
            start_tree_type = 'rand'
            if search_num > rand_trees:
                start_tree_type = 'pars'

            # print
            out_line_list = [log_fh_name, search_num, random_seed, model, start_tree_type]
            out_line_list.extend(search_to_info_dl[search_num])
            out_line = '\t'.join(map(str, out_line_list))
            # print(line)
            out_file_hdl.write(out_line + '\n')

    out_file_hdl.close()

    # -------------------------------------------------------------------------
    # DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time/60, ndigits=1)} minutes; {round(elapsed_time/60/60, ndigits=1)} hours)')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')

# --------------------------------------------------
if __name__ == '__main__':
    main()
