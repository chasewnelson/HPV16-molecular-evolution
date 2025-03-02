#!/usr/bin/env python3
"""
Purpose: Compute the mean values across all numeric cells for multiple identically formatted tables
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-03-10
"""

import argparse
import os
import numpy as np
import pandas as pd
import sys
from typing import NamedTuple, TextIO


usage = """# -----------------------------------------------------------------------------
mean_of_tables.py - Compute the mean values across all numeric cells for multiple identically formatted tables
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ mean_of_tables.py --help
    $ pydoc ./mean_of_tables.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ mean_of_tables.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    in_files: TextIO
    out_file: str
    key: str
    exclude: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Compute the mean values across all numeric cells for multiple identically formatted tables. ' 
                    'HELP: mean_of_tables.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--in_files',
                        metavar='FILE',
                        help='Two or more identically formatted TSV (TAB-separated values) tables [REQUIRED]',
                        required=True,
                        nargs='+',
                        type=argparse.FileType('rt'))

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of output file to contain table of mean values [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-k',
                        '--key',
                        metavar='str',
                        help='Name of column to use as the identifying KEY; auto-excluded from calculation [REQUIRED]',
                        required=True,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-e',
                        '--exclude',
                        metavar='str',
                        help='Column names to exclude; '
                             'multiple should be comma-separated (e.g., "name,color"). '
                             'Failure to exclude non-numeric columns will result in an ERROR [REQUIRED]',
                        required=False,
                        nargs=1,
                        type=str)

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # At least TWO input files
    if len(args.in_files) < 2:
        parser.error(f'\n### ERROR: must provide at least 2 input files, with names separated by spaces')

    # DIE if out_file already exists
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    # convert exclude to list with key
    exclude_list = [args.key]
    if args.exclude is not None:
        exclude_list.append(args.exclude[0].split(','))
    print(f'exclude_list={exclude_list}')

    # verify
    # 1) files readable as a pandas DataFrame
    # 2) files have identical headers
    # 3) files have identical shapes
    # 4) all non-excluded columns are NUMERIC
    header_list = None
    file_shape = None
    for fh in args.in_files:
        try:
            file_df = pd.read_table(fh, sep='\t')

            # check for identical shape
            if file_shape is None:
                file_shape = file_df.shape
            elif not file_shape == file_df.shape:
                parser.error(f'\n### ERROR: files do not have identical numbers of rows and columns: {file_shape} vs. '
                             f'{file_df.shape}')

            # check for identical header names
            if header_list is None:
                header_list = list(file_df.columns)
            elif not header_list == list(file_df.columns):
                parser.error(f'\n### ERROR: files do not use identical header names')

            # check for non-numeric columns
            for colname in header_list:
                if colname not in exclude_list:  # and colname != args.key:
                    if not pd.api.types.is_numeric_dtype(file_df[colname].dtypes):
                        parser.error(f'\n### ERROR: column is not numeric: file="{fh.name}", '
                                     f'colname={colname}, dtype={file_df[colname].dtypes}')

        except pd.errors.ParserError:  # didn't conform to TAB-delimited format
            parser.error(f'\n### ERROR: file="{fh.name}" cannot be read by pandas as a TAB-delimited table')

    # # make sure headers are present
    # if args.find[0] not in meta_df:
    #     parser.error(f'\n### ERROR: find="{args.find[0]}" not column in meta_file="{args.meta_file[0]}"')
    #
    # if bad_replace := [replacer for replacer in args.replace[0].split(',') if replacer not in meta_df]:
    #     parser.error(f'\n### ERROR: replace="{",".join(bad_replace)}" not column(s) in meta_file="{args.meta_file[0]}"')
    #
    # # TODO: delimiter only makes sense if --append or at least two values in --replace: WARNING?
    #
    # if re_NOT_header.search(args.replace[0]):
    #     parser.error(f'\n### ERROR: replace="{args.replace[0]}" contains invalid characters')

    return Args(in_files=args.in_files,
                out_file=args.out_file,
                key=args.key,
                exclude=args.exclude if args.exclude is None else args.exclude[0])


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    in_fhs = args.in_files
    out_file = args.out_file
    key = args.key
    exclude = args.exclude

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')

    # arguments received
    print(f'LOG:in_files="{",".join([fh.name for fh in in_fhs])}"')
    print(f'LOG:in_files_count={len([fh.name for fh in in_fhs])}')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:key="{key}"')
    print(f'LOG:exclude={exclude}')

    # convert exclude to list with key
    exclude_list = [args.key]
    if args.exclude is not None:
        exclude_list.append(args.exclude.split(','))
    print(f'exclude_list={exclude_list}')

    # -------------------------------------------------------------------------
    # LOOP FILES
    nfiles = 0

    # INITIALIZE dictionary
    col_row_list_dict = dict()  # colname -> rowname -> [list of values, 1 from each file]

    # FILES
    for fh in in_fhs:
        nfiles += 1
        filename = fh.name  # because already opened once?
        file_df = pd.read_table(filename, sep='\t')
        # file_df = pd.read_table(fh, sep='\t')

        if nfiles == 1:
            col_names = list(file_df.columns)

        # initialize keys
        if nfiles == 1:
            key_list = [str(item) for item in file_df[key]]

            for key_name in key_list:
                col_row_list_dict[key_name] = dict()

        # COLUMNS
        # ncols = 0
        for colname in file_df.columns:
            if colname not in exclude_list:
                # column = file_df[colname]
                # ncols += 1

                # if nfiles == 1:
                #     col_row_list_dict[colname] = dict()

                # ROWS
                for i in range(len(file_df[colname])):
                    key_name = file_df[key][i]
                    # print(f'nfiles={nfiles} colname={colname} i={i} key={key_name} save: {file_df[colname][i]}')  # ncols={ncols}
                    # print(f'key={key_name}')
                    if nfiles == 1:
                        # col_row_list_dict[key_name] = dict()
                        # print('here1')
                        # col_row_list_dict[colname] = dict()
                        col_row_list_dict[key_name][colname] = [file_df[colname][i]]
                    else:
                        # print('here2')
                        col_row_list_dict[key_name][colname].append(file_df[colname][i])

    # LOG number of files process
    print(f'LOG:files_processed={nfiles}')

    # -------------------------------------------------------------------------
    # COMPUTE and PRINT MEANS

    # STORE column, row (key) names for a consistent order
    row_names = []
    for row_name in col_row_list_dict.keys():
        row_names.append(row_name)

    col_names = []
    for col_name in col_row_list_dict[row_names[0]].keys():  # just use first row name
        col_names.append(col_name)

    # open out_file
    out_fh = open(out_file, 'wt')

    # form and write header
    header_list = [key]
    header_list.extend(col_names)
    out_fh.write('\t'.join(header_list) + '\n')

    # write values
    for row_name in row_names:  # col_row_list_dict.keys():
        this_line = [row_name]

        for col_name in col_names:  # col_row_list_dict[row_name].keys():
            this_col_list = col_row_list_dict[row_name][col_name]
            this_line.append(np.mean(this_col_list))
            # print(f'\nthis_list={this_col_list}\nthis_list_mean={np.mean(this_col_list)}\n')
        # print('\t'.join([str(item) for item in this_line]))
        out_fh.write('\t'.join([str(item) for item in this_line]) + '\n')

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
