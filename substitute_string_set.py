#!/usr/bin/env python3
"""
Purpose: Substitute all occurrences of one string set with another in a file
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-02-07

Details: Substitutes all occurrences of one set of strings (e.g., sequence
    IDs) with another set of strings (e.g., phenotype data) provided in a
    second file, where the strings to be replaced serve as keys. Original
    occurrences may be kept (appended to) or replaced (substituted). One or
    more strings may be provided for the substitution, to be combined with a
    user-provided delimiter.

    If replacements are to be drawn from multiple files, the script will need
    to be run once for each file. This means that after the first replacement
    is made, subsequent replacements will not take place. This should not be
    problematic if all sources are expected to agree, i.e., to have the same
    values or combinations of values when present for a given 'find' key.

Returns:
    - STDOUT: documentation
    - FILE: --out_file is a copy of the original file with the substitution
"""

import argparse
import os
import pandas as pd
import re
import sys
from collections import defaultdict
from typing import Dict, List, NamedTuple


usage = """# -----------------------------------------------------------------------------
substitute_string_set.py - Substitute all occurrences of one string set with another in a file
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ substitute_string_set.py --help
    $ pydoc ./substitute_string_set.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ substitute_string_set.py --seq_file=big_genomes.fa
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    in_file: str
    meta_file: str
    find: str
    replace: str
    out_file: str
    delimiter: str
    exclude: str
    end_flank: str
    append: bool
    unique: bool


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Substitute all occurrences of one string set with another in a file. ' 
                    'HELP: substitute_string_set.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--in_file',
                        metavar='FILE',
                        help='Input file with values to be replaced [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)

    parser.add_argument('-m',
                        '--meta_file',
                        metavar='FILE',
                        help='Metadata file with key/value pairs to be used in the replacement; must be TAB-delimited '
                             '[REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)

    parser.add_argument('-f',
                        '--find',
                        metavar='str',
                        help='Name of column in meta_file containing the key (string to be replaced) [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)

    parser.add_argument('-r',
                        '--replace',
                        metavar='str',
                        help='Name of column(s) in meta_file containing the value(s) (strings to be inserted); '
                             'multiple comma-separated values may be provided (e.g., "species,extinct") [REQUIRED]',
                        required=True,
                        nargs=1,  # one string, but multiple CSVs allowed
                        type=str)

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of output file, i.e., a copy of the original file with the substitution [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-d',
                        '--delimiter',
                        metavar='str',
                        help='Delimiter for multiple replacement values, if provided [OPTIONAL]',
                        required=False,
                        nargs=1,
                        type=str,
                        default=',')

    parser.add_argument('-e',
                        '--exclude',
                        metavar='str',
                        help='Values to exclude from consideration as potential replacements (e.g., "NA"); '
                             'multiple comma-separated values may be provided (e.g., "NA,nan" [OPTIONAL])',
                        required=False,
                        nargs=1,
                        type=str)

    parser.add_argument('-E',
                        '--end_flank',
                        metavar='str',
                        help='Character(s) required to flank the end of the search string (e.g., "," if value of '
                             '--find must be followed by a comma) [OPTIONAL])',
                        required=False,
                        nargs=1,
                        type=str,
                        default=[''])

    parser.add_argument('-a',
                        '--append',
                        help='Append to the original values rather than replacing them, using delimiter [OPTIONAL]',
                        action='store_true')

    parser.add_argument('-u',
                        '--unique',
                        help='Replace with only unique (distinct) values, i.e., do not use repeats [OPTIONAL]',
                        action='store_true')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # DIE if any input files don't exist
    if bad_files := [file for file in [args.in_file[0], args.meta_file[0]] if not os.path.isfile(file)]:
        parser.error(f'\n### ERROR: files do not exist: {",".join(bad_files)}')

    # DIE if out_file already exists
    if os.path.isfile(args.out_file[0]):
        parser.error(f'\n### ERROR: out_file="{args.out_file[0]}" already exists')

    # DIE if --find or --replace contains invalid characters for a header
    re_NOT_header = re.compile(r'[^\w\-\s\n,]')
    if re_NOT_header.search(args.find[0]):
        parser.error(f'\n### ERROR: find="{args.find[0]}" contains invalid characters')

    if re_NOT_header.search(args.replace[0]):
        parser.error(f'\n### ERROR: replace="{args.replace[0]}" contains invalid characters')

    # make sure meta_file is readable as a pandas DataFrame
    meta_df = None
    try:
        meta_df = pd.read_table(args.meta_file[0], sep='\t')
    except pd.errors.ParserError:  # didn't conform to TAB-delimited format
        parser.error(f'\n### ERROR: meta_file="{args.meta_file[0]}" cannot be read as a TAB-delimited table')

    # make sure headers are present
    if args.find[0] not in meta_df:
        parser.error(f'\n### ERROR: find="{args.find[0]}" not column in meta_file="{args.meta_file[0]}"')

    if bad_replace := [replacer for replacer in args.replace[0].split(',') if replacer not in meta_df]:
        parser.error(f'\n### ERROR: replace="{",".join(bad_replace)}" not column(s) in meta_file="{args.meta_file[0]}"')

    # TODO: delimiter only makes sense if --append or at least two values in --replace: WARNING?

    return Args(in_file=args.in_file[0],
                meta_file=args.meta_file[0],
                find=args.find[0],
                replace=args.replace[0],
                out_file=args.out_file[0],
                delimiter=args.delimiter[0],
                exclude=args.exclude if args.exclude is None else args.exclude[0],
                end_flank=args.end_flank[0],
                append=args.append,
                unique=args.unique)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    in_file = args.in_file
    meta_file = args.meta_file
    find = args.find
    replace = args.replace
    out_file = args.out_file
    delimiter = args.delimiter
    exclude = args.exclude
    end_flank = args.end_flank
    append = args.append
    unique = args.unique

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')

    # arguments received
    print(f'LOG:in_file="{in_file}"')
    print(f'LOG:meta_file="{meta_file}"')
    print(f'LOG:find="{find}"')
    print(f'LOG:replace="{replace}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:delimiter="{delimiter}"')
    print(f'LOG:exclude="{exclude}"')
    print(f'LOG:end_flank="{end_flank}"')
    print(f'LOG:append="{append}"')
    print(f'LOG:unique="{unique}"')

    # convert replace to list
    replace_list = replace.split(',')

    # -------------------------------------------------------------------------
    # meta_file: BUILD DICTIONARY of FIND/REPLACE values
    find_replace_dl: Dict[str, List[str]] = defaultdict(list)
    meta_df = pd.read_table(meta_file, sep='\t')
    for i, row in meta_df.iterrows():
        find_value = str(meta_df.at[i, find]) + end_flank  # append end_flank ('' by default)

        if find_value in find_replace_dl:
            sys.exit(f'\n### ERROR: value "{find_value}" should be unique but appears more than once in find="{find}"')

        for replacer in replace_list:
            replace_value = str(meta_df.at[i, replacer])

            if replace_value == 'nan':
                replace_value = 'NA'

            # append end_flank ('' by default)
            replace_value = replace_value + end_flank
            find_replace_dl[find_value].append(replace_value)

    # pprint(dict(find_replace_dl))

    # -------------------------------------------------------------------------
    # UNIQUE values and EXCLUSIONS, if specified
    for find_value in find_replace_dl.keys():
        if unique:
            find_replace_dl[find_value] = sorted(list(set(find_replace_dl[find_value])))

        if exclude:
            exclusions = [x + end_flank for x in exclude.split(',')]
            find_replace_dl[find_value] = [x for x in find_replace_dl[find_value] if x not in exclusions]

    # -------------------------------------------------------------------------
    # in_file: REPLACE all occurrences of find; write to out_file
    in_fh = open(in_file, 'rt')
    out_fh = open(out_file, 'wt')

    # LOOP INFILE, REPLACE, PRINT TO OUTFILE
    num_replacements = 0
    for line in in_fh:
        for find_value in sorted(find_replace_dl.keys(), reverse=True):  # reverse examines strings before substrings
            num_replaced = line.count(find_value)
            line = line.replace(find_value, delimiter.join(find_replace_dl[find_value]))

            # line, num_replaced = re.subn(find_value, delimiter.join(find_replace_dl[find_value]), line)
            # this version was too difficult because either value not interpreted as raw

            num_replacements += num_replaced
        out_fh.write(line)

    # CLOSE files
    in_fh.close()
    out_fh.close()

    # print num_replacements
    print(f'num_replacements={num_replacements}')

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
