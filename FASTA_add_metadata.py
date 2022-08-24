#!/usr/bin/env python3
"""
Purpose: Add tabular metadata to FASTA sequence headers
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-01-06

Details: Outputs a modified FASTA file with metadata added to headers

IN: /Users/cwnelson88/Desktop/NCI/research/HPV16
DO: FASTA_add_metadata.py --seq_file=seq/HPV16_PAP_20200813.N-30.fasta --meta_file=seq_metadata/Lisa_hpv16meth.20190723.tsv
"""

import argparse
import os
import pandas as pd
import re
import sys
from Bio import SeqIO
from collections import defaultdict
from datetime import datetime
from numpy import nan as NA
from typing import Dict, List, NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
FASTA_add_metadata.py - Add tabular metadata to FASTA sequence headers
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_add_metadata.py --help
    $ pydoc ./FASTA_add_metadata.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_add_metadata.py --seq_file=PAP.fasta --meta_file=metadata.tsv --join_key=MERGE_ID --out_file=PAP_meta.fasta > FASTA_add_metadata_PAP.out
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    seq_fh: TextIO
    meta_file: TextIO
    join_key: str
    out_file: str
    meta_label: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='PROGRAM: Add tabular metadata to FASTA sequence headers. HELP: FASTA_add_metadata.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    parser.add_argument('-i',
                        '--seq_file',
                        metavar='FILE',
                        help='FASTA file containing a multiple sequence alignment [REQUIRED]',
                        required=True,
                        nargs=1,  # makes it a list
                        type=argparse.FileType('rt'))

    parser.add_argument('-m',
                        '--meta_file',
                        metavar='FILE',
                        help='TSV file containing metadata [REQUIRED]',
                        required=True,
                        type=argparse.FileType('rt'))

    parser.add_argument('-j',
                        '--join_key',
                        metavar='str',
                        help='Column name from --meta_file containing key for data joining, expected at start of '
                             'FASTA header. May provide multiple comma-separated columns names to be searched in '
                             'order. Value cannot be "nan" or "NA" [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Name of output FASTA file, with metadata added to headers [REQUIRED]',
                        required=True,
                        type=str)

    parser.add_argument('-l',
                        '--meta_label',
                        metavar='str',
                        help='Label to be used for the new text block in the FASTA header [OPTIONAL]',
                        required=False,
                        type=str,
                        default='METADATA')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # DIE if output file already exists
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    return Args(seq_fh=args.seq_file[0],
                meta_file=args.meta_file,
                join_key=args.join_key,
                out_file=args.out_file,
                meta_label=args.meta_label)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    seq_fh = args.seq_fh
    meta_fh = args.meta_file
    join_key = args.join_key
    out_file = args.out_file
    meta_label = args.meta_label

    # -------------------------------------------------------------------------
    # REGEX

    # define regex for matching *BAD* PAP IDs, to be reported and excluded if encountered
    # regex_validID = re.compile(r'([ACDIPRSX]+)(\d+)_(\w+)')  # when there was a '*_HPV16' suffix to rec.id
    regex_validID = re.compile(r'([ACDIPRSX]+)(\d+)')
    # IRC: IRC204014_HPV16, ...
    # PAP: PAP0110_HPV16, ...
    # SBX|SCD: SBX1353_HPV16, SCD4176_HPV16, ...
    regex_date = re.compile(r'(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+).(\d+)')

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    date_string = str(datetime.now())
    date_string = regex_date.sub(r'\1-\2-\3', date_string)

    # Print output
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:seq_file="{seq_fh.name}"')
    print(f'LOG:meta_file="{meta_fh.name}"')
    print(f'LOG:join_key="{join_key}"')
    print(f'LOG:out_file="{out_file}"', flush=True)

    # -------------------------------------------------------------------------
    # IMPORT the HPV16 sequences
    recs = SeqIO.parse(seq_fh, 'fasta')

    # -------------------------------------------------------------------------
    # IMPORT HPV sequence metadata (one row => one sequence's metadata)
    meta_df = pd.read_table(meta_fh, sep='\t')
    print(f'LOG:meta_file_nrows={len(meta_df)}')

    """
               PID   MERGE_ID         PAP_ID  HPV16_SET  PATIENT_ID                   HPVSTRING  ... UPDATED_CASE_STAT CYTO_DX_CGR HISTOLOGY_CGR  CC_CGR TIME_TO_CASE  CLEARED
    0     37025735  PAP100004  PAP1000040030          2      236518  16,31,45,53,58,62,68,71,72  ...                 1           2           NaN     0.0          0.0      0.0
    1     31142965  PAP100016  PAP1000160030          2      691185                          16  ...                 1           0           NaN     0.0          0.0      0.0
    2     65728770  PAP100209  PAP1002090030          2      606979                          16  ...                 3           5           NaN     1.0         34.0      NaN
    3     30470610  PAP100223  PAP1002230030          2      404283                       16,66  ...                 1           1           NaN     0.0          0.0      0.0
    4     38235675  PAP100309  PAP1003090030          2      254310                     6,16,53  ...                 3           1           NaN     1.0       1425.0      NaN
    ...        ...        ...            ...        ...         ...                         ...  ...               ...         ...           ...     ...          ...      ...
    3574  67588670  PAP603828  PAP6038280050          2     1037979              16,31,52,54,61  ...                 0           0           NaN     0.0        256.0      1.0
    3575  62589300  PAP604398  PAP6043980050          2     1023646                          16  ...                 2           1           NaN     1.0         43.0      NaN
    3576  50601285  PAP604547  PAP6045470050          2      405059                          16  ...                 2           1           NaN     1.0         26.0      NaN
    3577  20766695    PAP2596  PAP6059930050          1     1031688                          16  ...                 3           3           NaN     1.0         59.0      NaN
    3578   5294620  PAP606638  PAP6066380050          2      405908                          16  ...                 0           0           NaN     0.0          0.0      0.0
    """

    # FORM INDICES based on column names provided in input
    join_key_list = join_key.split(',')
    # print(f'join_key_list={join_key_list}')
    indices_list: List[str] = []
    join_key_count_dict: Dict[str, int] = defaultdict(int)

    for i, row in meta_df.iterrows():
        found_key = False

        for key in join_key_list:
            key_value = meta_df.at[i, key]

            # If it's found, add it to indices and go to next row
            if pd.notnull(key_value) and key_value is not None and key_value != NA and key_value != '' \
                    and key_value != 'nan' and key_value != 'NA':
                # print(f'found key_value={key_value}')
                # print(f'type(key_value)={type(key_value)}')
                found_key = True
                indices_list.append(key_value)
                join_key_count_dict[key] += 1
                break

        if not found_key:
            sys.exit(f'\n### ERROR: no acceptable --join_key={join_key_list} found in metadata at row {i}')

    for key in join_key_list:
        print(f'LOG:times_used={join_key_count_dict[key]} for join_key="{key}"')
    # print(f'len(meta_df)={len(meta_df)}')
    # print(f'len(indices_list)={len(indices_list)}')

    # SET INDICES
    # meta_df.index = meta_df[join_key]  # this only necessary if using .at[] rather than .loc[] below?
    meta_df.index = indices_list  # this only necessary if using .at[] rather than .loc[] below?

    # -------------------------------------------------------------------------
    # CREATE new column containing metadata string for FASTA header

    # # -------------------------------------------------------------------------
    # # CREATE NEW COLUMN containing metadata string with semicolon-separated KEY=VALUE pairs
    # # -------------------------------------------------------------------------
    # # SOLUTION 1 FROM Twitter @least_recent (Jonathan)
    # # (1) df.insert()  <==  inserts column into DataFrame; parameters are (0-based loc to insert, new col name, Series)
    # # (2) len(df.columns)  <==  gets the number of columns (also works: df.shape[1])
    # # (3) 'meta_string'  <==  name of new col
    # # (4) []  <=  a list, which is an array-like object interpreted as a Series
    # # (5) df.at[i, c]  <==  retrieves a single value in a DataFrame, here the value at row number i, col name c
    # # (6) range(len(df))  <==  integers from 0 to number of rows
    # # (7) [<list comp>]  <==  creates semicolon-separated NAME=VALUE for all columns, once for each row i
    # meta_df.insert(len(meta_df.columns),
    #                'meta_string',
    #                [';'.join([f'{c}={meta_df.at[i, c]}' for c in meta_df.columns]) for i in range(len(meta_df))])
    # # run times: 1.914, 1.615, 1.688, 1.848 seconds

    # -------------------------------------------------------------------------
    # CREATE NEW COLUMN containing metadata string with semicolon-separated KEY=VALUE pairs
    # -------------------------------------------------------------------------
    # SOLUTION 2 FROM Twitter @guan (Guan)
    # (1) df.items()  <==  returns a generator of tuples with (colname, Series)
    # (2) Series.items()  <==  returns a generator of tuples with (index, value) [PERHAPS similar to df.iterrows()]]
    # (3) df.apply()  <==  applies a function across rows (axis=0 or 'rows') or across columns (axis=1 or 'columns')
    # (4) lambda x thus operates once for each row, across columns (axis=1)
    # (5) each row is a Series, so calling .items() in the lambda function gets all the col=value pairs
    meta_df['meta_string'] = meta_df.apply(lambda x: ';'.join(f'{k}="{v}"' for k, v in x.items()), axis='columns')
    # run times: 1.293, 0.859, 0.849, 1.213 seconds

    # # My own BULLSHIT
    # # meta_df['meta_string'] = ''
    # # for col_name in meta_df.columns:
    # #     print(f'col_name="{col_name}"')
    # #     # print(f'meta_df[col_name]="{meta_df[col_name]}"')  # outputs 1-column DF with indices
    # #     # print(f'meta_df[col_name].astype(str)="{meta_df[col_name].astype(str)}"')
    # #     # print(f'meta_df[col_name].values="{meta_df[col_name].values}"')
    # #     print(f'list(meta_df[col_name].values)="{list(meta_df[col_name].values)}"')
    #
    # for col_name in meta_df.columns:
    #     meta_df['meta_string'] += list(map(str, meta_df[col_name].values))
    #
    # # for col_name in ['MERGE_ID', 'ENRL_AGE', 'HPVSTRING']:
    # #     meta_df['meta_string'] += f'{col_name}={meta_df[col_name].astype(str)};'
    #
    # # # HPV16_metadata['meta_string'] = f'MERGE_ID=' + HPV16_metadata['MERGE_ID'].astype(str) + ';' + \
    # #     f'HPVSTRING={HPV16_metadata["HPVSTRING"].astype(str)};' + \
    # #     f'ENRL_AGE={HPV16_metadata["ENRL_AGE"].astype(str)};' + \
    # #     f'RACE_ETH=' + HPV16_metadata['RACE_ETH'].astype(str) + ';' + \
    # #     f'STUDY_WORST_CYTO_DX=' + HPV16_metadata['STUDY_WORST_CYTO_DX'].astype(str) + ';' + \
    # #     f'WORST_BX_STUDY=' + HPV16_metadata['WORST_BX_STUDY'].astype(str) + ';' + \
    # #     f'UPDATED_CASE_STAT=' + HPV16_metadata['UPDATED_CASE_STAT'].astype(str) + ';' + \
    # #     f'TIME_TO_CASE=' + HPV16_metadata['TIME_TO_CASE'].astype(str) + ';' + \
    # #     f'CLEARED=' + HPV16_metadata['CLEARED'].astype(str)

    # -------------------------------------------------------------------------
    # SAVE updated DataFrame - only for debugging
    # meta_df.to_csv('seq_metadata/Lisa_hpv16meth.20190723_wMetadataString.tsv', sep='\t', index=True)

    # -------------------------------------------------------------------------
    # JOIN METADATA to FASTA headers and PRINT
    # OPEN out_file for writing, to which SeqIO.write() will append
    out_fh = sys.stdout
    if out_file != '':
        out_fh = open(out_file, 'wt')

    # initialize counters
    nseqs = 0
    nseqs_wMeta = 0
    nseqs_woMeta = 0
    recs_present: List[str] = []
    recs_missing: List[str] = []
    # recs_updated = []
    for rec in recs:  # .id's of the form 'PAP0016_HPV16', 'PAP0037_HPV16', ...
        nseqs += 1
        # print(rec)  # Number of features: 0
        # print(f'rec.id="{rec.id}"')  # rec.id="PAP0016_HPV16"
        # print(f'rec.name="{rec.name}"')  # rec.name="PAP0016_HPV16"
        # print(f'rec.description="{rec.description}"')  # rec.description="PAP0016_HPV16"
        # print(f'rec.seq="{rec.seq}"')  # Seq('ACTACAATAATTCATGTATAAAACTAAGGGTGTAACCGAAATCGGTTGAACCGA...TAA')

        # ID_LIST = rec.id.split('_')
        # # print(f'ID_LIST="{ID_LIST}"')  # ID_LIST="['PAP237467', 'HPV16']"
        # ID = str(ID_LIST[0])
        # HPV_type = str(ID_LIST[1])

        # kill if unexpected format; otherwise extract ID
        if not regex_validID.search(rec.id):  # len(ID_LIST) != 2 or HPV_type != 'HPV16'  # regex_nonPAPID.search(ID):
            sys.exit(f"\n### TERMINATED: UNEXPECTED ID={rec.id}")

        if rec.id in meta_df.index:
            nseqs_wMeta += 1
            recs_present.append(rec.id)
            new_description = f'{meta_label}:source="{os.path.basename(meta_fh.name)}";date="{date_string}";{str(meta_df.at[rec.id, "meta_string"])}'  # assumes ID is the INDEX
            if rec.description == rec.id:
                rec.description = new_description
            else:
                rec.description = rec.description + ' ' + new_description
            # rec.description = f'{meta_label}}:source="{os.path.basename(meta_fh.name)}";date="{date_string}";{meta_df.loc[meta_df[join_key] == ID, "meta_string"].values[0]}'  # .values is a 1-elt list
            # LESSON: modifying will sadly NOT change the rec in the original FastaIterator
        else:
            nseqs_woMeta += 1
            # print(f'ID="{ID}" missing from metadata="{meta_fh.name}"')
            recs_missing.append(rec.id)
        # KeyError: 'PAP0097'

        # rec.id = ID + ' ' + str(meta_df.loc[meta_df['MERGE_ID'] == ID, 'meta_string'])  # WORKS but FORMAT!!!

        # rec.id = ID + ' ' + HPV16_metadata.get_value(ID, 'metadata_string')  # ERROR: 'DataFrame' object has no attribute 'get_value'
        # rec.id = ID + ' ' + str(HPV16_metadata.loc[HPV16_metadata['MERGE_ID'] == ID, 'metadata_string'].values[0])  # ERROR: index 0 is out of bounds for axis 0 with size 0
        # rec.id = ID + ' ' + str(HPV16_metadata.loc[HPV16_metadata['MERGE_ID'] == ID, 'metadata_string'].get_value())  # ERROR: 'Series' object has no attribute 'get_value'
        # rec.id = ID + ' ' + HPV16_metadata.loc[HPV16_metadata['MERGE_ID'] == ID, 'metadata_string'].get_value())  # ???

        # rec.id = ID + ' ' + str(meta_df.at(meta_df['MERGE_ID'] == ID, 'meta_string'])  # WORKS but FORMAT!!!

        # PRINT rec
        # recs_updated.append(rec)
        SeqIO.write(rec, out_fh, 'fasta')  # APPENDS to an open filehandle

    print(f'LOG:seq_file_nseqs={nseqs}')  # 3579
    print(f'LOG:seq_file_nseqs_withMetadata={nseqs_wMeta} ({round(100 * nseqs_wMeta / nseqs, 1)}%)')
    print(f'LOG:seq_file_nseqs_withoutMetadata={nseqs_woMeta} ({round(100 * nseqs_woMeta / nseqs, 1)}%)')
    print(f'LOG:recs_with_metadata={",".join(recs_present)}')
    # PRINT records missing from FASTA file
    if len(recs_missing) > 0:
        print(f'\n### WARNING: meta_file={meta_fh.name} is missing data for n={len(recs_missing)} sequences present in seq_file={seq_fh.name}:')

        # for i, ID in enumerate(recs_missing):
        #     print(f'{i + 1}={ID}')

        print(','.join(recs_missing))

    out_fh.close()  # close to flush

    # # SAVE FASTA with updated metadata-containing headers
    # # SeqIO.write(recs, 'seq/HPV16_PAP_20200813.N-30_mod1.fasta', 'fasta')
    # SeqIO.write(recs_updated, out_file, 'fasta')
    # # SeqIO.write(recs_updated, 'seq/HPV16_PAP_20200813.N-30_wMetadata.fasta', 'fasta')

    # -----------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
