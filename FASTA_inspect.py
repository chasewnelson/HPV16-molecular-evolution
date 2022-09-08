#!/usr/bin/env python3
"""
Purpose: Inspect a FASTA file for basic metrics of interest and lack of redundancy
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-01-29

Details: <FILL IN>

Returns:
    <FILL IN>
"""

import argparse
import os
import sys
from Bio import Align, AlignIO, SeqIO
from collections import Counter, defaultdict
from evobioinfo import GAPS, hamming, IUPAC, IUPAC_AMBIG, NUCS_DEFINED, NUCS_INDETERMINATE, summary_string
from numpy import nan as NA
from typing import Dict, List, NamedTuple


# Usage
usage = """# -----------------------------------------------------------------------------
FASTA_inspect.py - Inspect a FASTA file for basic metrics of interest and lack of redundancy
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_inspect.py --help
    $ pydoc ./FASTA_inspect.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_inspect.py --seq_file=seqs.fasta -p
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    seq_file: str
    p_dist: bool


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Inspect a FASTA file for basic metrics of interest and lack of redundancy. HELP: FASTA_inspect.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--seq_file',
                        metavar='FILE',
                        help='FASTA file containing sequence(s) [REQUIRED]',
                        required=True,
                        # nargs='+',
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-p',
                        '--p_dist',
                        help='Activate with alignments to calculate p-distance between sequences with identical IDs '
                             '[OPTIONAL]',
                        action='store_true')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments
    if not os.path.isfile(args.seq_file):
        parser.error(f'\n### ERROR: seq_file="{args.seq_file}" does not exist')

    return Args(seq_file=args.seq_file,
                p_dist=args.p_dist)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    seq_file = args.seq_file
    p_dist = args.p_dist

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT AND LOG
    print(usage)

    print('# -----------------------------------------------------------------------------')
    print('LOG:')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')

    # arguments received
    print(f'LOG:seq_file="{seq_file}"')
    print(f'LOG:p_dist="{p_dist}"')

    # -------------------------------------------------------------------------
    # REGEX & TUPLES
    # re_date = re.compile(r'(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+).(\d+)')

    # -------------------------------------------------------------------------
    # Determine whether US=unaligned sequences or MSA=multiple sequence alignment
    recs = None  # TODO determine whether this is necessary
    seq_file_type = 'US'
    MSA_length = NA
    nseqs = 0
    try:
        recs = AlignIO.read(seq_file, 'fasta')  # Bio.Align.MultipleSeqAlignment
        seq_file_type = 'MSA'  # only gets here if you don't throw error
        MSA_length = int(recs.get_alignment_length())
        nseqs = len(recs)
    except ValueError:  # thrown if AlignIO used for collection of sequences of differing lengths
        recs = SeqIO.parse(seq_file, 'fasta')  # Bio.SeqIO.FastaIO.FastaIterator

        # count sequences
        for rec in recs:
            nseqs += 1

    # change Bio.Align.MultipleSeqAlignment to Bio.SeqIO.FastaIO.FastaIterator if only one sequence
    if isinstance(recs, Align.MultipleSeqAlignment) and len(recs) == 1:
        recs = SeqIO.parse(seq_file, 'fasta')
        seq_file_type = 'US'  # change back

    # print seq file type
    print(f'LOG:seq_file_type="{seq_file_type}" ("US"=unaligned sequences; "MSA"=multiple sequence alignment)')
    print(f'LOG:MSA_length={MSA_length} sites')
    print(f'LOG:nseqs={nseqs} seqs')
    # print(type(recs))

    if p_dist and seq_file_type != 'MSA':
        print('\n### WARNING: ACTIVATED -p/--p_dist for unaligned sequence(s): WILL SKIP')

    # -------------------------------------------------------------------------
    # REOPEN as appropriate file type
    if seq_file_type == 'US':
        recs = SeqIO.parse(seq_file, 'fasta')
    else:
        recs = AlignIO.read(seq_file, 'fasta')

    # -------------------------------------------------------------------------
    # Print sequence characters and counts
    print('\n# -----------------------------------------------------------------------------')
    print('SEQUENCE CHARACTERS (*SUSPICIOUS):', flush=True)
    suspicious_rec_count = 0
    suspicious_counts_dict: Dict[str, int] = defaultdict(int)
    seq_length_dict: Dict[str, int] = defaultdict(int)
    gap_ambig_rec_count = 0
    gap_ambig_props_dict: Dict[str, int] = defaultdict(int)

    # initialize print keys
    seq_char_keys_list: List[str] = []
    seq_char_keys_list.extend(list(NUCS_DEFINED))
    seq_char_keys_list.extend(list(NUCS_INDETERMINATE))
    seq_char_keys_list.extend(list(GAPS))
    seq_char_keys_list.extend(sorted(list(set(IUPAC_AMBIG).difference(NUCS_INDETERMINATE))))
    # print(f'seq_char_keys_list={seq_char_keys_list}')

    # PRINT 'header'
    seq_char_header = '{:<5} {:<11}{:<9}'.format('num', 'name', 'length')
    seq_char_header += '{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}'.format(*seq_char_keys_list)
    seq_char_header += '{:>11}'.format('GAP_AMBIG')
    print(seq_char_header)

    # LOOP recs
    rec_num = 0
    for rec in recs:
        rec_num += 1  # first rec will be 1, 1-based
        seq_length = len(rec.seq)
        counts: Dict[str, int] = Counter(rec.seq)
        # counts_list: List[str] = []
        warning = False

        seq_length_dict[rec.id] = seq_length

        # max_digits = len(str(max(counts.values())))
        # print(f'max_digits={max_digits}')

        # PRINT the counts of recognized IUPAC characters
        print('{:<5} {:<11}{:<9}'.format(rec_num, rec.id, seq_length), end='')
        for seq_char in seq_char_keys_list:
            # print(f'{nuc}={counts[nuc]}', end='\t')
            print('{:>9d}'.format(counts[seq_char]), end='')

        gap_ambig_sum = 0
        for seq_char in sorted(counts.keys()):
            # counts_list.append(f'{char}={counts[char]}')
            if seq_char in GAPS or seq_char in IUPAC_AMBIG:
                gap_ambig_sum += counts[seq_char]

            if seq_char not in GAPS and seq_char not in IUPAC:
                warning = True
                suspicious_rec_count += 1
                suspicious_counts_dict[seq_char] += counts[seq_char]

        # print GAP_AMBIG percent and newline
        print('{:>11}'.format(str(round(100 * gap_ambig_sum / seq_length, 2)) + '%'), end='')
        print()

        gap_ambig_props_dict[rec.id] = gap_ambig_sum / seq_length
        if gap_ambig_sum > 0:
            gap_ambig_rec_count += 1

        if warning:
            # print(f'{rec.id}: ' + ','.join(counts_list))
            print(' <= *')

    # -------------------------------------------------------------------------
    # SEQUENCE SUMMARY, including SUSPICIOUS characters and counts
    print('\n# -----------------------------------------------------------------------------')
    print('SEQUENCE SUMMARY:', flush=True)
    # print(f'LOG:mean_seq_length={np.mean(list(seq_length_dict.values()))}')
    print(f'seq_length: {summary_string(list(seq_length_dict.values()))}')
    if suspicious_rec_count == 0:
        print('NO SUSPICIOUS CHARACTERS')
    else:
        print(f'### WARNING: {suspicious_rec_count} records with SUSPICIOUS CHARACTERS')
        for seq_char in sorted(suspicious_counts_dict.keys()):
            print(f'"{seq_char}"={suspicious_counts_dict[seq_char]}')

    # # -------------------------------------------------------------------------
    # # Print GAP and AMBIGUOUS counts
    # print('\n# -----------------------------------------------------------------------------')
    # print('GAPS and AMBIGUITIES:', flush=True)
    print(f'Records with gaps or ambiguities: {gap_ambig_rec_count} ({round(100 * gap_ambig_rec_count / nseqs, 2)}%)')
    print(f'gaps_ambiguities: {summary_string([round(100 * x, 2) for x in list(gap_ambig_props_dict.values())], suffix="%")}')

    # -------------------------------------------------------------------------
    # Detect and print duplicated IDs
    print('\n# -----------------------------------------------------------------------------')
    print('DUPLICATED SEQUENCE IDS:', flush=True)

    ID_counts: Dict[str, int] = defaultdict(int)
    for rec in recs:
        ID_counts[str(rec.id)] += 1

    # print(dict(ID_counts))

    # dict comprehension for repeats
    ID_repeats = {ID: count for ID, count in ID_counts.items() if count > 1}

    if ID_repeats:
        sorted_repeat_IDs = sorted(ID_repeats.keys())
        print(f'LOG:num_dup_IDs={len(sorted_repeat_IDs)}')
        print(f'LOG:dup_IDs="{",".join(sorted_repeat_IDs)}"')
        for dup in sorted_repeat_IDs:
            print(f'{dup} ({ID_repeats[dup]} seqs)')
    else:
        print('NONE')

    # -------------------------------------------------------------------------
    # Calculate p-distance between sequences having the same name
    print('\n# -----------------------------------------------------------------------------')
    print('P-DISTANCE BETWEEN DUPLICATED SEQUENCE IDS:', flush=True)

    if not p_dist:
        print('NOT ACTIVATED (-p/--p_dist not called)')
    elif seq_file_type != 'MSA':
        print('\n### WARNING: ACTIVATED -p/--p_dist for unaligned sequences (US): SKIPPING')
    elif ID_repeats:
        for dup_ID in sorted_repeat_IDs:
            dup_ID_seq_list: List[str] = []
            # p_dist_list: List[float] = []
            p_dist_string_list: List[str] = []
            if seq_file_type == 'MSA':
                recs = AlignIO.read(seq_file, 'fasta')
            else:
                recs = SeqIO.parse(seq_file, 'fasta')

            dup_ID_count = 0
            for rec in recs:
                if rec.id == dup_ID:
                    dup_ID_seq_list.append(rec.seq)
                    dup_ID_count += 1

            for i, seq1 in enumerate(dup_ID_seq_list):  # range(len(dup_ID_seq_list)):
                for j in range(i + 1, len(dup_ID_seq_list)):
                    seq2 = dup_ID_seq_list[j]
                    diffs = hamming(seq1, seq2)
                    prop_diff = diffs / MSA_length
                    # p_dist_list.append(prop_diff)
                    p_dist_string_list.append(f'{prop_diff} ({diffs} diffs)')

            print(f'{dup_ID} ({dup_ID_count} seqs) distances: {",".join(map(str, p_dist_string_list))}', flush=True)
    else:
        print('\nACTIVATED -p/--p_dist but there are NO DUPLICATE NAMES (a good thing!)')

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
