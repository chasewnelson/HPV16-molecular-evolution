#!/usr/bin/env python3
"""
Purpose: Determine the consensus sequence (majority allele at each site) for a MSA
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-10-15

Details: returns LOG information and a RESULT FASTA FILE containing the consensus sequence
"""

import argparse
import os
import re
import sys
from Bio import AlignIO
# from Bio.Align import AlignInfo  # for .dumb_consensus()
from collections import Counter, defaultdict
from evobioinfo import GAPS, NUCS_DEFINED, NUCS_INDETERMINATE
from typing import NamedTuple, TextIO

usage = """# -----------------------------------------------------------------------------
FASTA_to_consensus.py - Determine the consensus sequence (majority allele at each site) for a MSA
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_to_consensus.py --help
    $ pydoc ./FASTA_to_consensus.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ FASTA_to_consensus.py --aln_file my_alignment.fasta
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    aln_file: TextIO
    min_freq: float
    min_count: float
    out_file: str
    min_def_count: int
    min_def_prop: float
    exclude_seq: str
    exclude_sites: str
    ambig_char: str
    def_consensus: bool


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Determine the consensus sequence (majority allele at each site) for a MSA. HELP: FASTA_to_consensus.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--aln_file',
                        metavar='FILE',
                        help='FASTA file(s) containing a multiple sequence alignment(s) [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-m',
                        '--min_freq',
                        metavar='float',
                        help='Minimum frequency to require to consider an allele to be major [OPTIONAL]',
                        required=False,
                        type=float,
                        default=0.5)

    parser.add_argument('-c',
                        '--min_count',
                        metavar='int',
                        help='Minimum count required to determine the major allele at a site [OPTIONAL]',
                        required=False,
                        type=int,
                        default=1)

    parser.add_argument('-O',
                        '--out_file',
                        metavar='str',
                        help='Output file name [OPTIONAL]',
                        required=False,
                        type=str,
                        default='<input>_consensus.fasta')  # TODO: will have to detect this later

    parser.add_argument('-C',
                        '--min_def_count',
                        metavar='int',
                        help='Minimum number of defined (non-N/non-gap) alleles required to determine the ' + \
                             'major allele at a site [OPTIONAL]',
                        required=False,
                        type=int,
                        default=1)

    parser.add_argument('-M',
                        '--min_def_prop',
                        metavar='int',
                        help='Minimum proportion of defined (non-N/non-gap) alleles required to determine the ' + \
                             'major allele at a site [OPTIONAL]',
                        required=False,
                        type=int,
                        default=0)

    parser.add_argument('-e',
                        '--exclude_seq',
                        metavar='str',
                        help='Sequence ID(s) (comma-separated) to exclude (exact match up to first space) [OPTIONAL]',
                        required=False,
                        type=str,
                        default='')

    parser.add_argument('-E',
                        '--exclude_sites',
                        metavar='str/FILE',
                        help='String or file containing comma- or whitespace-delimited list of sites to exclude ' + \
                             '(mask as N) [OPTIONAL]',
                        required=False,
                        type=str,
                        default='')

    parser.add_argument('-a',
                        '--ambig_char',
                        metavar='str/FILE',
                        help='Character to use for ambiguous sites [OPTIONAL]',
                        required=False,
                        type=str,
                        default='N')

    parser.add_argument('-d',
                        '--def_consensus',
                        help='Limit the consensus sequence to defined characters, i.e., do not consider ambiguous '
                             'characters unless no defined characters exist at a site [OPTIONAL]',
                        action='store_true')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments
    if args.min_freq < 0.5 or args.min_freq > 1:
        parser.error(f'\n### ERROR: min_freq="{args.min_freq}" must be ≥0.5 and ≤1.0')

    # FORM the output file name if not given
    # if args.out_file == '<input>_consensus.fasta':
    #     args.out_file = os.path.splitext(args.aln_file)[0] + '_consensus.fasta'  # ??

    # FORM the output file name if not given
    if args.out_file == '<input>_consensus.fasta':
        args.out_file = os.path.splitext(args.aln_file[0].name)[0] + '_consensus.fasta'

    # ensure out_file not present
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    if args.min_def_prop < 0 or args.min_def_prop > 1:
        parser.error(f'\n### ERROR: min_def_prop="{args.min_def_prop}" must be ≥0 and ≤1.0')

    # if exclude_sites is file, convert to string
    if os.path.isfile(args.exclude_sites):
        args.exclude_sites = open(args.exclude_sites).read().rstrip()

    # die if exclude_sites contains anything besides integers, whitespace, and commas
    re_NOT_d_s_comma = re.compile(r'[^\d\s\n,]')
    if re_NOT_d_s_comma.search(args.exclude_sites):
        parser.error(f'\n### ERROR: exclude_sites="{args.exclude_sites[:20]}..." may only contain integers, whitespace, and/or commas(,)')

    return Args(aln_file=args.aln_file[0],
                min_freq=args.min_freq,
                min_count=args.min_count,
                out_file=args.out_file,
                min_def_count=args.min_def_count,
                min_def_prop=args.min_def_prop,
                exclude_seq=args.exclude_seq,
                exclude_sites=args.exclude_sites,
                ambig_char=args.ambig_char,
                def_consensus=args.def_consensus)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Taipei is bloody rainy """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    aln_fh = args.aln_file
    min_freq = args.min_freq
    min_count = args.min_count
    out_file = args.out_file
    min_def_count = args.min_def_count
    min_def_prop = args.min_def_prop
    exclude_seq = args.exclude_seq
    exclude_sites = args.exclude_sites
    ambig_char = args.ambig_char
    def_consensus = args.def_consensus

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # INITIALIZE lists from comma-separated input
    exclude_seq_list = exclude_seq.split(',')

    # SITE list for overlap checking
    # num_exclude_sites = 0
    # if exclude_sites != '':
    re_s = re.compile(r'\s+')  # whitespace
    exclude_sites = exclude_sites.replace(',', ' ')
    # exclude_sites = set(sorted(map(int, re_s.split(exclude_sites))))  # actually a set
    exclude_sites_list = re_s.split(exclude_sites)
    len_exclude_sites_list = len(exclude_sites_list)

    if len_exclude_sites_list == 1 and exclude_sites_list[0] == '':
        exclude_sites_list = []
        len_exclude_sites_list = 0
    else:
        exclude_sites_list = list(map(int, exclude_sites_list))
    exclude_sites = set(exclude_sites_list)  # actually a set - pointless to sort
    num_exclude_sites = len(exclude_sites)

    if len_exclude_sites_list != num_exclude_sites:
        sys.exit(f'\n### ERROR: Duplicates provided in --exclude_sites: exclude_sites={len_exclude_sites_list} but ' + \
                 f'only {num_exclude_sites} are unique.\n')

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # Nucleotides
    # defined_nucs = ('A', 'C', 'G', 'T', 'U')

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    # print(f'LOG:aln_file="{",".join([fh.name for fh in aln_fh])}"')
    print(f'LOG:aln_file="{aln_fh.name}"')
    print(f'LOG:min_freq="{min_freq}"')
    print(f'LOG:min_count="{min_count}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:min_def_count="{min_def_count}"')
    print(f'LOG:min_def_prop="{min_def_prop}"')
    print(f'LOG:exclude_seq="{exclude_seq}"')
    # print(f'LOG:exclude_sites="{exclude_sites}"', flush=True)
    # print(f'LOG:exclude_sites="{",".join(map(str, list(exclude_sites)))}"', flush=True)
    print(f'LOG:exclude_sites="{",".join(map(str, sorted(list(exclude_sites))))}"')
    print(f'LOG:num_exclude_sites={num_exclude_sites}')
    print(f'LOG:ambig_char={ambig_char}')
    print(f'LOG:def_consensus={def_consensus}')

    # -------------------------------------------------------------------------
    # INPUT FASTA SEQUENCE(S) as Biopython MSA
    # print(f'adding FASTA records from aln_file={aln_fh.name}')
    #recs = SeqIO.parse(aln_fh, 'fasta')
    recs = AlignIO.read(aln_fh, 'fasta')  # unit is ALIGNMENT, and .read() because we have only one ALN per file

    # CHECK FOR DUPLICATE SEQUENCE IDs - AlignIO does NOT do this
    recs_ids = [rec.id for rec in recs]
    if len(set(recs_ids)) != len(recs_ids):  # set may be LESS; this may catch more
        sys.exit(f'\n### ERROR: file={recs.name} has duplicate sequence names (IDs)')

    # # Store alignment length (number of SITES)
    # if aln_length is None:
    #     aln_length = recs.get_alignment_length()
    #     print(f'LOG:aln_length={aln_length}')
    # elif aln_length != recs.get_alignment_length():
    #     sys.exit(f'\n### ERROR: file={aln_fh.name} does not have alignment length (num sites) of {aln_length}')

    # Store original number of sequences
    nseqs_original = len(recs)

    # LOG nseqs_original
    print(f'LOG:nseqs_original={nseqs_original}')

    # Die if <2 sequences
    if nseqs_original < 2:
        sys.exit(f'\n### ERROR: nseqs={nseqs_original} must be >1 to determine a consensus')

    # recs = recs + SeqIO.parse(aln_fh, 'fasta')
    # print(f'type(recs)={type(recs)}')

    # build index-to-id mapper
    idx_to_id = defaultdict(str)
    for i, rec in enumerate(recs):
        idx_to_id[i] = rec.id

    # ELIMINATE records in reverse index order, so that coordinates don't shift
    excluded_seq_count = 0
    for idx in sorted(idx_to_id.keys(), reverse=True):
        if idx_to_id[idx] in exclude_seq_list:
            new_recs = recs[0:i]
            new_recs.extend(recs[i+1:])
            recs = new_recs
            del new_recs
            excluded_seq_count += 1

    # Store revised number of sequences
    nseqs = len(recs)

    # LOG nseqs_original
    print(f'LOG:nseqs={nseqs}')

    # Die if <2 sequences
    if nseqs < 2:
        sys.exit(f'\n### ERROR: nseqs={nseqs} after custom sequence exclusion; must be >1 to determine a consensus')

    # -------------------------------------------------------------------------
    # LOOP SITES and STORE CONSENSUS

    # Initialize consensus sequence and counters
    consensus_seq = ''
    # IUPAC_seq = ''
    excluded_site_count = 0

    for site_idx in range(recs.get_alignment_length()):
        site = site_idx + 1

        if site in exclude_sites:
            consensus_seq += ambig_char
            excluded_site_count += 1

        else:
            site_nucs = recs[:, site_idx]
            site_counts = Counter(site_nucs)

            # initialize counters
            undef_nuc_count = 0
            def_nuc_count = 0

            # count nucs; if def_consensus, delete ambiguous characters
            for nuc in site_counts.keys():
                if nuc not in NUCS_DEFINED:
                    undef_nuc_count += site_counts[nuc]

                    if def_consensus:
                        del site_counts[nuc]

                else:
                    # count DEFINED nucs
                    def_nuc_count += site_counts[nuc]

            # Determine major nuc
            maj_nuc, maj_nuc_count = site_counts.most_common(1)[0]

            # If second most common nuc has identical count, set position to undefined / AMBIGUOUS
            second_most_common = site_counts.most_common(2)
            if len(second_most_common) == 2:
                if maj_nuc_count == second_most_common[1][1]:
                    maj_nuc = ambig_char
                    # TODO: count and report ambiguous sites introduced

            # ASSIGN CONSENSUS NUC
            if def_nuc_count > 0:
                maj_nuc_freq = maj_nuc_count / def_nuc_count
                # maj_nuc_pct = round(100 * maj_nuc_freq)
            else:
                maj_nuc_freq = None

            def_nuc_prop = def_nuc_count / nseqs

            if maj_nuc_freq is not None and maj_nuc_freq >= min_freq and maj_nuc_count >= min_count and \
                    def_nuc_count >= min_def_count and def_nuc_prop >= min_def_prop:
                consensus_seq += maj_nuc
            else:
                consensus_seq += ambig_char

            # VALIDATE number of sequences
            total_nuc_count = def_nuc_count + undef_nuc_count

            if total_nuc_count != nseqs:
                sys.exit(f'\n### ERROR: Number of nucleotides not equal to number of sequences at site {site}\n')

    # -------------------------------------------------------------------------
    # PRINT SEQUENCE to OUTPUT FILE
    # group_out_file_hdl = open(f'{out_dir}_{freq}.tsv', "wt")
    out_file_hdl = open(out_file, "wt")
    out_file_hdl.write(f'>{os.path.splitext(aln_fh.name)[0] + "_consensus"}\n' +
                       consensus_seq + '\n')

    # -------------------------------------------------------------------------
    # Report excluded seq and site counts
    print(f'RESULT:excluded_site_count={excluded_site_count}')
    print(f'RESULT:excluded_seq_count={excluded_seq_count}')

    # -------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
