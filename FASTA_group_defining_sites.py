#!/usr/bin/env python3
"""
--- MODULE docstring as follows [DELETE THIS LINE] ---
Purpose: Determine group-defining variants for a range of within-group variant frequencies
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-01-09
"""

import argparse
# import dendropy
import math
import os
import re
import sys
from Bio import AlignIO  # , SeqIO
from collections import defaultdict  # Counter,
from ete3 import Tree
from ete3.parser import newick
# from numpy import nan as NA
from pprint import pprint
from typing import Any, Dict, List, NamedTuple, TextIO  # Set,

usage = """# -----------------------------------------------------------------------------
FASTA_group_defining_sites.py - Determine group-defining variants for a range of within-group variant frequencies
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ FASTA_group_defining_sites.py --help
    $ pydoc ./FASTA_group_defining_sites.py
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    aln_file: TextIO
    group_key: str
    min_freq: float
    max_freq: float
    step_size: float
    out_dir: str
    out_file: str
    min_def_count: int
    min_def_prop: float
    exclude_seq: str
    exclude_group: str
    custom_sites: str
    tree: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Determine group-defining variants for a range of within-group variant frequencies. HELP: FASTA_group_defining_sites.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-i',
                        '--aln_file',
                        metavar='FILE',
                        help='FASTA file(s) containing a multiple sequence alignment and headers with group ' + \
                             'metadata [REQUIRED]',
                        required=True,
                        nargs='+',
                        type=argparse.FileType('rt'))  # seq/HPV16_PAP_20200813.N-30.fasta

    parser.add_argument('-g',
                        '--group_key',
                        metavar='str',
                        help='Data key name containing unique name of group, expected within FASTA header. ' + \
                             'May provide multiple comma-separated keys to be searched in ' + \
                             'order [REQUIRED]',
                        required=True,
                        type=str)  # sublineage

    parser.add_argument('-m',
                        '--min_freq',
                        metavar='float',
                        help='Minimum frequency (≥0.5) to require to consider a variant to be group-defining (e.g., ' + \
                             '0.8) [REQUIRED]',
                        required=True,
                        type=float)

    parser.add_argument('-M',
                        '--max_freq',
                        metavar='float',
                        help='Maximum frequency (≤1.0) to require to consider a variant to be group-defining ' + \
                             '(e.g., 1.0) [REQUIRED]',
                        required=True,
                        type=float)

    parser.add_argument('-s',
                        '--step_size',
                        metavar='float',
                        help='Step size to use between --min_freq and --max_freq (e.g., 0.05). If 0, step size ' + \
                             'will be set to 1/nseqs (all empirical frequencies in the range examined) [REQUIRED]',
                        required=True,
                        type=float)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-o',
                        '--out_dir',
                        metavar='str',
                        help='Output directory name for variant data at group-defining sites (1 file per ' + \
                             'major frequency cutoff examined) [OPTIONAL]',
                        required=False,
                        type=str,
                        default='group_defining_sites')

    parser.add_argument('-O',
                        '--out_file',
                        metavar='str',
                        help='Output file prefix for reporting number of group-defining sites at each major allele ' + \
                             'frequency cutoff examined) [OPTIONAL]',
                        required=False,
                        type=str,
                        default='group_defining_site_counts_from<#>_to<#>_by<#>.tsv')

    parser.add_argument('-c',
                        '--min_def_count',
                        metavar='int',
                        help='Minimum number of defined (non-N/non-gap) alleles required to consider a ' + \
                             'variant group-defining at a site [OPTIONAL]',
                        required=False,
                        type=int,
                        default=1)

    parser.add_argument('-p',
                        '--min_def_prop',
                        metavar='int',
                        help='Minimum proportion of defined (non-N/non-gap) alleles required to consider a ' + \
                             'variant group-defining at a site [OPTIONAL]',
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
                        '--exclude_group',
                        metavar='str',
                        help='Group ID(s) (comma-separated) to exclude [OPTIONAL]',
                        required=False,
                        type=str,
                        default='')

    parser.add_argument('-C',
                        '--custom_sites',
                        metavar='str/FILE',
                        help='String or file containing comma- or whitespace-delimited list of sites to test for ' + \
                             'overlap with identified group-defining sites [OPTIONAL]',
                        required=False,
                        type=str,
                        default='')  # Smith_HPV16_lineage_def_sites.txt with n=162

    parser.add_argument('-t',
                       '--tree',
                       metavar='str',
                       help='Tree topology describing group relationships (Newick format) [OPTIONAL]',
                       required=False,
                       type=str,
                       default='')  # t = Tree('(((((C1,C3),C4),(D4,(D1,(D3,D2)))),(B1,B4)),(A4,(A3,(A2,A1))));')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments
    if args.min_freq < 0.5 or args.min_freq > 1:
        parser.error(f'\n### ERROR: min_freq="{args.min_freq}" must be ≥0.5 and ≤1.0')

    if args.max_freq < 0.5 or args.max_freq > 1:
        parser.error(f'\n### ERROR: max_freq="{args.max_freq}" must be ≥0.5 and ≤1.0')

    if args.max_freq < args.min_freq:  # CAN be equal if only testing one value
        parser.error(f'\n### ERROR: max_freq="{args.max_freq}" must be ≥ min_freq={args.min_freq}')

    if args.step_size < 0 or args.step_size > 0.5:
        parser.error(f'\n### ERROR: step_size="{args.step_size}" must be in the range [0,0.5]')

    # ensure out_dir not present, create it
    if os.path.isdir(args.out_dir):
        parser.error(f'\n### ERROR: out_dir="{args.out_dir}" already exists')
    else:
        os.makedirs(args.out_dir)

    # ensure out_file not present - OK to check for the default one too, even though it'll be changed (will check again)
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    # if custom_sites is file, convert to string
    if os.path.isfile(args.custom_sites):
        args.custom_sites = open(args.custom_sites).read().rstrip()

    # die if custom_sites contains anything besides integers, whitespace, and commas
    re_NOT_d_s_comma = re.compile(r'[^\d\s\n,]')
    if re_NOT_d_s_comma.search(args.custom_sites):
        parser.error(f'\n### ERROR: custom_sites="{args.custom_sites[:20]}..." may only contain integers, whitespace, and/or commas(,)')

    # if tree is file, convert to string
    if os.path.isfile(args.tree):
        args.tree = open(args.tree).read().rstrip()

    # die if tree contains anything NON-newick; should simply fail to initialize with Tree()
    # re_NOT_newick = re.compile(r'[^\w\d.:;()\s\n,]')
    # if re_NOT_newick.search(args.tree):
    #     parser.error(f'\n### ERROR: tree="{args.tree}..." may only contain Newick-compatible characters')
    if args.tree != '':
        try:
            _ = Tree(args.tree)
        except newick.NewickError:
            parser.error(f'\n### ERROR: tree="{args.tree}" does not conform to Newick format')

    return Args(aln_file=args.aln_file,
                group_key=args.group_key,
                min_freq=args.min_freq,
                max_freq=args.max_freq,
                step_size=args.step_size,
                out_dir=args.out_dir,
                out_file=args.out_file,
                min_def_count=args.min_def_count,
                min_def_prop=args.min_def_prop,
                exclude_seq=args.exclude_seq,
                exclude_group=args.exclude_group,
                custom_sites=args.custom_sites,
                tree=args.tree)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    aln_fh_list = args.aln_file
    group_key = args.group_key
    min_freq = args.min_freq
    max_freq = args.max_freq
    step_size = args.step_size
    out_dir = args.out_dir
    out_file = args.out_file
    min_def_count = args.min_def_count
    min_def_prop = args.min_def_prop
    exclude_seq = args.exclude_seq
    exclude_group = args.exclude_group
    custom_sites = args.custom_sites
    tree = args.tree

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # INITIALIZE lists from comma-separated input
    group_key_list = group_key.split(',')
    exclude_aln_list = exclude_seq.split(',')
    exclude_group_list = exclude_group.split(',')

    # SITE list for overlap checking
    # num_custom_sites = 0
    # if custom_sites != '':
    re_s = re.compile(r'\s+')  # whitespace
    custom_sites = custom_sites.replace(',', ' ')
    # custom_sites = set(sorted(map(int, re_s.split(custom_sites))))  # actually a set
    custom_site_list = re_s.split(custom_sites)
    len_custom_site_list = len(custom_site_list)

    if len_custom_site_list == 1 and custom_site_list[0] == '':
        custom_site_list = []
        len_custom_site_list = 0
    else:
        custom_site_list = list(map(int, custom_site_list))
    custom_sites = set(custom_site_list)  # actually a set - pointless to sort
    num_custom_sites = len(custom_sites)

    if len_custom_site_list != num_custom_sites:
        sys.exit(f'\n### ERROR: Duplicates provided in --custom_sites: custom_sites={len_custom_site_list} but ' + \
                 f'only {num_custom_sites} are unique.\n')

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # group ID(s)
    regex_group_list = []
    for this_group_key in group_key_list:  # group_key.split(','):
        regex_group = re.compile(f'{this_group_key}="([^"]+)"')
        regex_group_list.append(regex_group)

    # Nucleotides
    defined_nucs = ('A', 'C', 'G', 'T', 'U')

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:aln_file="{",".join([fh.name for fh in aln_fh_list])}"')
    print(f'LOG:group_key="{group_key}"')
    print(f'LOG:min_freq="{min_freq}"')
    print(f'LOG:max_freq="{max_freq}"')
    print(f'LOG:step_size="{step_size}"')
    print(f'LOG:out_dir="{out_dir}"')
    print(f'LOG:out_file="{out_file}"')
    print(f'LOG:min_def_count="{min_def_count}"')
    print(f'LOG:min_def_prop="{min_def_prop}"')
    print(f'LOG:exclude_seq="{exclude_seq}"')
    print(f'LOG:exclude_group="{exclude_group}"')
    # print(f'LOG:custom_sites="{custom_sites}"', flush=True)
    # print(f'LOG:custom_sites="{",".join(map(str, list(custom_sites)))}"', flush=True)
    print(f'LOG:custom_sites="{",".join(map(str, sorted(list(custom_sites))))}"')
    print(f'LOG:num_custom_sites={num_custom_sites}')
    print(f'LOG:tree="{tree}"', flush=True)

    # -------------------------------------------------------------------------
    # TODO: CHECK FOR DUPLICATE SEQUENCE IDs

    # -------------------------------------------------------------------------
    # INPUT FASTA SEQUENCE(S) and COUNT ALLELES FOR ALL GROUPS

    # NOTES from old aligned_fasta_group_diffs.pl
    # Now I want to go site-by-site and (1) determine the number of defined nucleotides at
    # each site; (2) calculate the consensus/majority nucleotide among those defined using
    # some threshold like 90%; (3) in cases where majority nucleotides differ between groups,
    # make note of the sites, nucleotide identities, and coverage.

    # recs: List[SeqIO.SeqRecord] = []
    # recs = SeqIO.FastaIO.FastaIterator(source='')
    # recs = []

    group_site_nuc_count_ddd: Dict[str, Dict[int, Dict[str, int]]] = defaultdict(dict)
    group_count_dict: Dict[str, int] = defaultdict(int)
    excluded_group_count_dict: Dict[str, int] = defaultdict(int)

    grouped_nseqs = 0
    excluded_group_nseqs = 0
    no_group_nseqs = 0

    # Keep track of alignment length
    aln_length = None

    for aln_fh in aln_fh_list:
        # print(f'adding FASTA records from aln_file={aln_fh.name}')
        #recs = SeqIO.parse(aln_fh, 'fasta')
        recs = AlignIO.read(aln_fh, 'fasta')  # unit is ALIGNMENT, and .read() because we have only one ALN per file
        if aln_length is None:
            aln_length = recs.get_alignment_length()
            print(f'LOG:aln_length={aln_length}')
        elif aln_length != recs.get_alignment_length():
            sys.exit(f'\n### ERROR: file={aln_fh.name} does not have alignment length (num sites) of {aln_length}')

        nseqs = 0
        nseqs_wMeta = 0
        nseqs_woMeta = 0
        recs_woMeta: List[str] = []
        # recs = recs + SeqIO.parse(aln_fh, 'fasta')
        # print(f'type(recs)={type(recs)}')

        for rec in recs:
            if rec.id not in exclude_aln_list:  # skip if an excluded sequence
                nseqs += 1
                group_name_match = None

                for regex_group in regex_group_list:
                    group_name_match = regex_group.search(rec.description)  # match object has .group() method  # TODO no match!?

                    if group_name_match is not None:
                        break

                if group_name_match is not None:
                    nseqs_wMeta += 1
                    group_name = group_name_match.group(1)

                    if group_name not in exclude_group_list:  # skip if an excluded group
                        # print(f'group_name="{group_name}"')
                        group_count_dict[group_name] += 1
                        grouped_nseqs += 1

                        for site, nuc in enumerate(rec.seq, start=1):
                            if site not in group_site_nuc_count_ddd[group_name]:
                                group_site_nuc_count_ddd[group_name][site] = defaultdict(int)

                            # count
                            group_site_nuc_count_ddd[group_name][site][nuc] += 1
                    else:
                        excluded_group_count_dict[group_name] += 1
                        excluded_group_nseqs += 1
                else:
                    no_group_nseqs += 1
            else:
                nseqs_woMeta += 1
                recs_woMeta.append(rec.id)
                # print(f'### WARNING: no group metadata in header of {rec.id}')

        # LOG number of recs
        print(f'LOG:nseqs={nseqs} in aln_file={aln_fh.name}')
        print(f'LOG:nseqs_wMeta={nseqs_wMeta} in aln_file={aln_fh.name}')
        print(f'LOG:nseqs_woMeta={nseqs_woMeta} in aln_file={aln_fh.name}: {",".join(recs_woMeta)}')

    # pprint(f'sorted(group_site_nuc_count_ddd.keys()):')
    # pprint(sorted(group_site_nuc_count_ddd.keys()))
    # pprint(group_site_nuc_count_ddd)

    # -------------------------------------------------------------------------
    # LOG sequence counts and groups
    # Number of sequences
    nseqs = grouped_nseqs + excluded_group_nseqs + no_group_nseqs
    print(f'LOG:nseqs={nseqs}')

    if nseqs < 2:
        sys.exit(f'\n### ERROR: nseqs={nseqs} must be >1 to determine group differences')

    # if step_size is 0, use empirical frequencies
    if step_size == 0:
        step_size = 1 / nseqs

    # LOG IT
    print(f'LOG:step_size_used="{step_size}"')

    # Groups
    groups = (sorted(group_site_nuc_count_ddd.keys()))
    print(f'LOG:groups={",".join(groups)}')

    # TODO: make sure all identified groups are within the tree, otherwise warning and just use what's there

    # groups: Set[str] = ()
    # set comprehensions use curly braces
    # print('group_count_dict:')
    # pprint(group_count_dict)

    # Number of sequences GROUPED and their counts
    print(f'LOG:grouped_nseqs={grouped_nseqs}')

    group_counts_print_list = []
    for group in sorted(group_count_dict.keys()):
        # print(f'{group}={group_count_dict[group]}', end='')
        group_counts_print_list.append(f'{group}={group_count_dict[group]}')
    print(f'LOG:group_counts:{",".join(group_counts_print_list)}')

    # Number of sequences EXCLUDED and their counts
    # print('excluded_group_count_dict:')
    # pprint(excluded_group_count_dict)
    print(f'LOG:excluded_group_nseqs={excluded_group_nseqs}')
    excluded_group_count_print_list = []
    for group in sorted(excluded_group_count_dict.keys()):
        # print(f'{group}={excluded_group_count_dict[group]}', end='')
        excluded_group_count_print_list.append(f'{group}={excluded_group_count_dict[group]}')
    print(f'LOG:excluded_group_counts:{",".join(excluded_group_count_print_list)}')

    # Number of sequences with group NA
    print(f'LOG:no_group_nseqs={no_group_nseqs}', flush=True)

    # -------------------------------------------------------------------------
    # TRAVERSE the SPECIFIED RANGE of allele frequency values

    # Determine values to loop
    freqs_to_examine: List[float] = []
    freq = min_freq
    while freq <= max_freq:
        # print(f'freq={freq}')
        freqs_to_examine.append(freq)

        # increment
        freq = round(freq + step_size, ndigits=len(str(nseqs)))  # nseqs is always an int count, so this works

    # print(f'LOG:freqs_to_examine={",".join(map(str, freqs_to_examine))}')

    # INITIALIZE headers
    group_header_string = '\t'.join(['freq', 'site', 'group', 'maj_nuc', 'maj_nuc_count', 'def_nuc_count', 'nseqs',
                                     'maj_nuc_freq', 'in_custom_sites', 'num_alleles', 'alleles', 'allele_counts',
                                     # 'phyly', 'monophyletic_nodes',  # incorrect
                                     'multiallelic', 'homoplasy', 'site_tree'])
    count_header_string = '\t'.join(['freq_cutoff', 'nsites', 'custom_sites_covered', 'prop_custom_sites_covered',
                                     'nsites_multiallelic', 'nsites_homoplasy'])

    # OPEN count file and write header
    # edit out_file if default
    if out_file == 'group_defining_site_counts_from<#>_to<#>_by<#>.tsv':
        out_file = f'group_defining_site_counts_from{min_freq}_to{max_freq}_by{round(step_size, ndigits=len(str(nseqs)))}.tsv'

    # CHECK IF IT EXISTS AGAIN (ensure out_file not present; OK if repeated)
    if os.path.isfile(args.out_file):
        sys.exit(f'\n### ERROR: out_file="{args.out_file}" already exists')

    # LOG IT
    print(f'LOG:out_file_used="{out_file}"')
    count_out_file_hdl = open(os.path.join(out_file), "wt")
    count_out_file_hdl.write(f'{count_header_string}\n')

    # INITIALIZE list of points
    freq_nsites_tuples: List[tuple] = []

    # INITIALIZE list of best site candidates, i.e., those detected at the lowest possible frequency threshold
    site_candidates: List[int] = []

    # LOOP FREQUENCY CUTOFFS
    for freq in freqs_to_examine:
        # INITIALIZE OUTPUT FILE for this freq and WRITE HEADER
        # group_out_file_hdl = open(f'{out_dir}_{freq}.tsv', "wt")
        group_out_file_hdl = open(os.path.join(out_dir, f'group_defining_sites_{freq}.tsv'), "wt")
        group_out_file_hdl.write(f'{group_header_string}\n')

        group_defining_site_count = 0
        custom_sites_covered = 0
        nsites_multiallelic = 0
        nsites_homoplasy = 0

        # LOOP SITES (COLUMNS) IN ALIGNMENT
        # TODO: could probably make MUCH more efficient by storing sites first, THEN looping frequency cutoffs
        # TODO: instead of group_site_nuc_count_ddd, make a group_site_MAJNUC_FREQ_[ddd], then loop THAT
        for site in range(1, aln_length + 1):
            # for this site only
            group_majNucData: Dict[str, Dict[str, Any]] = defaultdict(dict)
            print(f'site={site}')
            maj_nuc_list = []
            # maj_nuc_set = set()  # {} will think dict
            # maj_nuc_freq_list: List[float] = []
            maj_nuc_to_max_freq: Dict[str, float] = defaultdict(float)

            # LOOP GROUPS at this site (say, lineages A/B/C/D)
            for group in groups:
                print(f"group={group}")
                nucs = group_site_nuc_count_ddd[group][site].keys()
                maj_nuc = ''
                maj_nuc_count = 0
                def_nuc_count = 0
                undef_nuc_count = 0

                for nuc in nucs:
                    nuc_count = group_site_nuc_count_ddd[group][site][nuc]

                    if nuc in defined_nucs:
                        def_nuc_count += nuc_count

                        if nuc_count > maj_nuc_count:
                            maj_nuc = nuc
                            maj_nuc_count = nuc_count
                    else:
                        undef_nuc_count += nuc_count

                total_nuc_count = def_nuc_count + undef_nuc_count
                if group_count_dict[group] != total_nuc_count:
                    sys.exit(f'\nERROR: at site={site}, group nseqs={group_count_dict[group]} does not equal total_nuc_count={total_nuc_count}')

                def_nuc_prop = def_nuc_count / total_nuc_count

                # ADD to SET of major nucleotides at this site
                if maj_nuc in defined_nucs and def_nuc_count >= min_def_count and def_nuc_prop >= min_def_prop:
                    maj_nuc_list.append(maj_nuc)
                    # maj_nuc_set.add(maj_nuc)
                    # maj_nuc_freq = NA
                    # maj_nuc_pct = NA

                    maj_nuc_freq = maj_nuc_count / def_nuc_count
                    # maj_nuc_pct = round(100 * maj_nuc_freq)
                    # maj_nuc_freq_list.append(maj_nuc_freq)

                    if maj_nuc_freq > maj_nuc_to_max_freq[maj_nuc]:
                        maj_nuc_to_max_freq[maj_nuc] = maj_nuc_freq

                # print(f'group={group};maj_nuc={maj_nuc};maj_nuc_count={maj_nuc_count};maj_nuc_freq={maj_nuc_pct}%')
                print(f'group={group};maj_nuc={maj_nuc};maj_nuc_count={maj_nuc_count};maj_nuc_freq={maj_nuc_freq}%')

                # STORE FOR THIS SITE
                group_majNucData[group]['maj_nuc'] = maj_nuc
                group_majNucData[group]['maj_nuc_count'] = maj_nuc_count
                group_majNucData[group]['def_nuc_count'] = def_nuc_count
                group_majNucData[group]['maj_nuc_freq'] = maj_nuc_freq  # used later to filter

            # EXAMINE SITE for DIFFERENT maj_nuc, with at least TWO (2) DIFFERENT nucs (1 ALT) meeting FREQ CRITERION
            if len(maj_nuc_to_max_freq.keys()) > 1 and len([x for x in maj_nuc_to_max_freq.values() if x >= freq]) >= 2:
                # pprint(group_majNucData)
                group_defining_site_count += 1

                # add to candidates
                site_candidates.append(site)  # I know this is less than efficient

                # IN CUSTOM SITES?
                in_custom_sites = 'NA'
                if num_custom_sites > 0:
                    in_custom_sites = site in custom_sites
                    if in_custom_sites:
                        custom_sites_covered += 1

                # maj_nuc_set = set(maj_nuc_list)
                # print(f'site={site};maj_nuc_set={maj_nuc_set};maj_nuc_to_max_freq={dict(maj_nuc_to_max_freq)}')

                # MULTIALLELIC
                site_num_alleles = len(set(maj_nuc_to_max_freq.keys()))
                multiallelic = site_num_alleles > 2

                if multiallelic:
                    nsites_multiallelic += 1

                site_alleles_list = sorted(list(maj_nuc_to_max_freq.keys()))
                site_alleles_set_list = sorted(list(set(site_alleles_list)))
                site_alleles = ','.join(site_alleles_set_list)

                # COUNTS for each allele (num groups with major nucleotide)
                site_allele_counts_list = [maj_nuc_list.count(i) for i in site_alleles_set_list]
                site_allele_counts = ','.join(map(str, site_allele_counts_list))

                # DETERMINE MONO/POLY/PARA-PHYLETIC
                # http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#checking-the-monophyly-of-attributes-within-a-tree
                # replace each group's name with its nucleotide
                phyly_list = []
                monophyletic_nodes_list = []
                homoplasy = 'NA'
                site_tree = 'NA'
                # print(f'site={site}')
                if tree != '':
                    print("group_majNucData:")
                    pprint(group_majNucData)

                    print(f'tree={tree}')

                    # # Copy and replace names with alleles
                    # tree_copy = str(tree)  # maybe it's a mistake to make it a string here? No
                    # # tree_copy = tree.copy()  # doesn't work for strings; and tree already IS a string, from input
                    # print(f'tree_copy BEFORE={tree_copy}')
                    # for group in groups:
                    #     print(f"this_group={group}")
                    #     tree_copy = tree_copy.replace(group, group_majNucData[group]['maj_nuc'])
                    # print(f'tree_copy AFTER={tree_copy}')
                    #
                    # # NEW APPROACH --
                    tree_copy = Tree(tree)
                    print(f'tree_copy BEFORE={tree_copy}')

                    # Iterate over each leaf and replace the name with the nucleotide allele
                    for leaf in tree_copy:
                        if leaf.name in group_majNucData:
                            leaf.name = group_majNucData[leaf.name]['maj_nuc']

                    # back to Newick format
                    tree_copy = tree_copy.write(format=0)
                    print(f'tree_copy AFTER={tree_copy}')
                    # DONE --

                    # save site's tree
                    site_tree = str(tree_copy)

                    # initialize tree object
                    tree_copy = Tree(tree_copy)
                    print("Leaves")
                    print([leaf.name for leaf in tree_copy.get_leaves()])  # DEBUG

                    print(f"site_alleles_set_list: {site_alleles_set_list}")
                    print(f"len(site_alleles_set_list): {len(site_alleles_set_list)}")

                    if len(site_alleles_set_list) > 1:
                        for nuc in site_alleles_set_list:
                            # print(f'nuc={nuc}')
                            # print(site)  # DEBUG
                            print([leaf.name for leaf in tree_copy.get_leaves() if leaf.name in [nuc]])  # DEBUG
                            tree_monophyly = tree_copy.check_monophyly(values=[nuc], target_attr='name')  # added ignore_missing=True
                            phyly = tree_monophyly[1]
                            phyly_list.append(phyly)

                            tree_monophyletic_nodes = tree_copy.get_monophyletic(values=[nuc], target_attr='name')
                            monophyletic_nodes = len(list(tree_monophyletic_nodes))
                            monophyletic_nodes_list.append(monophyletic_nodes)

                if len(phyly_list) == 0:
                    phyly_list = ['NA']

                if len(monophyletic_nodes_list) == 0:
                    monophyletic_nodes_list = ['NA']

                # CHECK for HOMOPLASY
                if len(monophyletic_nodes_list) > 0:
                    monophyletic_nodes_list_wo1s = [x for x in monophyletic_nodes_list if x != 1]  # .remove only first

                    # if multiallelic:
                    #     homoplasy = 'NA'
                    if len(monophyletic_nodes_list_wo1s) > 1:
                        homoplasy = True
                    else:
                        homoplasy = False

                if homoplasy is True:
                    nsites_homoplasy += 1

                # remove branch lengths and support values
                site_tree_simple = re.sub(r':\d+(\.\d+)?', '', site_tree)
                site_tree_simple = re.sub(r'\)\d+', '', site_tree_simple)

                # LOOP GROUPS to PRINT DATA for THIS SITE
                for group in groups:
                    # OUTPUT FILE
                    line_list = [freq, site, group,
                                 group_majNucData[group]['maj_nuc'],
                                 group_majNucData[group]['maj_nuc_count'],
                                 group_majNucData[group]['def_nuc_count'],
                                 nseqs,
                                 group_majNucData[group]['maj_nuc_freq'],
                                 in_custom_sites,
                                 site_num_alleles,
                                 site_alleles,
                                 site_allele_counts,
                                 # ','.join(phyly_list),  # incorrect
                                 # ','.join(map(str, monophyletic_nodes_list)),  # incorrect
                                 multiallelic,
                                 homoplasy,
                                 site_tree_simple]
                    # "\t".join(map(str, [freq, site, group, maj_nuc, maj_nuc_count, def_nuc_count, maj_nuc_freq, nseqs]))
                    line_string = '\t'.join(map(str, line_list))
                    # print(line_string)
                    group_out_file_hdl.write(f'{line_string}\n')

        # CLOSE group file
        # print(f'freq={freq};group_defining_site_count={group_defining_site_count}')
        group_out_file_hdl.close()

        # SAVE point
        freq_nsites_tuples.append((freq, group_defining_site_count))

        # PRINT count file
        # print(f'LOG:freq={freq};group_defining_site_count={group_defining_site_count}')
        prop_custom_sites_covered = 0
        if num_custom_sites > 0:
            prop_custom_sites_covered = custom_sites_covered / num_custom_sites
        count_out_file_hdl.write(f'{freq}\t{group_defining_site_count}\t{custom_sites_covered}\t' + \
                                f'{prop_custom_sites_covered}\t{nsites_multiallelic}\t{nsites_homoplasy}\n')

    # CLOSE count file
    count_out_file_hdl.close()

    # -------------------------------------------------------------------------
    # DETERMINE BEST FREQUENCY CUTOFF and that which minimizes normalized (to 1) distance to (xmax,ymax)
    # print(f'freq_nsites_tuples:')
    # print(freq_nsites_tuples)
    xmin = 0.5  # TODO: hardcoded theoretical limit, but consider flexibility later
    xmax = 1  # TODO: hardcoded theoretical limit, but consider flexibility later
    ymin = 0  # TODO: hardcoded theoretical limit, but consider flexibility later
    ymax = max([i for _, i in freq_nsites_tuples])
    print(f'LOG:max_nsites={ymax}')

    best_dist = 1
    best_nsites = ymax
    best_freq_cutoff = 0.5

    # LOOP normalized (x,y) values to find the one closest to (1,1)
    for xval, yval in freq_nsites_tuples:
        if (xmax - xmin) > 0:
            xnorm = (xval - xmin) / (xmax - xmin)
        else:
            xnorm = None

        if (ymax - ymin) > 0:
            ynorm = (yval - ymin) / (ymax - ymin)
        else:
            ynorm = None

        if xnorm is not None and ynorm is not None:
            this_dist = math.dist((xnorm, ynorm), (1, 1))
        else:
            this_dist = None

        if this_dist is not None and best_dist is not None:
            if this_dist < best_dist:
                best_freq_cutoff = xval
                best_nsites = yval
                best_dist = this_dist

    # LOG best frequency cutoff
    nsites_string = f'LOG:best_nsites={best_nsites}'

    # DETERMINE candidate sites
    site_candidates_str = None
    site_candidates_list = list(set(site_candidates))
    if len(site_candidates_list) > 0:
        site_candidates_str = ",".join(map(str, site_candidates_list))

    # NOTE to user if there were none
    if best_nsites == 0:
        nsites_string += ' <== NOTE: there were NO group-defining sites by any of the criteria examined'

    print(f'LOG:site_candidates={site_candidates_str}')
    print(f'LOG:best_dist={best_dist}')
    print(nsites_string)
    print(f'LOG:best_freq_cutoff={best_freq_cutoff}')

    # -------------------------------------------------------------------------
    # DONE message
    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
