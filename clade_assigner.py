#!/usr/bin/env python3
"""
Purpose: Given a set of clade representative IDs, determine the mutually exclusive members of each clade in the tree
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2022-11-10
"""

import argparse
import numpy as np
import os
# import random
import re
import sys
import time
# from Bio import AlignIO
# from Bio.Align import AlignInfo  # for .dumb_consensus()
# from Bio.Align import MultipleSeqAlignment
from collections import Counter, defaultdict
from ete3 import Tree
# from ete3.parser import newick
# from evobioinfo import GAPS, hamming, NUCS_DEFINED, NUCS_INDETERMINATE
from pprint import pprint
from typing import Dict, List, NamedTuple, Optional, TextIO

usage = """# -----------------------------------------------------------------------------
clade_assigner.py - Given representatives, determine the mutually exclusive members of each clade in the tree
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ clade_assigner.py --help
    $ pydoc ./clade_assigner.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $ clade_assigner.py --help
# -----------------------------------------------------------------------------
"""


class Args(NamedTuple):
    """ Command-line arguments """
    tree: TextIO
    reps: TextIO
    prefix: str
    permutations: int
    seed: Optional[int]  # TODO: add --root for specific OR random?
    recursion: int
    threshold: float
    maxNA: float
    begin_num: int
    normalize: bool
    suppress: bool
    verbose: bool
    ebullient: bool


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Determine the mutually exclusive members of each clade in the tree. HELP: clade_assigner.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-t',
                        '--tree',
                        metavar='FILE',
                        help='NEWICK file containing the tree(s) [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    parser.add_argument('-r',
                        '--reps',
                        metavar='FILE|str',
                        help='Sequence ID(s) (comma-separated) of clade representatives [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=str)

    # -------------------------------------------------------------------------
    # OPTIONAL
    parser.add_argument('-o',
                        '--prefix',
                        metavar='str',
                        help='Output file name PREFIX [OPTIONAL]',
                        required=False,
                        type=str,
                        default='<NEWICK>_clades')  # detect this later

    parser.add_argument('-p',
                        '--permutations',
                        metavar='int',
                        help='Number of permutations (random root and pruning orders) per tree [OPTIONAL]',
                        required=False,
                        type=int,
                        default=100)

    parser.add_argument('-s',
                        '--seed',
                        help='RANDOM NUMBER SEED to be used [OPTIONAL]',
                        required=False,
                        metavar='int',
                        type=int,
                        default=None)

    parser.add_argument('-R',
                        '--recursion',
                        help='System RECURSION LIMIT (may be exceeded for large trees); use with caution [OPTIONAL]',
                        required=False,
                        metavar='int',
                        type=int,
                        default=1000)

    parser.add_argument('-T',
                        '--threshold',
                        help='Minimum proportion of supporting permutations required to classify a sequence [OPTIONAL]',
                        required=False,
                        metavar='float',
                        type=float,
                        default=0.95)

    parser.add_argument('-m',
                        '--maxNA',
                        help='Maximum proportion of NA permutations allowed for a sequence to be assigned [OPTIONAL]',
                        required=False,
                        metavar='float',
                        type=float,
                        default=0.5)

    parser.add_argument('-b',
                        '--begin_num',
                        help='Beginning tree number, useful for generating unique tree numbers for downstream '
                             'analyses [OPTIONAL]',
                        required=False,
                        metavar='int',
                        type=int,
                        default=1)

    parser.add_argument('-n',
                        '--normalize',
                        help='Normalize evolutionary distance to largest leaf-to-rep distance [OPTIONAL]',
                        action='store_true')

    parser.add_argument('-S',
                        '--suppress',
                        help='Suppress most output [OPTIONAL]',
                        action='store_true')

    parser.add_argument('-v',
                        '--verbose',
                        help='Verbose output [OPTIONAL]',
                        action='store_true')

    parser.add_argument('-e',
                        '--ebullient',
                        help='Positively ebullient output [OPTIONAL]',
                        action='store_true')

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # FORM the output file PREFIX if not given
    if args.prefix == '<NEWICK>_clades':
        args.prefix = os.path.splitext(args.tree[0].name)[0] + '_clades_p' + str(args.permutations)

    # ensure THREE out_file are not present
    out_file_counts = args.prefix + '_counts.tsv'
    out_file_props = args.prefix + '_props.tsv'
    out_file_means = args.prefix + '_means.tsv'
    if os.path.isfile(out_file_counts) or os.path.isfile(out_file_props) or os.path.isfile(out_file_means):
        parser.error(f'\n### ERROR: out_file="{out_file_counts}"/"{out_file_props}"/"{out_file_means}" already exist')

    # if reps is file, convert to string
    # print(f'reps:')
    # print(args.reps)
    if type(args.reps) is list and os.path.isfile(args.reps[0]):
        args.reps = open(args.reps[0]).read().rstrip()

    # die if reps contains anything besides integers, words, whitespace, and commas
    # re_NOT_d_w_s_comma = re.compile(r'[^\d\w\s\n,]')
    re_NOT_FASTA_name = re.compile(r'[^\d\w\s\n,.|]')
    if re_NOT_FASTA_name.search(args.reps):
        # parser.error(f'\n### ERROR: reps="{args.reps[:20]}..." may only contain integers, whitespace, and/or commas(,)')
        parser.error(f'\n### ERROR: reps="{args.reps[:20]}..." may only contain words, numbers, decimals, whitespace, '
                     'commas(,), or pipes (|)')

    if args.seed is not None and args.seed < 1:
        parser.error(f'\n### ERROR: seed="{args.seed}" must be > 0')

    if not 1000 <= args.recursion <= 10000:
        parser.error(f'\n### ERROR: recursion="{args.recursion}" must be >= 1000 and <= 10000')

    if not 0 < args.threshold <= 1:
        parser.error(f'\n### ERROR: threshold="{args.threshold}" must be > 0 and <= 1')

    if not args.begin_num >= 1:
        parser.error(f'\n### ERROR: begin_num="{args.begin_num}" must be >= 1')

    return Args(tree=args.tree[0],
                reps=args.reps,
                prefix=args.prefix,
                permutations=args.permutations,
                seed=args.seed,
                recursion=args.recursion,
                threshold=args.threshold,
                maxNA=args.maxNA,
                begin_num=args.begin_num,
                normalize=args.normalize,
                suppress=args.suppress,
                verbose=args.verbose,
                ebullient=args.ebullient)


# -----------------------------------------------------------------------------
def main() -> None:
    """ The one day the sun shines in Taipei I'm coding """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    tree_fh = args.tree
    reps = args.reps
    prefix = args.prefix
    permutations = args.permutations
    seed = args.seed
    recursion = args.recursion
    threshold = args.threshold
    maxNA = args.maxNA
    begin_num = args.begin_num
    normalize = args.normalize
    suppress = args.suppress
    verbose = args.verbose
    ebullient = args.ebullient

    # THREE OUTPUT FILES from prefix
    out_file_counts = prefix + '_counts.tsv'
    out_file_props = prefix + '_props.tsv'
    out_file_means = prefix + '_means.tsv'

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # PASS

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:tree="{tree_fh.name}"')
    print(f'LOG:reps="{reps}"')
    print(f'LOG:out_file_counts="{out_file_counts}"')
    print(f'LOG:out_file_counts="{out_file_props}"')
    print(f'LOG:out_file_counts="{out_file_means}"')
    print(f'LOG:permutations={permutations}')
    print(f'LOG:seed={seed}')
    print(f'LOG:recursion={recursion}')
    print(f'LOG:threshold={threshold}')
    print(f'LOG:maxNA={maxNA}')
    print(f'LOG:begin_num={begin_num}')
    print(f'LOG:normalize={normalize}')
    print(f'LOG:suppress={suppress}')
    print(f'LOG:verbose={verbose}')
    print(f'LOG:ebullient={ebullient}')

    # RECURSION LIMIT
    sys.setrecursionlimit(recursion)

    # -------------------------------------------------------------------------
    # INITIALIZE lists from comma-separated input
    re_s = re.compile(r'\s+')  # whitespace

    # replace commas with whitespace
    reps = reps.replace(',', ' ')
    rep_seq_names = tuple(re_s.split(reps))
    rep_seq_names_count = len(rep_seq_names)
    rep_seq_names_uniq_count = len(set(rep_seq_names))

    # LOG RESULT
    print(f'RESULT:rep_seq_names_count={rep_seq_names_count}')

    if verbose:
        print('rep_seq_names:')
        print(rep_seq_names)

    if rep_seq_names_count != rep_seq_names_uniq_count:
        sys.exit(f'\n### ERROR: Duplicates provided in --reps: count={rep_seq_names_count} but ' + \
                 f'only {rep_seq_names_uniq_count} are unique.\n')

    # -------------------------------------------------------------------------
    # Calculate or record random number seed for random permutations of reps
    if seed is None:  # emulates Slatkin Item 24
        seed = np.random.randint(2 ** 32 - 1)

    # Set random number seed
    np.random.seed(seed)  # np.random does not affect random module

    print(f'LOG:seed={seed}')
    print(f'RESULT:seed_rand_result={np.random.rand()}')

    # -------------------------------------------------------------------------
    # LOOP TREES (one per line)

    # PREPARE a dict[seqID] -> list
    ID_rep_count_ddl: Dict[str, Dict[str, list]] = defaultdict(dict)  # ID => REP_ID => LIST OF COUNTS
    ID_prop_ddl: Dict[str, Dict[str, list]] = defaultdict(dict)  # ID => REP_ID => LIST OF PROPS
    ID_rep_dist_ddl: Dict[str, Dict[str, list]] = defaultdict(dict)  # ID => REP_ID => LIST OF PROPS

    ID_count_total_dl: Dict[str, list] = defaultdict(list)  #
    ID_count_exp_dl: Dict[str, list] = defaultdict(list)
    ID_count_assn_dl: Dict[str, list] = defaultdict(list)  # ID => LIST OF COUNT SUMS

    rep_seq_names_list_sorted: List[str] = []
    rep_seq_names_list_sorted_wNA : List[str] = []

    # initialize tree_num
    tree_num = begin_num - 1  # tree_num = 0
    for tree in tree_fh:
        tree = tree.rstrip()
        tree_num += 1

        # # Debugging break
        # if tree_num > 2:
        #     break

        print('\n# -------------------------------------------------------------------------')
        print(f'Tree {tree_num}')
        if not suppress:
            print('permutations', flush=True)

        # ---------------------------------------------------------------------
        # INPUT TREE using ete3
        tree_node = Tree(tree)

        # # Determine topology height of each node
        # rep_height_topology = defaultdict(int)
        # for this_rep_name in rep_seq_names:
        #     node = tree_node.search_nodes(name=this_rep_name)[0]
        #     rep_height_topology[this_rep_name] = tree_node.get_distance(node, topology_only=True)
        #
        # print('\nrep_height_topology:')
        # pprint(rep_height_topology)

        # ---------------------------------------------------------------------
        # STORE DISTANCES between each sequence and each rep, ignoring self comparisons
        all_leaf_names = tree_node.get_leaf_names()
        leaf_to_rep_dists: Dict[str, Dict[str, float]] = defaultdict(dict)
        leaf_to_closest_rep_name: Dict[str, str] = defaultdict(str)
        leaf_to_closest_rep_dist: Dict[str, float] = defaultdict(float)
        largest_dist = 0

        for leaf_name in all_leaf_names:
            # initialize
            leaf_to_rep_dists[leaf_name]: Dict[str, float] = defaultdict(float)
            closest_rep_name = ''
            closest_rep_dist = 100000.  # something arbitrary, could do better

            for rep_name in rep_seq_names:
                if leaf_name != rep_name:
                    rep_node = tree_node&rep_name
                    leaf_node = tree_node&leaf_name
                    this_dist = leaf_node.get_distance(rep_node)

                    # store
                    leaf_to_rep_dists[leaf_name][rep_name] = this_dist

                    # check if closest
                    if this_dist < closest_rep_dist:
                        closest_rep_name = str(rep_name)
                        closest_rep_dist = float(this_dist)

                    # check if largest
                    if this_dist > largest_dist:
                        largest_dist = this_dist

                else:
                    closest_rep_name = 'REP'
                    closest_rep_dist = 0.

            # store closest rep for this leaf
            leaf_to_closest_rep_name[leaf_name] = closest_rep_name
            leaf_to_closest_rep_dist[leaf_name] = closest_rep_dist

        # NORMALIZE to LARGEST DISTANCE
        if normalize:
            for leaf_name in all_leaf_names:
                for rep_name in rep_seq_names:
                    leaf_to_rep_dists[leaf_name][rep_name] = leaf_to_rep_dists[leaf_name][rep_name] / largest_dist
                    leaf_to_closest_rep_dist[leaf_name] = leaf_to_closest_rep_dist[leaf_name] / largest_dist

        # print('leaf_to_rep_dists:')
        # pprint(dict(leaf_to_rep_dists))
        # print('leaf_to_closest_rep_name:')
        # pprint(dict(leaf_to_closest_rep_name))
        # print('leaf_to_closest_rep_dist:')
        # pprint(dict(leaf_to_closest_rep_dist))

        # NOTES ---
        # these dictionaries have VALUES that correspond to this ONE TREE
        # this is the SAME as the counts/props for the permutations, which also corresponds to THIS ONE TREE
        # JUST LIKE THE MEANS of the trees -- same deal

        # ---------------------------------------------------------------------
        # STORE number of steps required to maximize clade members without including another representative
        rep_steps_dl = defaultdict(list)
        # rep_height_dl = defaultdict(list)
        # rep_height_topology_dl = defaultdict(list)

        # STORE assignments of each ID
        seq_rep_count_dd = defaultdict(dict)

        # initialize each dict of each ID
        for leaf_name in all_leaf_names:
            seq_rep_count_dd[leaf_name] = defaultdict(int)

        # print(f'LENGTH of all_leaf_names={len(all_leaf_names)}')  # LENGTH of all_leaf_names=6052 Q.E.D.
        # initiate list of rep seq names that can be shuffled
        rep_seq_names_list = list(rep_seq_names)
        rep_seq_names_list_sorted = sorted(rep_seq_names_list)
        # rep_seq_names_list_sorted = sorted(list(rep_seq_names))

        # version with NA
        rep_seq_names_list_sorted_wNA = sorted(list(rep_seq_names_list_sorted))
        rep_seq_names_list_sorted_wNA.append('NA')

        # -------------------------------------------------------------------------
        # INITIALIZE the dl and ddl if FIRST TREE
        if tree_num == begin_num:  # if tree_num == 1:
            # ID_rep_count_ddl: ID => REP_ID => LIST OF COUNTS
            # ID_rep_count_ddl: Dict[str, Dict[str, list]] = defaultdict(dict)

            for leaf_name in all_leaf_names:
                ID_count_total_dl[leaf_name] = []
                ID_count_exp_dl[leaf_name] = []
                ID_count_assn_dl[leaf_name] = []

                for rep_name in rep_seq_names_list_sorted_wNA:  # rep_seq_names_list_sorted:
                    ID_rep_count_ddl[leaf_name][rep_name] = []
                    ID_prop_ddl[leaf_name][rep_name] = []
                    ID_rep_dist_ddl[leaf_name][rep_name] = []

                # # NA
                # ID_rep_count_ddl[leaf_name]['NA'] = []
                # ID_prop_ddl[leaf_name]['NA'] = []

        # -------------------------------------------------------------------------
        # LOOP permutations
        for i in range(permutations):
            # Make a COPY of the tree for detaching
            tree_node_copy = tree_node.copy()

            # ROOT the TREE COPY using a RANDOMLY SELECTED representative node
            random_rep_name = np.random.choice(rep_seq_names_list)  # str() didn't work
            # print(f'\nrandom_rep_name="{random_rep_name}" type {type(random_rep_name)}')

            # PRINT PROGRESS
            if not suppress:
                print('{:<7}'.format(i+1) + f'root {random_rep_name} - prune ', end='')

            random_rep_node = tree_node_copy.search_nodes(name=random_rep_name)[0]
            tree_node_copy.set_outgroup(random_rep_node)

            # DEBUG
            # print(f'rep_seq_names_list="{rep_seq_names_list}"')
            # Make shuffled list of reps with random root subtracted
            # rep_seq_names_shuffled = list(set(rep_seq_names_list).difference(set(random_rep_name)))
            rep_seq_names_shuffled = list(rep_seq_names_list)
            # print(f'rep_seq_names_shuffled1="{rep_seq_names_shuffled}" type="{type(rep_seq_names_shuffled)}"')
            # rep_seq_names_shuffled = [str(item) for item in rep_seq_names_shuffled]
            # print(f'rep_seq_names_shuffled2="{rep_seq_names_shuffled}"')
            rep_seq_names_shuffled.remove(random_rep_name)
            # print(f'rep_seq_names_shuffled3="{rep_seq_names_shuffled}"')
            np.random.shuffle(rep_seq_names_shuffled)
            # print(f'rep_seq_names_shuffled4="{rep_seq_names_shuffled}"')

            # Loop all reps and store their nodes
            # if ebullient:
            #     rooted_tree_name = str(out_file)  # TODO
            #     rooted_tree_name.rstrip(".tsv")
            #     rooted_tree_name = rooted_tree_name + '_root' + random_rep_name + ".tree"
            #     tree_node_copy.write(format=1, outfile=rooted_tree_name)
            leaf_names_assigned = []
            random_rep_order = []

            # for count, this_rep_name in enumerate(rep_seq_names_list, start=1):  # SHUFFLED; count is 1, 2, ..., length
            for count, this_rep_name in enumerate(rep_seq_names_shuffled, start=1):  # doesn't include the root
                if count < len(rep_seq_names_shuffled):
                    # if this_rep_name != random_rep_name and count < len(rep_seq_names_list) - 1:
                    # if this_rep_name == random_rep_name:
                    #     count -= 1  # don't count the matching one
                    # elif count < len(rep_seq_names_list) - 1:

                    node = tree_node_copy.search_nodes(name=this_rep_name)[0]
                    random_rep_order.append(this_rep_name)
                    # # get distance from root node
                    # rep_height_dl[this_rep_name].append(tree_node_copy.get_distance(node))
                    #
                    # # get topology distance from root node
                    # rep_height_topology_dl[this_rep_name].append(tree_node_copy.get_distance(node, topology_only=True))

                    # COUNTER for number of nodes visited in search for this rep's clade
                    steps = 0

                    while node:
                        leaf_names = node.get_leaf_names()

                        if verbose:
                            print(f'steps={steps} (n={len(leaf_names)} leaves)')

                        if ebullient:
                            print(f'steps={steps} Leaves: {leaf_names} (n={len(leaf_names)})')

                        # Check whether leaf names contain any other representatives
                        leaves_in_reps = [leaf for leaf in leaf_names if leaf in rep_seq_names]

                        # if verbose:
                        #     print(f'leaves_in_reps: {leaves_in_reps} (n={len(leaves_in_reps)})')
                        #     print("\n================================================================================")

                        # WENT TOO FAR; die when matching other ref
                        if len(leaves_in_reps) > 1:
                            if verbose:
                                print(f'At steps={steps}, included {len(leaves_in_reps)} representatives; store steps={steps - 1}')
                                print(f'Found clade containing the {this_rep_name} representative, {steps} nodes from leaf')
                            # rep_steps[this_rep_node_name] = steps
                            rep_steps_dl[this_rep_name].append(steps - 1)  # because we went too far

                            if verbose:
                                print("\n================================================================================")
                                print("================================================================================")
                            break
                        else:
                            # Move up the tree
                            node = node.up
                            steps += 1

                    # FIND CLADE NODE & DETACH
                    node = tree_node_copy.search_nodes(name=this_rep_name)[0]
                    step_counter = 0
                    for i in range(steps - 1):  # we stored the the right number of steps; for 2 steps, want len([0,1]] steps
                        step_counter += 1
                        node = node.up

                    if verbose:
                        print(f'Traversed step_counter={step_counter} steps')

                    # GET LEAVES and SAVE DATA
                    leaf_names = node.get_leaf_names()
                    if verbose:
                        print(f'leaf_names (n={len(leaf_names)}):')
                        print(f'n={len(leaf_names)} leaves assigned to clade={this_rep_name}')
                    # leaf_names_counter_dict = dict(Counter(leaf_names))

                    if ebullient:
                        print(leaf_names)

                    for leaf_name in leaf_names:
                        if verbose:
                            print(f'name={leaf_name}, this_rep_name={this_rep_name}')
                        seq_rep_count_dd[leaf_name][this_rep_name] += 1

                    leaf_names_assigned.extend(leaf_names)

                    # DETACH
                    node.detach()

            # TODO: REMAINING leaves are classified as NA
            leaf_names_NA = sorted(list(set(all_leaf_names).difference(set(leaf_names_assigned))))
            # print(f'len(leaf_names_NA)={len(leaf_names_NA)}')
            for leaf_name in leaf_names_NA:
                seq_rep_count_dd[leaf_name]['NA'] += 1

            # Print order of pruning
            if not suppress:
                print(','.join(random_rep_order))  # + f' (n={len(random_rep_order)})')

            # print("\n================================================================================")
            # print('RESULTS')
            #
            # print('\nrep_steps:')
            # pprint(rep_steps)
            #
            # print('\nrep_height:')
            # pprint(rep_height)
            #
            # print('\nrep_height_topology:')
            # pprint(rep_height_topology)
            #
            # # Sort by topological height
            # rep_height_topology_sorted_ids = [id for dist, id in sorted(zip(rep_height_topology.values(), rep_height_topology.keys()), reverse=True)]
            # print('\nrep_height_topology_sorted_ids:')
            # pprint(rep_height_topology_sorted_ids)

        # print('\nmax_height_node:')

        print("\n================================================================================")
        print('RESULTS')

        if verbose:
            print('\nrep_steps_dl:')
            print(rep_steps_dl)

        # COUNT times each number of steps observed
        if not suppress:
            print('\nrep_steps_dl - max steps toward root before another clade representative included (STEPS: PERMUTS):')
            for name in rep_seq_names:
                step_counts = Counter(rep_steps_dl[name])
                print(f'{name}:')
                pprint(step_counts)

            if verbose:
                # print('\nrep_height_dl:')
                # print(rep_height_dl)

                print('\nseq_rep_count_dd:')
                pprint(seq_rep_count_dd)

                # print('\nrep_height_topology_dl:')
                # print(rep_height_topology_dl)

        # -------------------------------------------------------------------------
        # Calculate the PROPORTIONAL assignment of each sequence to each clade

        # verify correct number of sequences
        print(f'RESULT:num_seqs_logged={len(seq_rep_count_dd)}')

        # PRINT proportions to output file
        # initialize output file
        # out_file_hdl = open(out_file, "wt")
        out_file_tree_count = f'{prefix}_tree{tree_num}_counts.tsv'
        out_file_tree_count_hdl = open(out_file_tree_count, "wt")

        out_file_tree_prop = f'{prefix}_tree{tree_num}_props.tsv'
        out_file_tree_prop_hdl = open(out_file_tree_prop, "wt")

        out_file_tree_dist = f'{prefix}_tree{tree_num}_dists.tsv'
        out_file_tree_dist_hdl = open(out_file_tree_dist, "wt")

        # FORM header for counts and props
        header_list = ['seq_name', 'permut_total', 'permut_assn_exp', 'permut_assn_obs', 'permut_assn_NA']
        header_list.extend(rep_seq_names_list_sorted_wNA)

        # FORM heade for dists
        header_dists_list = ['seq_name', 'permut_total', 'permut_assn_exp', 'permut_assn_obs', 'permut_assn_NA',
                             'closest_rep', 'closest_rep_dist']

        # additional columns for DISTANCES to each rep
        rep_seq_names_list_sorted_dist = []
        for this_rep in rep_seq_names_list_sorted:
            rep_seq_names_list_sorted_dist.append(this_rep + '_dist')
            # header_dists_list.append(this_rep + '_dist')

        header_dists_list.extend(rep_seq_names_list_sorted_dist)

        # list all the reps twice?
        # leaf_to_rep_dists
        # leaf_to_closest_rep_name
        # leaf_to_closest_rep_dist

        # header_list.append('NA')
        header = '\t'.join(header_list)
        header_dists = '\t'.join(header_dists_list)

        # PRINT header
        # print(header)
        out_file_tree_count_hdl.write(header + '\n')
        out_file_tree_prop_hdl.write(header + '\n')
        out_file_tree_dist_hdl.write(header_dists + '\n')

        # TODO: I guess one solution is simply to make LISTS at THIS step instead of single values
        # Then a subsequent step would compute the means (or whatever) of those lists
        # This would also allow a printing of two versions: one of the, say, comma-separated lists, the other their means

        # ID_count_total_dl[leaf_name] = []
        # ID_count_exp_dl[leaf_name] = []
        # ID_count_assn_dl[leaf_name] = []

        # LOOP ALL seqs, print proportional assignment to each representative (column)
        for name in sorted(seq_rep_count_dd.keys()):
            # if name in rep_seq_names:  # was thinking to treat the representatives themselves specially, but naw
            # count_list: List = []
            count_dict: Dict[str, int] = defaultdict(int)

            # LOOP REP seqs (output columns), append number of times seq assigned
            for rep_name in rep_seq_names_list_sorted_wNA:  # rep_seq_names_list_sorted:
                # line_list.append(seq_rep_count_dd[name][rep_name])
                # count_list.append(seq_rep_count_dd[name][rep_name])  # NUMBER of times this SEQ assigned to this REP
                count_dict[rep_name] = seq_rep_count_dd[name][rep_name]
                ID_rep_count_ddl[name][rep_name].append(seq_rep_count_dd[name][rep_name])

            # # SPECIALLY ADD NA
            # count_list.append(seq_rep_count_dd[name]['NA'])  # NUMBER of times this SEQ assigned to this REP
            # ID_rep_count_ddl[name]['NA'].append(seq_rep_count_dd[name]['NA'])

            exp_permut = permutations * (len(rep_seq_names_list) - 2) / len(rep_seq_names_list)  # NOT with NA
            # line_list = [name, exp_permut]

            # convert counts to proportions and append to line_list
            # count_sum = sum(count_list)
            count_total_sum = sum(count_dict.values())
            count_assn_sum = sum(count_dict.values()) - count_dict['NA']

            # INITIATE line lists
            tree_count_line_list = [name, count_total_sum, exp_permut, count_assn_sum]
            tree_prop_line_list = [name, count_total_sum, exp_permut, count_assn_sum]
            tree_dist_line_list = [name, count_total_sum, exp_permut, count_assn_sum]

            ID_count_total_dl[name].append(count_total_sum)
            ID_count_exp_dl[name].append(exp_permut)
            ID_count_assn_dl[name].append(count_assn_sum)

            # add NA first for ALL
            if count_total_sum > 0:  # always true
                this_prop = count_dict['NA'] / count_total_sum
                tree_count_line_list.append(count_dict['NA'])
                ID_prop_ddl[name]['NA'].append(this_prop)
                tree_prop_line_list.append(this_prop)
                tree_dist_line_list.append(this_prop)
            else:
                tree_count_line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
                ID_prop_ddl[name]['NA'].append('NA')
                tree_prop_line_list.append('NA')
                tree_dist_line_list.append('NA')

            # TODO: add the closest rep and its distance
            tree_dist_line_list.append(leaf_to_closest_rep_name[name])
            tree_dist_line_list.append(leaf_to_closest_rep_dist[name])

            # add all reps (not NA)
            for rep_name in rep_seq_names_list_sorted:  # no NA
                count = count_dict[rep_name]
                if count_assn_sum > 0:
                    this_prop = count / count_assn_sum
                    tree_count_line_list.append(count)
                    tree_prop_line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
                    ID_prop_ddl[name][rep_name].append(this_prop)
                else:
                    tree_count_line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
                    ID_prop_ddl[name][rep_name].append('NA')
                    tree_prop_line_list.append('NA')

                # TODO REGARDLESS, the distance
                # leaf_to_rep_dists
                # leaf_to_closest_rep_name
                # leaf_to_closest_rep_dist

                # append dist to line
                tree_dist_line_list.append(leaf_to_rep_dists[name][rep_name])

                # append dist to master list
                ID_rep_dist_ddl[name][rep_name].append(leaf_to_rep_dists[name][rep_name])

            # Add NA LAST for COUNTS ONLY
            if count_total_sum > 0:
                this_prop = count_dict['NA'] / count_total_sum
                tree_count_line_list.append(count_dict['NA'])
                tree_prop_line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
                # ID_prop_ddl[name]['NA'].append(this_prop)
            else:
                tree_count_line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
                tree_prop_line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
                # ID_prop_ddl[name]['NA'].append('NA')

            # for i, count in enumerate(count_list):
            #     if count_sum > 0:
            #         this_prop = count / count_sum
            #         line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
            #         ID_prop_ddl[name][rep_seq_names_list_sorted_wNA[i]].append(this_prop)
            #     else:
            #         line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
            #         ID_prop_ddl[name][rep_seq_names_list_sorted_wNA[i]].append('NA')

            # ASSEMBLE LINE
            tree_count_line = '\t'.join(map(str, tree_count_line_list))
            tree_prop_line = '\t'.join(map(str, tree_prop_line_list))
            tree_dist_line = '\t'.join(map(str, tree_dist_line_list))
            # print(line)

            # WRITE
            out_file_tree_count_hdl.write(tree_count_line + '\n')
            out_file_tree_prop_hdl.write(tree_prop_line + '\n')
            out_file_tree_dist_hdl.write(tree_dist_line + '\n')

        # CLOSE
        out_file_tree_count_hdl.close()
        out_file_tree_prop_hdl.close()
        out_file_tree_dist_hdl.close()

        # ---------------------------------------------------------------------
        # TODO check for convergence WITHIN tree

    # FINISH LAST TREE

    # -------------------------------------------------------------------------
    # MEAN RESULTS OVER ALL TRESS

    # PRINT results to THREE FILES
    # out_file_counts = prefix + '_counts.tsv'
    # out_file_props = prefix + '_props.tsv'
    # out_file_means = prefix + '_means.tsv'
    out_file_counts_hdl = open(out_file_counts, "wt")
    out_file_props_hdl = open(out_file_props, "wt")
    out_file_means_hdl = open(out_file_means, "wt")

    # header
    header_start_list = ['seq_name', 'permut_total', 'permut_assn_exp', 'permut_assn_obs', 'permut_assn_NA']
    # header_start_list.extend(rep_seq_names_list_sorted_wNA)
    header_start = '\t'.join(header_start_list)
    # print(header)

    out_file_counts_hdl.write(header_start + '\t' + '\t'.join(rep_seq_names_list_sorted) + '\n')
    out_file_props_hdl.write(header_start + '\t' + '\t'.join(rep_seq_names_list_sorted) + '\n')
    out_file_means_hdl.write(header_start + '\t' +
                             'closest_rep\tclosest_rep_dist\t' +
                             f'top_rep\ttop_rep_prop\tassignment_{threshold}\t' +
                             '\t'.join(rep_seq_names_list_sorted) + '\t' +
                             '\t'.join(rep_seq_names_list_sorted_dist) + '\n')

    ##
    # line_list = [name, count_total_sum, exp_permut, count_assn_sum, count_dict['NA']]
    #
    # ID_count_total_dl[name].append(count_total_sum)
    # ID_count_exp_dl[name].append(exp_permut)
    # ID_count_assn_dl[name].append(count_assn_sum)
    #
    # # add NA first
    # if count_total_sum > 0:  # always true
    #     this_prop = count_dict['NA'] / count_total_sum
    #     line_list.append(count_dict['NA'])
    #     ID_prop_ddl[name]['NA'].append(this_prop)
    # else:
    #     line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
    #     ID_prop_ddl[name]['NA'].append('NA')
    #
    # # add all reps (not NA)
    # for rep_name, count in count_dict.items():
    #     if rep_name != 'NA':
    #         if count_assn_sum > 0:
    #             this_prop = count / count_assn_sum
    #             line_list.append(count)  # PROPORTION of times this SEQ assigned to this REP
    #             # line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
    #             ID_prop_ddl[name][rep_name].append(this_prop)
    #         else:
    #             line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
    #             ID_prop_ddl[name][rep_name].append('NA')
    #
    # for name in sorted(seq_rep_count_dd.keys()):
    #     # if name in rep_seq_names:  # was thinking to treat the representatives themselves specially, but naw
    #     # count_list: List = []
    #     count_dict: Dict[str, int] = defaultdict(int)
    #
    #     # LOOP REP seqs (output columns), append number of times seq assigned
    #     for rep_name in rep_seq_names_list_sorted_wNA:  # rep_seq_names_list_sorted:
    #         # line_list.append(seq_rep_count_dd[name][rep_name])
    #         # count_list.append(seq_rep_count_dd[name][rep_name])  # NUMBER of times this SEQ assigned to this REP
    #         count_dict[rep_name] = seq_rep_count_dd[name][rep_name]  # TODO ended here
    #         ID_rep_count_ddl[name][rep_name].append(seq_rep_count_dd[name][rep_name])
    #
    #     # # SPECIALLY ADD NA
    #     # count_list.append(seq_rep_count_dd[name]['NA'])  # NUMBER of times this SEQ assigned to this REP
    #     # ID_rep_count_ddl[name]['NA'].append(seq_rep_count_dd[name]['NA'])
    #
    #     exp_permut = permutations * (len(rep_seq_names_list) - 2) / len(rep_seq_names_list)  # NOT with NA
    #     line_list = [name, exp_permut]
    #     ID_exp_permut_dl[name].append(exp_permut)
    #
    #     # convert counts to proportions and append to line_list
    #     # count_sum = sum(count_list)
    #     count_sum = sum(count_dict.values()) - count_dict['NA']
    #     line_list.append(count_sum)
    #     ID_count_sum_dl[name].append(count_sum)
    #
    #     for rep_name, count in count_dict.items():
    #         if count_sum > 0:
    #             this_prop = count / count_sum
    #             line_list.append(count)  # PROPORTION of times this SEQ assigned to this REP
    #             # line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
    #             ID_prop_ddl[name][rep_name].append(this_prop)
    #         else:
    #             line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
    #             ID_prop_ddl[name][rep_name].append('NA')
    #
    #     # for i, count in enumerate(count_list):
    #     #     if count_sum > 0:
    #     #         this_prop = count / count_sum
    #     #         line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
    #     #         ID_prop_ddl[name][rep_seq_names_list_sorted_wNA[i]].append(this_prop)
    #     #     else:
    #     #         line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
    #     #         ID_prop_ddl[name][rep_seq_names_list_sorted_wNA[i]].append('NA')
    #
    #     line = '\t'.join(map(str, line_list))
    #     # print(line)
    #     out_file_hdl.write(line + '\n')
    #
    #     # for name in sorted(seq_rep_count_dd.keys()):
    #     #     # if name in rep_seq_names:  # was thinking to treat the representatives themselves specially, but naw
    #     #     # count_list: List = []
    #     #     count_dict: Dict[str, int] = defaultdict(int)
    #     #
    #     #     # LOOP REP seqs (output columns), append number of times seq assigned
    #     #     for rep_name in rep_seq_names_list_sorted_wNA:  # rep_seq_names_list_sorted:
    #     #         # line_list.append(seq_rep_count_dd[name][rep_name])
    #     #         # count_list.append(seq_rep_count_dd[name][rep_name])  # NUMBER of times this SEQ assigned to this REP
    #     #         count_dict[rep_name] = seq_rep_count_dd[name][rep_name]  # TODO ended here
    #     #         ID_rep_count_ddl[name][rep_name].append(seq_rep_count_dd[name][rep_name])
    #     #
    #     #     # # SPECIALLY ADD NA
    #     #     # count_list.append(seq_rep_count_dd[name]['NA'])  # NUMBER of times this SEQ assigned to this REP
    #     #     # ID_rep_count_ddl[name]['NA'].append(seq_rep_count_dd[name]['NA'])
    #     #
    #     #     exp_permut = permutations * (len(rep_seq_names_list) - 2) / len(rep_seq_names_list)  # NOT with NA
    #     #     line_list = [name, exp_permut]
    #     #     ID_exp_permut_dl[name].append(exp_permut)
    #     #
    #
    #     #
    #     #     # for i, count in enumerate(count_list):
    #     #     #     if count_sum > 0:
    #     #     #         this_prop = count / count_sum
    #     #     #         line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
    #     #     #         ID_prop_ddl[name][rep_seq_names_list_sorted_wNA[i]].append(this_prop)
    #     #     #     else:
    #     #     #         line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
    #     #     #         ID_prop_ddl[name][rep_seq_names_list_sorted_wNA[i]].append('NA')
    #     #
    #     #     line = '\t'.join(map(str, line_list))
    #     #     # print(line)
    #     #     out_file_hdl.write(line + '\n')
    #     #
    #     #
    #     #

    # OUTPUT:
    # ID_rep_count_ddl: Dict[str, Dict[str, list]] = defaultdict(dict)  # ID => REP_ID => LIST OF COUNTS
    # ID_prop_ddl: Dict[str, Dict[str, list]] = defaultdict(dict)  # ID => REP_ID => LIST OF PROPS
    # ID_count_total_dl[leaf_name] = []
    # ID_count_exp_dl[leaf_name] = []
    # ID_count_assn_dl[leaf_name] = []  # ID => LIST OF COUNT SUMS

    for name in sorted(ID_rep_count_ddl.keys()):
        count_list = []  # each item is a comma-separated string
        prop_list = []
        mean_list = []
        dist_list = []

        # count_dict: Dict[str, str] = defaultdict(str)
        # prop_dict: Dict[str, str] = defaultdict(str)
        # mean_dict: Dict[str, float] = defaultdict(float)

        # NA first
        count_list.append(','.join(map(str, ID_rep_count_ddl[name]['NA'])))  # NUM times SEQ assigned to REP
        prop_list.append(','.join(map(str, ID_prop_ddl[name]['NA'])))

        this_NA_mean = np.mean(ID_prop_ddl[name]['NA'])  # this is of TOTAL, not of ASSIGNED
        # mean_list.append(this_NA_mean)  # put it in another place later

        # Track top assn and value
        top_rep_choice = 'NA'
        top_rep_prop = 0.
        rep_assigned = 'NA'

        # Track closest rep and dist
        closest_rep = 'NA'
        closest_rep_d = 100000.  # something arbitrary, could do better

        # TODO: here, get the mean distance for each rep_name and save in a dict
        # then, the mean dist to each rep can be calculated and the lowest mean determined

        # INITIALIZE mean distances; this refers to ONLY THIS SEQ NAME
        # rep_dist_dl: Dict[str, list] = defaultdict(list)  # REP_ID => LIST OF DISTS

        for rep_name in rep_seq_names_list_sorted:  # no NA
            # -----------------------------------------------------------------
            # COUNTS and PROPS

            # count_dict[rep_name] = ','.join(map(str, ID_rep_count_ddl[name][rep_name]))  # NUM times SEQ assigned to REP
            # prop_dict[rep_name] = ','.join(map(str, ID_prop_ddl[name][rep_name]))
            # mean_dict[rep_name] = np.mean(ID_prop_ddl[name][rep_name])
            count_list.append(','.join(map(str, ID_rep_count_ddl[name][rep_name])))  # NUM times SEQ assigned to REP
            prop_list.append(','.join(map(str, ID_prop_ddl[name][rep_name])))

            # this_rep_mean = np.mean(ID_prop_ddl[name][rep_name])
            # this_rep_mean = np.mean(map(float, ID_prop_ddl[name][rep_name]))
            this_rep_numer = 0
            this_rep_denom = 0

            for this_num in ID_prop_ddl[name][rep_name]:
                # print(f'name={name} | rep_name={rep_name} | this_num={this_num}')
                if this_num == 'NA':
                    this_rep_numer = 'NA'
                    this_rep_denom = 'NA'
                else:
                    this_rep_numer += this_num
                    this_rep_denom += 1

            this_rep_mean = 'NA'

            if this_rep_denom != 'NA' and this_rep_denom > 0:
                this_rep_mean = this_rep_numer / this_rep_denom

            mean_list.append(this_rep_mean)

            # update top choice if beats previous
            if this_rep_mean == 'NA':
                top_rep_prop = 'NA'
            elif this_rep_mean > top_rep_prop:
                top_rep_prop = this_rep_mean
                top_rep_choice = rep_name

            # -----------------------------------------------------------------
            # DISTANCES
            # rep_dist_dl[rep_name].append(leaf_to_rep_dists[name][rep_name])
            # print('ID_rep_dist_ddl[name][rep_name]')
            # pprint(ID_rep_dist_ddl[name][rep_name])

            this_dist_numer = 0
            this_dist_denom = 0

            # print(f'ID_rep_dist_ddl[name][rep_name] where name={name} rep_name={rep_name}')
            # pprint(ID_rep_dist_ddl[name][rep_name])

            for this_dist in ID_rep_dist_ddl[name][rep_name]:
                if this_dist == 'NA':
                    this_dist_numer = 'NA'
                    this_dist_denom = 'NA'
                else:
                    this_dist_numer += this_dist
                    this_dist_denom += 1

            this_dist_mean = 'NA'

            if this_dist_denom != 'NA' and this_dist_denom > 0:
                this_dist_mean = this_dist_numer / this_dist_denom

            dist_list.append(this_dist_mean)  # TODO: is this used for anything? need dist file?

            # update top choice if beats previous
            if this_dist_mean == 'NA':
                closest_rep_d = 'NA'
            elif this_dist_mean < closest_rep_d:
                closest_rep = rep_name
                closest_rep_d = this_dist_mean

        # TODO: comeback and devise a different strategy for assignment?

        # ASSIGN clade if threshold met
        if this_NA_mean < maxNA and top_rep_prop >= threshold:
            rep_assigned = top_rep_choice

        # ASSIGN but CHECK that all count_total and count_exp values are the SAME
        count_total = list(set(ID_count_total_dl[name]))[0]
        count_exp = list(set(ID_count_exp_dl[name]))[0]

        if len(set(ID_count_total_dl[name])) != 1:
            print('\n### WARNING! There are conflicting numbers of total and expected permutations. BUG?')
            count_total = ','.join(map(str, ID_count_total_dl[name]))
            count_exp = ','.join(map(str, ID_count_exp_dl[name]))

        count_assn = ','.join(map(str, ID_count_assn_dl[name]))

        # TODO COMEBACK - why?
        prop_assn_list = []
        for i, count in enumerate(ID_count_assn_dl[name]):  # for rep_name in rep_seq_names_list_sorted:  # No NA;
            # count = ID_count_assn_dl[name][rep_name]
            this_prop_assn = count / count_total
            prop_assn_list.append(this_prop_assn)

            # # add all reps (not NA)
            # for rep_name in rep_seq_names_list_sorted:  # no NA
            #     count = count_dict[rep_name]
            #     if count_assn_sum > 0:
            #         this_prop = count / count_assn_sum
            #         line_list.append(count)  # PROPORTION of times this SEQ assigned to this REP
            #         # line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
            #         ID_prop_ddl[name][rep_name].append(this_prop)
            #     else:
            #         line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
            #         ID_prop_ddl[name][rep_name].append('NA')
            #
            # # Add NA LAST for COUNTS ONLY
            # if count_total_sum > 0:
            #     # this_prop = count_dict['NA'] / count_total_sum
            #     line_list.append(count_dict['NA'])  # PROPORTION of times this SEQ assigned to this REP
            #     # line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
            #     # ID_prop_ddl[name]['NA'].append(this_prop)
            # else:
            #     line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
            #     # ID_prop_ddl[name]['NA'].append('NA')

        prop_assn = ','.join(map(str, prop_assn_list))
        prop_assn_mean = np.mean(prop_assn_list)
        # count_assn_mean = np.mean(ID_count_assn_dl[name])

        # BEGIN line list
        line_start_list = [name, count_total, count_exp]

        # for rep_name, count in count_dict.items():
        #     if count_sum > 0:
        #         this_prop = count / count_sum
        #         line_list.append(count)  # PROPORTION of times this SEQ assigned to this REP
        #         # line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
        #         ID_prop_ddl[name][rep_name].append(this_prop)
        #     else:
        #         line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
        #         ID_prop_ddl[name][rep_name].append('NA')
        #
        #     for i, count in enumerate(count_list):
        #         if count_sum > 0:
        #             this_prop = count / count_sum
        #             line_list.append(this_prop)  # PROPORTION of times this SEQ assigned to this REP
        #             ID_prop_ddl[name][rep_seq_names_list_sorted[i]].append(this_prop)
        #         else:
        #             line_list.append('NA')
        #             ID_prop_ddl[name][rep_seq_names_list_sorted[i]].append('NA')

        line_start = '\t'.join(map(str, line_start_list)) + '\t'

        # # add NA first for ALL
        # if count_total_sum > 0:  # always true
        #     this_prop = count_dict['NA'] / count_total_sum
        #     line_list.append(count_dict['NA'])
        #     ID_prop_ddl[name]['NA'].append(this_prop)
        # else:
        #     line_list.append('NA')  # this is for a numerical VALUE, not assigned clade name
        #     ID_prop_ddl[name]['NA'].append('NA')

        # PRINT
        out_file_counts_hdl.write(line_start + count_assn + '\t' + '\t'.join(map(str, count_list)) + '\n')
        out_file_props_hdl.write(line_start + prop_assn + '\t' + '\t'.join(map(str, prop_list)) + '\n')
        out_file_means_hdl.write(line_start + f'{prop_assn_mean}\t{this_NA_mean}\t' +
                                 f'{closest_rep}\t{closest_rep_d}\t' +
                                 f'{top_rep_choice}\t{top_rep_prop}\t{rep_assigned}\t' +
                                 '\t'.join(map(str, mean_list)) + '\t' +
                                 '\t'.join(map(str, dist_list)) +
                                 '\n')

        # out_file_means_hdl.write(line_start + str(prop_assn_mean) + '\t' + '\t'.join(map(str, mean_list)) + '\n')
        # out_file_counts_hdl.write(line_start + '\t'.join(map(str, count_list)) + '\n')
        # out_file_props_hdl.write(line_start + '\t'.join(map(str, prop_list)) + '\n')
        # out_file_means_hdl.write(line_start + '\t'.join(map(str, mean_list)) + '\n')

    # CLOSE output files
    out_file_counts_hdl.close()
    out_file_props_hdl.close()
    out_file_means_hdl.close()

    # -------------------------------------------------------------------------
    # TODO check for convergence between trees

    # -------------------------------------------------------------------------
    # DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time/60, ndigits=1)} minutes; {round(elapsed_time/60/60, ndigits=1)} hours)')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
