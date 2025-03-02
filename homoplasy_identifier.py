#!/usr/bin/env python3

"""
Purpose: Determine nuc or aa sites with homoplasies from a tree and an alignment of sequences
Author : Chase W. Nelson <chase.nelson@nih.gov>
Cite   : https://github.com/chasewnelson/
Date   : 2024-11-09
"""

import argparse
import os
import pandas as pd
import re
import sys
import time
from ete3 import Tree
from ete3.parser import newick
from Bio import SeqIO
from collections import defaultdict
from pprint import pprint
from typing import Any, Dict, List, NamedTuple, TextIO  # Set,

usage = """# -----------------------------------------------------------------------------
homoplasy_identifier.py - Determine nuc or aa sites with homoplasies from a tree and an alignment of sequences
# -----------------------------------------------------------------------------
For DOCUMENTATION, run:
    $ homoplasy_identifier.py --help
    $ pydoc ./homoplasy_identifier.py
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
EXAMPLE:
    $  homoplasy_identifier.py --fasta_file seqs.fasta --tree_file raxml.bestTree --out_file raxml.bestTree-homoplasies.txt --seq_type n > raxml.bestTree-homoplasies.log
# ----
"""


class Args(NamedTuple):
    """ Command-line arguments """
    fasta_file: TextIO
    tree_file: TextIO
    seq_type: str
    out_file: str


# -----------------------------------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Determine nuc or aa sites with homoplasies from a tree and an alignment of sequences. HELP: homoplasy_identifier.py --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Rename "optional" arguments
    parser._optionals.title = 'Named arguments'

    # -------------------------------------------------------------------------
    # REQUIRED
    parser.add_argument('-f',
                        '--fasta_file',
                        metavar='FILE',
                        help='Multiple sequence alignment, nucleotide or amino acid [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))  # seq/HPV16_PAP_20200813.N-30.fasta

    parser.add_argument('-t',
                        '--tree_file',
                        metavar='FILE',
                        help='Tree in Newick format [REQUIRED]',
                        required=True,
                        nargs=1,
                        type=argparse.FileType('rt'))

    parser.add_argument('-s',
                        '--seq_type',
                        metavar='str',
                        help='Sequence type (n=nucleotide; a=amino acid) [REQUIRED]',
                        required=True,
                        choices=['n', 'a'],  # (n)ucleotide or (a)mino acid
                        type=str)

    parser.add_argument('-o',
                        '--out_file',
                        metavar='str',
                        help='Output file name [REQUIRED]',
                        required=True,
                        type=str)

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # VALIDATE arguments

    # ensure out_file not present - OK to check for the default one too, even though it'll be changed (will check again)
    if os.path.isfile(args.out_file):
        parser.error(f'\n### ERROR: out_file="{args.out_file}" already exists')

    # # if tree is file, convert to string
    # if os.path.isfile(args.tree):
    #     args.tree = open(args.tree).read().rstrip()

    # die if tree contains anything NON-newick; should simply fail to initialize with Tree()
    # re_NOT_newick = re.compile(r'[^\w\d.:;()\s\n,]')
    # if re_NOT_newick.search(args.tree):
    #     parser.error(f'\n### ERROR: tree="{args.tree}..." may only contain Newick-compatible characters')
    # if args.tree != '':
    #     try:
    #         _ = Tree(args.tree)
    #     except newick.NewickError:
    #         parser.error(f'\n### ERROR: tree="{args.tree}" does not conform to Newick format')

    return Args(fasta_file=args.fasta_file,
                out_file=args.out_file,
                seq_type=args.seq_type,
                tree_file=args.tree_file)


# -----------------------------------------------------------------------------
def main() -> None:
    """ Tell them they are walking around shining like the sun """
    start_time = time.time()

    # -------------------------------------------------------------------------
    # GATHER arguments
    args = get_args()
    # fasta_fh_list = args.fasta_file
    # tree_fh_list = args.tree_file
    fasta_fh = args.fasta_file[0]
    tree_fh = args.tree_file[0]
    seq_type = args.seq_type
    out_file = args.out_file

    # -------------------------------------------------------------------------
    # INITIALIZE OUTPUT
    print(usage)

    # -------------------------------------------------------------------------
    # REGEX & TUPLES

    # Move to check_homoplasy() function
    # # Nucleotides
    # DEFINED_NUCLEOTIDES = set("ACGTU")
    #
    # # Define valid nucleotide or amino acid characters
    # DEFINED_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
    #
    # # CHOOSE THE APPROPRIATE SET
    # DEFINED_CHARACTERS = None
    #
    # if seq_type == "n":
    #     DEFINED_CHARACTERS = DEFINED_NUCLEOTIDES
    # elif seq_type == "a":
    #     DEFINED_CHARACTERS = DEFINED_AMINO_ACIDS

    # -------------------------------------------------------------------------
    # INITIALIZE LOG
    print('# -----------------------------------------------------------------------------')
    print(f'LOG:command="{" ".join(sys.argv)}"')
    print(f'LOG:cwd="{os.getcwd()}"')
    print(f'LOG:fasta_file="{fasta_fh.name}"')
    print(f'LOG:tree_file="{tree_fh.name}"')
    # print(f'LOG:fasta_file="{",".join([fh.name for fh in fasta_fh_list])}"')
    # print(f'LOG:tree_file="{tree}"', flush=True)
    print(f'LOG:out_file="{out_file}"', flush=True)

    # -------------------------------------------------------------------------
    # TODO: CHECK FOR DUPLICATE SEQUENCE IDs?

    # -------------------------------------------------------------------------
    # INPUT FASTA SEQUENCE(S) and COUNT ALLELES FOR ALL GROUPS

    # Example usage:
    # tree_file = "path/to/tree_file.newick"
    # fasta_file = "path/to/alignment_file.fasta"
    check_homoplasy(tree_fh, fasta_fh, seq_type, out_file)

    # -------------------------------------------------------------------------
    # DONE message
    end_time = time.time()
    elapsed_time = end_time - start_time

    print('\n# -----------------------------------------------------------------------------')
    print(f'TIME ELAPSED: {round(elapsed_time, ndigits=1)} seconds '
          f'({round(elapsed_time / 60, ndigits=1)} minutes; {round(elapsed_time / 60 / 60, ndigits=1)} hours)')

    print('\n# -----------------------------------------------------------------------------')
    print('DONE')


# -----------------------------------------------------------------------------
def load_alignment(fasta_file):
    # Load the alignment and return a dictionary where each key is a taxon name and each value is a sequence
    alignment = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        alignment[record.id] = str(record.seq)
    return alignment


# -----------------------------------------------------------------------------
def check_homoplasy(tree_fh, fasta_file, seq_type, out_file):
    # Load the tree
    base_tree = Tree(tree_fh.name, format=1)

    # Load the alignment
    alignment = load_alignment(fasta_file)

    # Define valid characters based on seq_type
    VALID_CHARACTERS = None

    if seq_type == "n":
        VALID_CHARACTERS = set("ACGTU")  # nucleotides
    elif seq_type == "a":
        VALID_CHARACTERS = set("ACDEFGHIKLMNPQRSTVWY")  # amino acids

    # Ensure all taxa in the alignment are in the tree
    alignment_taxa = set(alignment.keys())
    tree_taxa = set(leaf.name for leaf in base_tree.get_leaves())
    if not alignment_taxa.issubset(tree_taxa):
        raise ValueError("All taxa in the alignment must be present in the tree.")

    # Get the length of the alignment (all sequences should be the same length)
    alignment_length = len(next(iter(alignment.values())))

    # Prepare a DataFrame to store results
    homoplasy_results = pd.DataFrame(index=range(1, alignment_length + 1), columns=list(VALID_CHARACTERS))

    # # Track positions with homoplasy
    # homoplasy_sites = []
    # position_to_num_homoplasies = defaultdict(int)
    # position_to_homoplasies = defaultdict(str)
    # position_to_states = defaultdict(int)

    # Iterate over each site in the alignment
    for position in range(alignment_length):
        # DEBUG
        # print(f'position={position + 1}')

        tree = base_tree.copy()

        # Collect nucleotides or amino acids at this position, excluding gaps and ambiguous characters
        site_data = defaultdict(list)
        invalid_taxa = []
        for taxon, sequence in alignment.items():
            character = sequence[position]
            if character in VALID_CHARACTERS:
                site_data[character].append(taxon)
            else:
                invalid_taxa.append(taxon)  # Mark taxa with invalid characters for pruning

        # Determine taxa to keep in the tree (those not in invalid_taxa)
        taxa_to_keep = [taxon for taxon in tree.get_leaf_names() if taxon not in invalid_taxa]

        # Check if all taxa are invalid (i.e., no valid taxa remain)
        if not taxa_to_keep:
            # If all taxa are invalid, set all columns to "NA" for this site and skip further checks
            homoplasy_results.loc[position + 1] = "NA"
            continue

        # Prune the tree to keep only the valid taxa
        tree.prune(taxa_to_keep)

        # Check homoplasy for each valid character at this site
        for character in VALID_CHARACTERS:
            if character not in site_data:
                homoplasy_results.at[position + 1, character] = "NA"  # Character not present
            else:
                taxa = site_data[character]
                if len(taxa) > 1:  # Only check monophyly for groups with more than one taxon
                    is_monophyletic, _, _ = tree.check_monophyly(values=taxa, target_attr="name")
                    homoplasy_results.at[position + 1, character] = not is_monophyletic  # True if homoplasy
                else:
                    homoplasy_results.at[position + 1, character] = False  # Single taxon is trivially monophyletic

        # Save results to a TSV file
        homoplasy_results.to_csv(out_file, sep='\t', index_label="site")

        # # # Skip monotypic sites (only one nucleotide type across all taxa)
        # # if len(site_data) == 1:
        # #     continue  # No homoplasy possible if only one nucleotide at this site
        #
        # # Check monophyly for each nucleotide group
        # has_homoplasy = False
        # num_homoplasies = 0
        # for character, taxa in site_data.items():
        #     # print(f"character={character}")
        #     # pprint(taxa)
        #
        #     position_to_states[position + 1] += 1
        #
        #     if len(taxa) > 1:  # Only check monophyly for groups with more than one taxon
        #         # Check if the group of taxa with this nucleotide is monophyletic
        #         # is_monophyletic, _, _ = tree.check_monophyly(values=taxa, target_attr="name")
        #         is_monophyletic, _, _ = tree.check_monophyly(values=taxa, target_attr="name")
        #
        #         if not is_monophyletic:
        #             has_homoplasy = True
        #             num_homoplasies += 1
        #             position_to_num_homoplasies[position + 1] += 1
        #             position_to_homoplasies[position + 1] += character
        #             # break  # If one nucleotide group is non-monophyletic, mark the site as having homoplasy
        #
        # # Record the site if homoplasy is detected
        # if has_homoplasy:
        #     homoplasy_sites.append(position + 1)  # Use 1-based index for positions

        # DEBUG: stop early
        # if position > 84:
        #     break

    # # Output the result
    # if homoplasy_sites:
    #     print(f"Sites with homoplasy: {homoplasy_sites}")
    # else:
    #     print("No homoplasy detected at any site.")

    # # Print the results
    # for position in sorted(position_to_states.keys()):
    #     print(f"{position}: " +
    #           f"states={position_to_states[position]} | " +
    #           f"num_homoplasies={position_to_num_homoplasies[position]} | " +
    #           f"homoplasies={position_to_homoplasies[position]}")


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()

