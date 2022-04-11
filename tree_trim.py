#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import dendropy
import pathlib
# path_root = pathlib.Path(__file__).parents[0]
# sys.path.append(str(path_root) + '/modules')
# from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='remove taxa from tree', \
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--taxa_list', default=False, help='String that contains the names of each taxon to keep. (Example: "name1 name2 name3")')
    # parser.add_argument('--br_csv', default='branch_length.csv', help='output branch length CSV file path and name.')

    return parser.parse_args()

def main():
    args = parse_args()

    # names_to_keep = []

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    #
    # print(in_tree)
    #
    split_names = args.taxa_list.split()

    taxa_to_retain = set([taxon for taxon in in_tree.taxon_namespace if taxon.label in split_names])

    print(taxa_to_retain)


    filtered_tree = in_tree.extract_tree_with_taxa(taxa=taxa_to_retain)
    #
    #
    #
    print(filtered_tree.as_ascii_plot())




if __name__ == '__main__':
    main()
