#!/usr/bin/env python3

import os
import dendropy
import argparse
from dendropy.calculate import treecompare
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree_dir', help='Directory of phylogeny files in newick format')
    parser.add_argument('--tree_file', help='Phylogeny file in newick format')
    parser.add_argument('--br_l', default='unweighted', help='Toggle for weighted or unweighted RF distance. DEFAULT: unweighted')
    parser.add_argument('--output_csv', default='unweighted_rf_distances.csv', help='csv table of output RF distances. change name if using Weighted RF.')
    return parser.parse_args()

def prune_trees_to_match(t1, t2):
    """requires two trees with a common taxon namespace"""
    tree_1_taxa = set()
    tree_2_taxa = set()

    for tip in t1.leaf_node_iter():
        tree_1_taxa.add(tip.taxon)

    for tip in t2.leaf_node_iter():
        tree_2_taxa.add(tip.taxon)

    shared_taxa = tree_1_taxa.intersection(tree_2_taxa)

    assert(len(shared_taxa) >= 1)
    # print("These two tree have {s} shared taxa".format(s=len(shared_taxa)))

    t1.retain_taxa(shared_taxa)
    t2.retain_taxa(shared_taxa)
    return(t1, t2)

def main():
    args = parse_args()

    list_of_trees = os.listdir(args.tree_dir)

    rf_values = []

    output_df = pd.DataFrame(columns=list_of_trees)
    # print(output_df)

    for tree_1 in list_of_trees:

        tns = dendropy.TaxonNamespace()
        # # ensure all trees loaded use common namespace
        t1 = dendropy.Tree.get_from_path(
        src=args.tree_dir + '/' + tree_1,
        schema='newick',
        taxon_namespace=tns, preserve_underscores=True)
        t2 = dendropy.Tree.get_from_path(
        src=args.tree_file,
        schema='newick',
        taxon_namespace=tns, preserve_underscores=True)

        t1, t2 = prune_trees_to_match(t1, t2)
        rf = 100000000
        if args.br_l == 'unweighted':
            rf = treecompare.symmetric_difference(t1, t2)
            # print("Unweighted symmetric difference is {}".format(treecompare.symmetric_difference(t1, t2)))
        elif args.br_l == 'weighted':
            rf = treecompare.weighted_robinson_foulds_distance(t1, t2)
            # print("Weighted symmetric difference is {}".format(treecompare.weighted_robinson_foulds_distance(t1, t2)))

        # output_df.at[file_name, tree_1] = rf
        rf_values.append(rf)

    output_df.loc['true_tree'] = rf_values
    print(output_df)
    output_df.to_csv(args.output_csv)
    stats = output_df.describe(include='all')
    print(stats)
    print(output_df.max().max())
    print(output_df.min().min())
    print(output_df.mean(axis=1))


if __name__ == '__main__':
    main()
