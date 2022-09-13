#!/usr/bin/env python3

import os
import argparse
import dendropy
import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm


def parse_args():
    parser = argparse.ArgumentParser(prog='tree_length_figure.py', \
        description='generates a figure based on the length of the true tree compared to a set of reference based trees.')
    parser.add_argument('--ref_trees_dir', default=False, help='input ref-based phylogeny dir.')
    parser.add_argument('--true_tree_file', default=False, help='input "true" phylogeny option.')
    parser.add_argument('--true_aln_file', default=False, help='input "true" alignment option.')

    return parser.parse_args()

def main():
    args = parse_args()

    list_of_trees = os.listdir(args.ref_trees_dir)

    list_of_taxa = [name.replace('_ref_tre','') for name in list_of_trees]
    # print(list_of_taxa)

    columns = ['ref_taxon', 'tree_length', 'true_tree_length', 'tree_length_difference']

    dist_df = pd.DataFrame(columns=columns, index=list_of_taxa)

    true_tree = dendropy.Tree.get(path=args.true_tree_file, schema="newick")
    print(true_tree.length())

    for tree_file in list_of_trees:
        # tree_file = tree_file.replace('_ref_tre','')
        tree = dendropy.Tree.get(path=args.ref_trees_dir + '/' + tree_file, schema="newick")
        tree_file = tree_file.replace('_ref_tre','')
        dist_df.at[tree_file, 'ref_taxon'] = '_'.join(tree_file.split('_',2)[:2])
        dist_df.at[tree_file, 'tree_length'] = tree.length()
        dist_df.at[tree_file, 'true_tree_length'] = true_tree.length()
        dist_diff = tree.length() - true_tree.length()
        dist_df.at[tree_file, 'tree_length_difference'] = dist_diff

    #print(dist_df)

    updated_df = get_missing_data(args.true_aln_file, dist_df)

    sorted_length_diffs = updated_df.sort_values(by = ['tree_length_difference'])
    # print(sorted_length_diffs)

    sorted_length_diffs.to_csv('missing_data_vs_ref_tree_length_change.csv')


    # plt.figure(figsize=(10,10), dpi=100)
    #
    # sorted_length_diffs.reset_index(inplace=True)
    #
    # plt.hlines(y=sorted_length_diffs.ref_taxon, xmin = 0, xmax = sorted_length_diffs.tree_length_difference, alpha=0.4, linewidth=5)
    #
    # plt.gca().set(ylabel='Reference_taxon', xlabel='Length difference to true tree')
    #
    # plt.tight_layout()
    #
    # plt.show()

    #plt.savefig('tree_length_differences_to_true.png', bbox_inches='tight')


def get_missing_data(align, df):
    if align != False:

        missing_data = ['N', '-']

        df['missing_data_count'] = 0
        # print(df)

        read_align = open(align, 'r').read()

        split_align = read_align.split('>')

        for name_and_seq in split_align:
            if len(name_and_seq) > 0:
                split_name_and_seq = name_and_seq.split('\n', 1)

                name = split_name_and_seq[0]
                seq = split_name_and_seq[1]
                # print(name)

                missing_data_count = 0
                for nuc in seq:
                    if nuc.upper() in missing_data:
                        missing_data_count+=1

                df.at[name, 'missing_data_count'] = missing_data_count

    # print(df)
    return df











if __name__ == '__main__':
    main()
