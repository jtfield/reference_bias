#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import dendropy
import pathlib
path_root = pathlib.Path(__file__).parents[0]
sys.path.append(str(path_root) + '/modules')
from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='bias analysis', \
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--align_file', default=False, help='input starting alignment option.')
    parser.add_argument('--align_dir', default=False, help='input directory of alignments option.')
    parser.add_argument('--suffix', default='_removed.aln', help='input suffix of new aligments option.')
    parser.add_argument('--output_csv', default='ref_bias_comparison.csv', help='output CSV file path and name.')
    parser.add_argument('--br_csv', default='branch_length.csv', help='output branch length CSV file path and name.')

    return parser.parse_args()

def main():
    args = parse_args()

    # Build matrix of branch lengths between taxa in the tree
    df = build_branch_length_matrix(args.tree_file)

    # # print(df.loc['Aphrodroma brevirostris Nunn KGP 1'])
    # print('Aphrodroma brevirostris Nunn KGP 1', 'Puffinus nativitatis USNM 613922')
    # print(df.columns.tolist())
    # print(df.index.tolist())
    # test_check = df.at['Aphrodroma brevirostris Nunn KGP 1', 'Puffinus nativitatis USNM 613922']
    # print(test_check)

    df.to_csv(args.br_csv)

    # Make sequence comparisons
    final_data_list = prepare_seq_comparison(args.align_file, args.align_dir, args.suffix, df)


    filled_final_df = build_basic_comparison_df(final_data_list, args.output_csv)




if __name__ == '__main__':
    main()
