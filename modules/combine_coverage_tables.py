#!/usr/bin/env python3
"""example usage:
python diff_counter.py seq1.fas seq2.fas
Both sequences much be in fasta
"""

import os
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(prog='combine_coverage_tables.py', \
        description='Combines the multiple coverage tables into a single table for analysis.')
    parser.add_argument('--csv_dir', help='Directory of multiple coverage vcf files.')
    # parser.add_argument('--output_coverage_dir', default='coverage_analysis_results', help='directory of updated unambiguous differences files that now include coverage data for each ref_base.')
    return parser.parse_args()

def main():
    args = parse_args()

    list_of_table_files = os.listdir(args.csv_dir)

    list_of_dfs = []

    for table in list_of_table_files:
        current_df = pd.read_csv(args.csv_dir + '/' + table)
        list_of_dfs.append(current_df)
        # print(current_df)

    combined_table = pd.concat(list_of_dfs).groupby('coverage').sum().reset_index()

    combined_table.drop(combined_table.columns[[1]], axis=1, inplace=True)

    combined_table.to_csv('combined_unambig_coverage_counts.csv')




if __name__ == '__main__':
    main()
