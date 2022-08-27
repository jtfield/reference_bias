#!/usr/bin/env python3
"""example usage:
python diff_counter.py seq1.fas seq2.fas
Both sequences much be in fasta
"""

import os
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(prog='calculate_sum_identical.py', \
        description='takes the combined summary file and a directory of coverage error counts.')
    parser.add_argument('--csv_dir', help='Directory of multiple coverage vcf files.')
    parser.add_argument('--comparison_summary', help='file of combined comparison values.')
    parser.add_argument('--coverage_summary', help='file of combined coverage and error values.')
    return parser.parse_args()

def main():
    args = parse_args()

    counts = []

    list_of_csvs = os.listdir(args.csv_dir)

    read_comparison = pd.read_csv(args.comparison_summary)

    read_coverage = pd.read_csv(args.coverage_summary)

    for csv_file in list_of_csvs:

        split_file_on_ref = csv_file.split('_ref_')
        query = split_file_on_ref[0].replace('coverage_results_unambig_differences_query_', '')
        ref = split_file_on_ref[1].replace('.csv','')
        # print(query)
        # print(ref)
        # print("###")

        row_with_count = read_comparison.loc[(read_comparison['taxon_name'] == query) & (read_comparison['ref_name'] == ref)]
        ident_bases = row_with_count['all_unambig_identical'].tolist()
        counts.append(ident_bases[0])

    total_ident = sum(counts)
    updated_df = read_coverage.assign(total_ident=total_ident)
    updated_df = updated_df.assign(error_rate=(updated_df['count'] / updated_df['total_ident']))

    updated_df.to_csv('full_counts_combined_unambig_coverage_counts.csv')


if __name__ == '__main__':
    main()
