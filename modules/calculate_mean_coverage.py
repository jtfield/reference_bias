#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='calculate_mean_coverage.py', \
        description='calculates the average coverage of a sequence.')
    parser.add_argument('--csv_dir', help='directory of total coverage count csv files.')
    parser.add_argument('--distance_table', help='table of pairwise distances of taxa.')
    parser.add_argument('--output_csv', default='average_query_and_ref_coverage.csv', help='output csv name.')
    return parser.parse_args()

def main():
    args = parse_args()

    list_of_cov_counts = os.listdir(args.csv_dir)

    distance_table = pd.read_csv(args.distance_table)

    distance_table.columns = distance_table.columns.str.replace(' ', '_')

    distance_table['tax_names'] = distance_table['tax_names'].str.replace(' ', '_')

    distance_table = distance_table.set_index('tax_names')

    taxa_list = distance_table.columns

    columns = ['query_taxon', 'ref_taxon', 'dist_to_ref', 'avg_coverage']

    output_df = pd.DataFrame(columns=columns)

    for num, cov_file in enumerate(list_of_cov_counts):
        print(cov_file)
        split_file_name = cov_file.split('_query_')
        ref = split_file_name[0].replace('total_unambig_counts_ref_','')
        query = split_file_name[1].replace('.csv','')
        # print(ref)
        # print(query)
        # print(num)
        avg_cov = calc_cov_avg(args.csv_dir + '/' + cov_file)

        distance = distance_table.at[ref, query]

        data_row = [query, ref, distance, avg_cov]

        row = pd.Series(data_row, index=output_df.columns)

        output_df = output_df.append(row, ignore_index=True)

    # print(output_df)
    output_df.to_csv(args.output_csv)




def calc_cov_avg(cov_file):
    coverage_df = pd.read_csv(cov_file)

    coverage_df= coverage_df[coverage_df['count'] != 0]

    coverage_df['cov_times_count'] = coverage_df['coverage'] * coverage_df['count']
    sum_cov_bases = coverage_df['cov_times_count'].sum()
    # print(sum_cov_bases)

    seq_len = coverage_df['count'].sum()
    # print(seq_len)

    cov_avg = sum_cov_bases / seq_len
    print(cov_avg)
    return cov_avg



if __name__ == '__main__':
    main()
