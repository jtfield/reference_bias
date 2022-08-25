#!/usr/bin/env python3
"""example usage:
python diff_counter.py seq1.fas seq2.fas
Both sequences much be in fasta
"""

import os
import argparse
import re
import pathlib
import subprocess
import shutil
import pandas as pd
from itertools import islice

def parse_args():
    parser = argparse.ArgumentParser(prog='coverage_analysis.py', \
        description='Run after combine_cov_vcfs.py to replace the sra numbers in the file names with the sample IDs.')
    parser.add_argument('--vcf_dir', help='Directory of multiple coverage vcf files.')
    parser.add_argument('--diffs_dir', help='directory of differences files.')
    parser.add_argument('--gap_files_dir', help='directory of files tracking gap positions in each ref files.')
    parser.add_argument('--output_csv', default='unambigous_errors_coverage.csv', help='table of coverage of each unambiguous error analyzed.')
    # parser.add_argument('--output_coverage_dir', default='coverage_analysis_results', help='directory of updated unambiguous differences files that now include coverage data for each ref_base.')
    return parser.parse_args()

def main():
    args = parse_args()

    vcf_list = os.listdir(args.vcf_dir)
    diffs_list = os.listdir(args.diffs_dir)
    missing_diff_files = []

    # output_dict = {}
    # for i in range(0,500):
    #     output_dict[i] = 0

    for file_name in vcf_list:
        split_file_name = file_name.split("_query_")
        ref = split_file_name[0].replace('ref_','')
        query = split_file_name[1].replace('.vcf','')

        unfound_unambig_diffs_file = 'unambig_differences_query_' + query + '_ref_' + ref + '.csv'


        if unfound_unambig_diffs_file in diffs_list:
            print(unfound_unambig_diffs_file)
            ref_gaps = pd.read_csv(args.gap_files_dir + '/' + ref + '_gap_tracker.csv')
            read_diffs = pd.read_csv(args.diffs_dir + '/' + unfound_unambig_diffs_file)
            # read_vcf = pd.read_table(args.vcf_dir + '/' + file_name, sep='\s+', skiprows=21, escapechar='#')
            read_vcf = smart_read_vcf(args.vcf_dir + '/' + file_name)

            # print(ref_gaps.columns)
            # print(read_diffs.columns)
            # print(read_vcf.columns)

            #get_real_position(read_vcf, ref_gaps, read_diffs, output_dict, unfound_unambig_diffs_file)
            get_real_position(read_vcf, ref_gaps, read_diffs, unfound_unambig_diffs_file)

            # exit()

        elif unfound_unambig_diffs_file not in diffs_list:
            missing_diff_files.append(unfound_unambig_diffs_file)

    print(missing_diff_files)
    # print(output_dict)

    # output_df = pd.DataFrame(list(output_dict.items()), columns = ['coverage', 'count'])
    # output_df.to_csv('unambigous_diffs_coverage_counts.csv')


#def get_real_position(vcf, gaps_df, bai_df, output_dict, file_name):
def get_real_position(vcf, gaps_df, bai_df, file_name):

    output_dict = {}
    for i in range(0,500):
        output_dict[i] = 0

    vcf_columns = vcf.columns
    bai_columns = bai_df.columns
    # print(bai_columns)

    checked_base_count = 0
    for idx, row in bai_df.iterrows():
        bai_pos = row['positions']
        # print(bai_pos)

        earlier_gap_rows = gaps_df.loc[gaps_df['pos'] <= bai_pos]
        # print(earlier_gap_rows)
        num_gaps = len(earlier_gap_rows)

        updated_vcf_position = (bai_pos - num_gaps) + 1
        # correct_vcf_row = vcf.loc[vcf['POS'] == updated_vcf_position]
        correct_vcf_row = vcf.index[vcf['POS'] == updated_vcf_position].tolist()

        try:
            vcf_base = vcf.at[correct_vcf_row[0], 'REF']
            error_base = row[bai_columns[3]]
            assert vcf_base.upper() == error_base.upper()
            info_column = vcf.at[correct_vcf_row[0], 'INFO']

            split_info = info_column.split(';', 1)
            dp_cov = split_info[0]
            split_cov = dp_cov.split('=')
            coverage_num = split_cov[1]
            # print(coverage_num)
            output_dict[int(coverage_num)] +=1

        except IndexError:
            print("error at vcf indexing")
            print(correct_vcf_row)
            print(correct_vcf_row[0])

    output_df = pd.DataFrame(list(output_dict.items()), columns = ['coverage', 'count'])
    output_df.to_csv('coverage_results_' + file_name)



def smart_read_vcf(unread_vcf_obj):
    number_double_hash = 0
    with open(unread_vcf_obj) as vcf_file:
        head = list(islice(vcf_file, 50))
    for line in head:
        if "##" in line:
            number_double_hash+=1

    read_vcf = pd.read_table(unread_vcf_obj, sep='\s+', skiprows=number_double_hash, escapechar='#')

    return read_vcf


if __name__ == '__main__':
    main()
