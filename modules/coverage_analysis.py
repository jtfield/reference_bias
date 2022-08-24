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

def parse_args():
    parser = argparse.ArgumentParser(prog='coverage_analysis.py', \
        description='Run after combine_cov_vcfs.py to replace the sra numbers in the file names with the sample IDs.')
    parser.add_argument('--vcf_dir', help='Directory of multiple coverage vcf files.')
    parser.add_argument('--diffs_dir', help='directory of differences files.')
    parser.add_argument('--gap_files_dir', help='directory of files tracking gap positions in each ref files.')
    # parser.add_argument('--output_coverage_dir', default='coverage_analysis_results', help='directory of updated unambiguous differences files that now include coverage data for each ref_base.')
    return parser.parse_args()

def main():
    args = parse_args()

    vcf_list = os.listdir(args.vcf_dir)
    diffs_list = os.listdir(args.diffs_dir)
    missing_diff_files = []

    for file_name in vcf_list:
        split_file_name = file_name.split("_query_")
        ref = split_file_name[0].replace('ref_','')
        query = split_file_name[1].replace('.vcf','')

        unfound_unambig_diffs_file = 'unambig_differences_query_' + query + '_ref_' + ref + '.csv'


        if unfound_unambig_diffs_file in diffs_list:
            print(unfound_unambig_diffs_file)
            ref_gaps = pd.read_csv(args.gap_files_dir + '/' + ref + '_gap_tracker.csv')
            read_diffs = pd.read_csv(args.diffs_dir + '/' + unfound_unambig_diffs_file)
            read_vcf = pd.read_table(args.vcf_dir + '/' + file_name, sep='\s+', skiprows=21, escapechar='#')

            print(ref_gaps.columns)
            print(read_diffs.columns)
            print(read_vcf.columns)

            # coverage_check(read_vcf, read_diffs)

        elif unfound_unambig_diffs_file not in diffs_list:
            missing_diff_files.append(unfound_unambig_diffs_file)

    print(missing_diff_files)


def coverage_check(cov_vcf, diff_csv):
    for idx, diff_line in diff_csv.iterrows():
        print("check_position##############################################")
        diff_position = diff_line['positions']
        print(diff_position)
        print(diff_line)

        vcf_pos = cov_vcf.loc[cov_vcf['POS'] == diff_position]
        print(vcf_pos['POS'])
        print(vcf_pos['REF'])
        print(vcf_pos['ALT'])




if __name__ == '__main__':
    main()
