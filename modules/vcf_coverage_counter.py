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
    parser.add_argument('--output_csv', default='unambigous_coverage_levels.csv', help='table of coverage of each unambiguous error analyzed.')
    # parser.add_argument('--output_coverage_dir', default='coverage_analysis_results', help='directory of updated unambiguous differences files that now include coverage data for each ref_base.')
    return parser.parse_args()

def main():
    args = parse_args()

    vcf_list = os.listdir(args.vcf_dir)

    output_dict = {}
    for i in range(0,500):
        output_dict[i] = 0

    for file_name in vcf_list:
        print(file_name)

        read_vcf = smart_read_vcf(args.vcf_dir + '/' + file_name)

        get_real_position(read_vcf, output_dict)

        print(output_dict)


def get_real_position(vcf, cov_dict):
    vcf_columns = vcf.columns

    for idx, row in vcf.iterrows():

        ambig_bases = ['N', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X']

        ref = row['REF']
        alt = row['ALT']

        if 'N' not in alt and ref not in ambig_bases:

            info = row['INFO']

            split_info = info.split(';', 1)
            dp_cov = split_info[0]
            split_cov = dp_cov.split('=')
            coverage_num = split_cov[1]
            # print(coverage_num)
            cov_dict[int(coverage_num)] +=1





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
