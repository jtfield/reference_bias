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
from os.path import exists

def parse_args():
    parser = argparse.ArgumentParser(prog='vcf_renamer.py', \
        description='Run after combine_cov_vcfs.py to replace the sra numbers in the file names with the sample IDs.')
    parser.add_argument('--vcf_dir', help='Directory of multiple vcf files you wish to rename.')
    parser.add_argument('--id_csv', default='coverage_vcf_files', help='csv with the IDs and sample names for use in renaming the vcf_files.')
    return parser.parse_args()

def main():
    args = parse_args()

    list_files = os.listdir(args.vcf_dir)
    # print(list_files)

    sample_df = pd.read_csv(args.id_csv)
    # print(sample_df.columns)
    # print(sample_df['Run'])
    # print(sample_df['SampleName'])

    for file_name in list_files:
        split_file = file_name.split('query_')
        front_of_file = split_file[0]
        query_sra = split_file[1].replace('.vcf','')
        print(query_sra)

        sra_loc = sample_df.loc[sample_df['Run'] == query_sra]
        found_sample_name = sra_loc['SampleName'].values[0]
        print(found_sample_name)

        new_file_name = front_of_file + 'query_' + found_sample_name + '.vcf'
        print(new_file_name)

        old_file_path = args.vcf_dir + '/' + file_name
        new_file_path = args.vcf_dir + '/' + new_file_name

        print(old_file_path)
        print(new_file_path)

        os.rename(old_file_path, new_file_path)






if __name__ == '__main__':
    main()
