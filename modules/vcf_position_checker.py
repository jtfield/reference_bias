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
    parser = argparse.ArgumentParser(prog='vcf_position_checker.py', \
        description='Run after combine_cov_vcfs.py to replace the sra numbers in the file names with the sample IDs.')
    parser.add_argument('--aln_file', help='vcf_file.')
    return parser.parse_args()

def main():
    args = parse_args()

    #read_vcf = pd.read_table(args.vcf, sep='\s+', skiprows=21, escapechar='#')

    read_file = open(args.aln_file, 'r').read()
    columns = ['pos', 'prior_gaps']
    # gap_position_df = pd.DataFrame(columns=columns)

    split_file = read_file.split('>')
    for chunk in split_file:
        if len(chunk) > 1:
            gap_position_df = pd.DataFrame(columns=columns)
            split_name_and_seq = chunk.split('\n', 1)
            name = split_name_and_seq[0]
            seq = split_name_and_seq[1]
            print(name)


            gap_count = 0
            for num, char in enumerate(seq):
                if char == '-':
                    gap_count+=1
                    gap_position_df.loc[len(gap_position_df.index)] = [num, gap_count]

            output_table_name = name + '_gap_tracker.csv'
            gap_position_df.to_csv(output_table_name)

    # print(gap_position_df)




    # print(read_vcf.columns)
    #
    # rows_to_print = [2321299,2321300,2321301]
    # test_rows = [2,3,4]
    # focus_row = 2321300
    #
    # pos_row_1 = read_vcf.loc[read_vcf['POS'] == rows_to_print[0]]
    # pos_row_2 = read_vcf.loc[read_vcf['POS'] == rows_to_print[1]]
    # pos_row_3 = read_vcf.loc[read_vcf['POS'] == rows_to_print[2]]
    #
    # print(pos_row_1)
    # print(pos_row_2)
    # print(pos_row_3)
    #
    # row_count = 0
    # # print(read_vcf.iloc[[2321299]])
    # # print(read_vcf.iloc[[2321300]])
    # # print(read_vcf.iloc[[2321301]])
    #
    # print(len(read_vcf))
    #
    # for idx, row in read_vcf.iterrows():
    #     row_count+=1
    #     if row['POS'] != row_count:
    #         print(row)
    #         row_count == row['POS']






if __name__ == '__main__':
    main()
