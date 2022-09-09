#!/usr/bin/env python3

import os
import sys
import argparse
import pathlib
import pandas as pd
import dendropy
import re


def parse_args():
    parser = argparse.ArgumentParser(prog='tree_renamer.py', \
        description='Renames all alignments in a dir using, replacing the SRA run numbers used by Extensiphy. \
        Numbers are replaced with scientific names from a csv.')
    parser.add_argument('--csv_file', default=False, help='input phylogeny option.')
    parser.add_argument('--tree_dir', default=False, help='input phylogeny.')
    # parser.add_argument('--output_file', default='taxa_renamed.tre', help='output tree file name.')
    parser.add_argument('--taxon_name_column', default="SampleName", help='the name of the column that has the sample names used in the csv.')

    return parser.parse_args()

def main():
    args = parse_args()

    # read csv
    metadata_csv = pd.read_csv(args.csv_file)

    tree_dir = os.listdir(args.tree_dir)

    # read_tree_str = open(args.tree_file, 'r').read()

    for tree_file in tree_dir:

        tree_path = args.tree_dir + '/' + tree_file

        # new_tree = trim_non_sra_names(tree_path)

        fix_names(tree_path, metadata_csv, args.taxon_name_column)


def fix_names(tree_file, df, sample_name_column):
    """
    Find and replace names if the names in the tree are the SRA numbers
    """
    tree = open(tree_file, 'r').read()

    for idx, row in df.iterrows():
        # print(row['SampleName'])
        # print(row['Run'])
        # print("###")

        sra_regex = row['Run']
        sra_compile = re.compile(sra_regex)
        find_sra = re.findall(sra_compile, tree)
        if find_sra:
            tree = tree.replace(row['Run'], row['SampleName'])

        else:
            print("Problems on regex in: ", tree_file)

    tree_file_name = tree_file.rsplit('/', 1)
    updated_tree_file_name = 'names_updated_' + tree_file_name[1]
    output_tree = open(updated_tree_file_name, 'w')
    output_tree.write(tree)

    output_tree.close()


# def trim_non_sra_names(dp_tree):


            # # Split taxon name from seq
            # split_chunk = chunk.split('\n', 1)
            #
            # #assign name and seq to variables
            # name = split_chunk[0]
            # seq = split_chunk[1]
            #
            # # check if name is in DF under the Run column
            # if name in df['Run'].values:
            #
            #     # If the name(sra number) is in the run column, get the row number
            #     df_row_of_sra = df.index[df.Run == name].tolist()[0]
            #     # print(df_row_of_sra)
            #
            #     # Use the column name and the row number (index) to find the sample name
            #     # or whatever the name is that you want to replace the SRA number with
            #     sci_name = df.at[df_row_of_sra, sample_name_column]
            #
            #     # print(sci_name)
            #
            #     # Replace SRA number with new name and write to file
            #     output.write('>' + sci_name)
            #     output.write('\n')
            #     output.write(seq)
            #     output.write('\n')









if __name__ == '__main__':
    main()
