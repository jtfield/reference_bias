#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import pathlib
import pandas as pd
# path_root = pathlib.Path(__file__).parents[0]
# sys.path.append(str(path_root) + '/modules')
# from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='sequence renamer', \
        description='Renames all alignments in a dir using, replacing the SRA run numbers used by Extensiphy. \
        Numbers are replaced with scientific names from a csv.')
    parser.add_argument('--csv_file', default=False, help='input phylogeny option.')
    parser.add_argument('--align_dir', default=False, help='input directory of alignments option.')
    parser.add_argument('--output_dir', default=False, help='output directory of alignments option.')
    parser.add_argument('--taxon_name_column', default="SampleName", help='the name of the column that has the sample names used in the csv.')

    return parser.parse_args()

def main():
    args = parse_args()

    # read csv
    metadata_csv = pd.read_csv(args.csv_file)

    # List alignments in alignment directory
    list_aligns = os.listdir(args.align_dir)

    # Loop over aligns, reading them and finding names
    for align in list_aligns:
        fix_names(align, args.align_dir, metadata_csv, args.taxon_name_column, args.output_dir)

def fix_names(alignment, align_dir, df, sample_name_column, outdir):
    """
    Find and replace names if the names in the alignment are the SRA numbers
    """

    # Read alignment
    align_data = open(align_dir + '/' + alignment, 'r').read()
    print(alignment)

    is_dir = os.path.isdir(outdir)

    if is_dir == False:
        os.mkdir(outdir)

    # Open output file
    output = open(outdir + '/' + alignment, 'w')

    # split alignment on the carrot
    # This means each chunk of the file contains a taxon name and a sequence
    # separated by a newline character
    split_align = align_data.split('>')

    # loop over each chunk and do stuff as long as the chunk actually has data
    for chunk in split_align:
        if len(chunk) > 0:

            # Split taxon name from seq
            split_chunk = chunk.split('\n', 1)

            #assign name and seq to variables
            name = split_chunk[0]
            seq = split_chunk[1]

            # check if name is in DF under the Run column
            if name in df['Run'].values:

                # If the name(sra number) is in the run column, get the row number
                df_row_of_sra = df.index[df.Run == name].tolist()[0]
                # print(df_row_of_sra)

                # Use the column name and the row number (index) to find the sample name
                # or whatever the name is that you want to replace the SRA number with
                sci_name = df.at[df_row_of_sra, sample_name_column]

                # print(sci_name)

                # Replace SRA number with new name and write to file
                output.write('>' + sci_name)
                output.write('\n')
                output.write(seq)
                output.write('\n')









if __name__ == '__main__':
    main()
