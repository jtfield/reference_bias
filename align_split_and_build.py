#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib
import shutil
import pandas as pd
import subprocess
import datetime
import dateutil
# print(sys.path)
path_root = pathlib.Path(__file__).parents[0]
sys.path.append(str(path_root) + '/modules')
# print(sys.path)

from modules.find_seqs import *

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_dir', default=False, help='input directory of alignments option.')
    parser.add_argument('--output_dir', default=False, help='output directory of alignments option.')
    parser.add_argument('--specific_taxon', default=False, help='taxon of interest to pull from the alignments in the align_dir.')
    return parser.parse_args()

def main():
    args = parse_args()

    list_of_files = os.listdir(args.align_dir)

    os.mkdir(args.output_dir)

    if args.specific_taxon == False:
        # Make directory to hold new alignments
        for num, file in enumerate(list_of_files):

            # suffix = '_removed.aln'
            # pattern = '>\w+\s'

            # get names from the first file
            if num == 1:

                read_file_seqs = open(args.align_dir + '/' + file).read()
                # print(read_file_seqs)

                matches = re.findall(pattern, read_file_seqs)

                if matches:

                    for num2, name in enumerate(matches):

                        # LINE ONLY FOR TESTING
                        # Prevent building an alignment for every taxon
                        if num2 < 3:

                            name = name.strip('\n').strip('>')
                            print(name)

                            find_specific_seqs(args.align_dir, list_of_files, args.output_dir, name)

    else:

        find_specific_seqs(args.align_dir, list_of_files, args.output_dir, args.specific_taxon)





if __name__ == '__main__':
    main()
