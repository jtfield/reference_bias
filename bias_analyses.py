#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import dendropy
import pathlib
path_root = pathlib.Path(__file__).parents[0]
sys.path.append(str(path_root) + '/modules')
from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='bias analysis', \
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')

    return parser.parse_args()

def main():
    args = parse_args()

    # Build matrix of branch lengths between taxa in the tree
    df = build_branch_length_matrix(args.tree_file)

    





if __name__ == '__main__':
    main()
