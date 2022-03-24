#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
# import dendropy
# import pathlib
# path_root = pathlib.Path(__file__).parents[0]
# sys.path.append(str(path_root) + '/modules')
# from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='ref bias results figure', \
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--csv_file', default=False, help='input csv option.')

    return parser.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.csv_file)

    # ref_bases_to_br = df
    scatter_plot(df)




def scatter_plot(table):
    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_neither_ref_nor_orig', color='red')
    # plt.show()
    #
    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_orig_nucs_only', color='red')
    # plt.show()
    #
    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_orig_only', color='red')
    # plt.show()
    # #
    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_neither_nucs_only', color='red')
    # plt.show()

    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_ref_nucs_only', color='red')
    # plt.show()
    #
    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_ref_only', color='red')
    # plt.show()

    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_both', color='red')
    # plt.show()
    #
    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='matches_both_nucs_only', color='red')
    # plt.show()

    # table.plot(kind='scatter', x='patristic_distance_to_ref', y='dont_match_ref', color='red')
    # plt.show()

    table.plot(kind='scatter', x='patristic_distance_to_ref', y='ref_orig_identical_dont_match_either', color='red')
    plt.show()

    table.plot(kind='scatter', x='patristic_distance_to_ref', y='ref_orig_nonidentical_dont_match_either', color='red')
    plt.show()


if __name__ == '__main__':
    main()
