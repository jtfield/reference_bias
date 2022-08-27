#!/usr/bin/env python3
"""example usage:
python diff_counter.py seq1.fas seq2.fas
Both sequences much be in fasta
"""

import os
import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression


def parse_args():
    parser = argparse.ArgumentParser(prog='add_distances_to_diff_data.py', \
        description='Combine summary difference files and add distance metrics.')
    parser.add_argument('--data_csv', default=False, help='differing bases and distance to reference per sequence.')

    return parser.parse_args()

def main():
    args = parse_args()

    data = pd.read_csv(args.data_csv)

    print(data.columns)

    print(data['all_diffs'].max())

    print(data['all_diffs'].min())

    print(data['diffs_to_true'].max())

    print(data['diffs_to_true'].min())








    # # model = LinearRegression()
    #
    # diffs_to_true_data = np.array(data.diffs_to_true.tolist())
    #
    # dists_data = np.array(data.distance.tolist()).reshape((-1,1))
    #
    # model = LinearRegression().fit(dists_data, diffs_to_true_data)
    #
    # # print(help(model))
    #
    # r_sq = model.score(dists_data, diffs_to_true_data)
    # print(f"coefficient of determination: {r_sq}")
    #
    # print(model.coef_)
    #
    # print(model.intercept_)





if __name__ == '__main__':
    main()
