#!/usr/bin/env python3

import os
import argparse
import pathlib
import numpy as np
import pandas as pd
import subprocess
import datetime
# from multiprocessing import Pool, freeze_support
# import multiprocessing as mp

def find_specific_seqs(input_path, align_list, outdir, name_to_find):
    """
    loop over files, find the sequence and name you want.
    Add it to the new alignment file.
    """

    suffix = '_removed.aln'

    output_file = outdir + '/' + name_to_find + '_separated.aln'

    output = open(output_file,'w')

    for file in align_list:

        ref_name = file.strip(suffix)

        path_to_file = input_path + '/' + file

        read_current_file = open(path_to_file).read()

        split_file = read_current_file.split('>')

        for chunk in split_file:
            if name_to_find in chunk:

                split_chunk = chunk.split('\n', 1)

                name = split_chunk[0]
                seq = split_chunk[1]

                seq_label = '>' + name + '--' + ref_name

                # Write selected sequences to output file
                output.write(seq_label)
                output.write('\n')
                output.write(seq)
                output.write('\n')
