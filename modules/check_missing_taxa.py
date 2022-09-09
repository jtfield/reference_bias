#!/usr/bin/env python3
"""example usage:

"""

import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(prog='tree_combine_to_nexus.py', \
        description='Run after combine_cov_vcfs.py to replace the sra numbers in the file names with the sample IDs.')
    parser.add_argument('--taxa_list_file', help='taxa list file.')
    #parser.add_argument('--output_nexus', default='combined_phylos.nexus', help='Nexus file of multiple phylogenies.')
    # parser.add_argument('--output_coverage_dir', default='coverage_analysis_results', help='directory of updated unambiguous differences files that now include coverage data for each ref_base.')
    return parser.parse_args()

def main():
    args = parse_args()

    read_file = open(args.taxa_list_file, 'r').read()

    # strip_unnecessary = read_file.replace('combine_and_infer', '').replace('ep_dev_log.txt', '').replace('intermediate_files','').replace('RESULTS','')
    storage_dict = {}

    split_lists = read_file.split('###')
    for tax_list in split_lists:
        split_out_ref = tax_list.strip().split('\n', 1)
        if len(split_out_ref) > 1:
            ref_name = split_out_ref[0]
            seq_list = []
            present_taxa = split_out_ref[1].split('\n')
            for seq_name in present_taxa:
                stripped_name = seq_name.rsplit('/', 1)
                align_file = stripped_name[1]
                seq_list.append(align_file)
            storage_dict[ref_name] = seq_list

    counter = 0
    for key1, value1 in storage_dict.items():
        if len(value1) == 54:
            if counter == 0:
                for key2, value2 in storage_dict.items():
                    if key1 != key2:
                        # print("###")
                        # print("REF: ", key2)
                        diffs = set(value1).difference(value2)
                        if len(diffs) != 0:
                            print("###")
                            print("REF: ", key2)
                            print(diffs)
                        #print(set(value2).difference(value1))
        counter+=1



def common(a, b):
    c = [value for value in a if value in b]
    return c






if __name__ == '__main__':
    main()
