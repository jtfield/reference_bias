#!/usr/bin/env python3
"""example usage:

"""

import os
import argparse
import dendropy

def parse_args():
    parser = argparse.ArgumentParser(prog='tree_combine_to_nexus.py', \
        description='Run after combine_cov_vcfs.py to replace the sra numbers in the file names with the sample IDs.')
    parser.add_argument('--tree_dir', help='Directory of multiple phylogeny files.')
    parser.add_argument('--output_nexus', default='combined_phylos.nexus', help='Nexus file of multiple phylogenies.')
    # parser.add_argument('--output_coverage_dir', default='coverage_analysis_results', help='directory of updated unambiguous differences files that now include coverage data for each ref_base.')
    return parser.parse_args()

def main():
    args = parse_args()

    tree_list = dendropy.TreeList()
    file_list = os.listdir(args.tree_dir)

    for tree_file in file_list:
        # print(tree_file)
        tree_file = args.tree_dir + '/' + tree_file
        # print(tree_file)

        tree_list.read(path=tree_file, schema="newick", preserve_underscores=True)

    tree_list.write(path=args.output_nexus, schema="nexus")







if __name__ == '__main__':
    main()
