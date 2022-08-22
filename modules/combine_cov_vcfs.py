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
from os.path import exists

def parse_args():
    parser = argparse.ArgumentParser(prog='multi_diff_counter.py', \
        description='Compare original sequences to ref_based sequences. One alignment of original sequences and one of ref based sequences.')
    parser.add_argument('--multi_ep_dir', help='Directory of multiple Extensiphy outputs.')
    parser.add_argument('--output_dir', default='coverage_vcf_files', help='Directory of relabelled coverage vcf files.')
    return parser.parse_args()

def main():
    args = parse_args()

    list_of_alignment_dirs = os.listdir(args.multi_ep_dir)
    print(list_of_alignment_dirs)

    for align_dir in list_of_alignment_dirs:
        print('###')
        print(align_dir)
        ref_taxon = align_dir
        path_to_align_dir = args.multi_ep_dir + '/' + align_dir
        ep_dir_container = os.listdir(path_to_align_dir)
        #print(ep_dir_container)

        ep_dir_contents = path_to_align_dir + '/ep_output'

        ep_run_outputs = os.listdir(ep_dir_contents)

        for taxon_dir in ep_run_outputs:
            if taxon_dir.endswith('output_dir'):
                taxon = taxon_dir.replace('output_dir','')
                #print(taxon)
                path_to_taxon_dir = ep_dir_contents + '/' + taxon_dir
                taxon_dir_contents = os.listdir(path_to_taxon_dir)
                #print(taxon_dir)
                #print(taxon_dir_contents)

                path_to_vcf = path_to_taxon_dir + '/dupes_removed_best_cns.vcf'
                exists_check = exists(path_to_vcf)

                if exists_check:

                    new_ref_and_taxon_vcf_name = 'ref_' + ref_taxon + '_query_' + taxon + '.vcf'
                    path_to_new_vcf = args.output_dir + '/' + new_ref_and_taxon_vcf_name
                    #print(path_to_new_vcf)

                    shutil.copyfile(path_to_vcf, path_to_new_vcf)

                else:
                    print("MISSING TAXON VCF FILE: {t}, ref {r}".format(t=taxon_dir, r=align_dir))







if __name__ == '__main__':
    main()
