#! /usr/bin/python3

import os
import argparse
import subprocess
from modules.alignment_splitter import *


def parse_args():
    parser = argparse.ArgumentParser(prog='Multi_EP', \
        description='Split an input alignment into separate sequences. Run EP using the same input data on each sequence.')
    parser.add_argument('--align_file', default=False, help='input alignment option.')
    parser.add_argument('--reads_dir', default=False, help='directory of reads.')
    parser.add_argument('--out_dir', default='multi_ep_output', help='Directory of multi_ep.py outputs.')
    parser.add_argument('--ep_path', default=False)
    parser.add_argument('--ep_runs', default=2)
    parser.add_argument('--cores_per_run', default=2)
    parser.add_argument('--tail_1', default="_1.fastq")
    parser.add_argument('--tail_2', default="_2.fastq")

    return parser.parse_args()

def main():
    args = parse_args()

    # Make output dir
    os.mkdir(args.out_dir)

    # Move into output dir
    os.chdir(args.out_dir)

    # Get absolute path for the output dir
    absolute_output_dir_path = os.path.abspath(os.getcwd())

    # Split alignment into separate sequences
    split_alignment(args.align_file, absolute_output_dir_path)

    list_of_taxon_dirs = os.listdir(absolute_output_dir_path)

    for taxon_dir in list_of_taxon_dirs:

        os.chdir(taxon_dir)

        current_abs_path = os.path.abspath(os.getcwd())

        sequence = os.listdir(current_abs_path)[0]

        # path_to_align = current_abs_path + "/ep_output/RESULTS/extended.aln"
        path_to_ref_removed_align = current_abs_path + "/" + taxon_dir + "_removed.aln"

        # Run EP using each sequence as an input
        subprocess.run([args.ep_path + "/extensiphy.sh", "-a", sequence, "-d", args.reads_dir, "-i", "CLEAN", "-p", str(args.ep_runs) ,"-c", str(args.cores_per_run), "-1", args.tail_1, "-2", args.tail_2, "-o", current_abs_path + "/ep_output"])

        # Read alignment and make a copy of the alignment with the reference sequence removed
        print("Building the extended EP alignment without the reference sequence.")

        # align = open(path_to_align, 'r').read()
        #
        # split_align = align.split(">")
        #
        # output = open(path_to_ref_removed_align, 'w')
        #
        # for chunk in split_align:
        #     if len(chunk) > 0:
        #         if taxon_dir not in chunk:
        #             output.write(">")
        #             output.write(chunk)
        #             output.write("\n")
        #
        # output.close()

        path_to_new_seqs = absolute_output_dir_path + "/" + taxon_dir + "/ep_output/combine_and_infer"

        list_of_new_seqs = os.listdir(path_to_new_seqs)

        output = open(path_to_ref_removed_align, 'w')

        for file in list_of_new_seqs:
            if file.endswith(".fas"):
                seq = open(path_to_new_seqs + "/" + file, 'r').read()
                output.write(seq)
                output.write("\n")


        os.chdir(absolute_output_dir_path)





if __name__ == '__main__':
    main()
