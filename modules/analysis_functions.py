#!/usr/bin/env python3

import os
import argparse
import pathlib
import pandas as pd
import subprocess
import datetime
import dateutil
import re
import dendropy



def build_branch_length_matrix(phylogeny):

    taxa = dendropy.TaxonNamespace()

    tree = dendropy.Tree.get(path=phylogeny, schema="newick", taxon_namespace=taxa)

    pdm = tree.phylogenetic_distance_matrix()

    taxon_list = []

    # Populate taxon list
    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        # print(dir(taxon1))
        # print(type(taxon1.label))
        # print(type(taxon1))
        str_taxon1 = str(taxon1.label)
        # print(type(str_taxon))
        taxon_list.append(str_taxon1)


    br_df = build_all_to_all_table(taxon_list)
    # print(br_df)

    # print(dir(taxa))
    # print(taxa.all_taxa_bitmask)

    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        # print(taxon1)
        # print(idx1)
        str_taxon1 = str(taxon1.label)
        for taxon2 in tree.taxon_namespace:
            str_taxon2 = str(taxon2.label)
            # mrca = pdm.mrca(taxon1, taxon2)
            weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)

            br_df.at[str_taxon1, str_taxon2] = weighted_patristic_distance

            # print(str_taxon1)
            # print(str_taxon2)
            # # print(mrca)
            # print(weighted_patristic_distance)
            # print("#####")

    return br_df

def build_all_to_all_table(taxon_list):
    """Input a list of taxon names and build an empty table matching each taxon
        pairwise"""

    df = pd.DataFrame(columns=taxon_list, index=taxon_list)

    return df

def build_basic_comparison_df(data, file_name):
    """
    Builds a table with columns matching:
    """

    columns = ['target_taxon', \
    'reference', \
    'patristic_distance_to_ref', \
    'total_identical_nucleotides', \
    'proportion_ident_nucleotides', \
    'identical_to_original_and_ref_seq', \
    'proportion_ident_to_orig_and_ref', \
    'total_gaps_and_degens', \
    'proportion_gaps_and_degens']

    # # output.append("ident_to_origin_seq")
    # output.append(total_ident)
    # output.append(ident_nucs)
    # output.append(ident_gaps)
    # output.append(ident_degens)
    # # output.append(non_ident_bases)
    #
    # # output.append("ident_to_ref_seq")
    # output.append(total_ref_ident)
    # output.append(ident_to_ref_nucs)
    # output.append(ident_to_ref_gaps)
    # output.append(ident_to_ref_degens)
    #
    # output.append("not_ident_to_any")
    # output.append(not_ident_to_any)

    df = pd.DataFrame(data, columns=columns)

    df.to_csv(file_name)

    return df


def prepare_seq_comparison(orig_align, dir_of_new_aligns, suffix, br_df):
    """
    Make comparisons between the original alignment/sequences and the alignments/sequences based on each reference"
    """

    is_orig_align = os.path.exists(orig_align)

    list_of_aligns = os.listdir(dir_of_new_aligns)

    list_of_taxa = [file_n.replace(suffix,'') for file_n in list_of_aligns]

    prior_to_df_list = []

    # df = build_basic_comparison_df()

    # Just a check to make sure you input an actual file
    if is_orig_align:

        # print(orig_align)
        orig_align_contents = open(orig_align, 'r').read()

        for align in list_of_aligns:
            # print(align)

            ref_name = align.replace(suffix, '')

            ref_name_regex = ref_name + '\s[\w | -]+'

            ref_regex_compile = re.compile(ref_name_regex)

            read_new_align = open(dir_of_new_aligns + '/' + align, 'r').read()

            find_ref_seq = re.findall(ref_regex_compile, orig_align_contents)
            if find_ref_seq:

                split_ref_name_and_seq = find_ref_seq[0].split('\n', 1)

                ref_name_check = split_ref_name_and_seq[0]
                ref_seq = split_ref_name_and_seq[1]
                # print("ref name check: ", ref_name_check)
                # print("found ref")
                # print(find_ref_seq)

                loop_over_aligns(prior_to_df_list, orig_align_contents, read_new_align, ref_name, ref_seq, br_df)

                # print(prior_to_df_list)

    # print(prior_to_df_list)
    return prior_to_df_list

def loop_over_aligns(prior_to_df_list, orig_align, new_align, ref_name, ref_seq, br_df):
    """
    Fills dataframe with simple comparison info
    """

    # print(ref_name)

    # split original align into chunks

    # prior_to_df_list = []

    split_orig_align = orig_align.split('>')

    # Split new alignment
    split_new_align = new_align.split('>')

    # Loop over the sequences in new align and find the correct seq in the original align
    for chunk_1 in split_new_align:
        if len(chunk_1) > 0:

            split_name_and_seq_new = chunk_1.split('\n', 1)
            name_1 = split_name_and_seq_new[0]
            seq_1 = split_name_and_seq_new[1]


            for chunk_2 in split_orig_align:
                if name_1 in chunk_2:
                    # print("found name: ", name_1)

                    split_name_and_seq_orig = chunk_2.split('\n', 1)

                    name_2 = split_name_and_seq_orig[0]
                    seq_2 = split_name_and_seq_orig[1]

                    analysis_results = make_comparisons(seq_2, seq_1, ref_seq, br_df, prior_to_df_list, name_1, ref_name)
                    # print(prior_to_df_list)


def make_comparisons(orig_seq, new_seq, ref_seq, br_df, list_of_lists, new_seq_name, ref_seq_name):
    """
    Function to make comparison between sequences
    """
    # print("waffle")

    standard_nucleotides = ['A', 'C', 'G', 'T']
    gaps = ['-']
    degen_nucleotides = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']

    output = []

    ident_nucs = 0
    ident_to_ref_seq = 0
    ident_gaps_and_degens = 0


    # total_ident = 0
    # ident_nucs = 0
    # ident_gaps = 0
    # ident_degens = 0
    # ident_gaps_and_degens = 0
    #
    #
    # total_ref_ident = 0
    # ident_to_ref_nucs = 0
    # ident_to_ref_gaps = 0
    # ident_to_ref_degens = 0
    #
    # not_ident_to_any = 0

    br_to_ref = ''
    new_seq_name = new_seq_name.replace('_', ' ')
    ref_seq_name = ref_seq_name.replace('_', ' ')

    len_orig = len(orig_seq)
    len_ref = len(ref_seq)
    len_new = len(new_seq)

    # print(new_seq_name)
    # print(ref_seq_name)
    try:
        br_to_ref = br_df.at[new_seq_name, ref_seq_name]
        # print(br_to_ref)
        # print(type(br_to_ref))
    except KeyError:
        print("combination not found: ", new_seq_name + ' ' + ref_seq_name)


    for num, orig_nuc in enumerate(orig_seq):
        if orig_nuc != None and new_seq[num] != None:
            # print(nuc)
            # print(new_seq[num])
            # print("*****")

            orig_nuc = orig_nuc.upper()
            new_nuc = new_seq[num].upper()

            if new_nuc in gaps or new_nuc in degen_nucleotides:
                ident_gaps_and_degens+=1


            if orig_nuc in standard_nucleotides and new_nuc in standard_nucleotides:
                if orig_nuc == new_nuc:
                    ident_nucs+=1

                    ref_nuc = ''

                    if num < len(ref_seq):
                        ref_nuc = ref_seq[num].upper()

                        if len(ref_nuc) > 0 and new_nuc == ref_nuc:
                            ident_to_ref_seq+=1

    calc_proportion_ident = ident_nucs / len_new
    calc_proportion_ident_to_ref = ident_to_ref_seq / len_new
    calc_proportion_degens = ident_gaps_and_degens / len_new

    output.append(new_seq_name)
    output.append(ref_seq_name)
    output.append(br_to_ref)

    output.append(ident_nucs)
    output.append(calc_proportion_ident)
    output.append(ident_to_ref_seq)
    output.append(calc_proportion_ident_to_ref)
    output.append(ident_gaps_and_degens)
    output.append(calc_proportion_degens)

    # output.append(new_seq_name)
    # output.append(ref_seq_name)
    # output.append(br_to_ref)
    #
    # output.append(total_ident)
    # output.append(ident_nucs)
    # output.append(ident_gaps)
    # output.append(ident_degens)
    #
    # output.append(total_ref_ident)
    # output.append(ident_to_ref_nucs)
    # output.append(ident_to_ref_gaps)
    # output.append(ident_to_ref_degens)
    #
    # output.append(not_ident_to_any)

    list_of_lists.append(output)

    print(output)

    return list_of_lists
