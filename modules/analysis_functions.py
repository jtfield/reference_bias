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
    'matches_both', \
    'matches_orig_only', \
    'matches_ref_only', \
    'matches_neither_ref_nor_orig', \
    'matches_both_nucs_only', \
    'matches_orig_nucs_only', \
    'matches_ref_nucs_only', \
    'matches_neither_nucs_only']

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

            ref_name_regex = '('+ ref_name + '\s.+)\s'

            ref_regex_compile = re.compile(ref_name_regex)

            read_new_align = open(dir_of_new_aligns + '/' + align, 'r').read()

            find_ref_seq = re.findall(ref_regex_compile, orig_align_contents)
            if find_ref_seq:

                split_ref_name_and_seq = find_ref_seq[0].split('\n', 1)

                ref_name_check = split_ref_name_and_seq[0]
                ref_seq = split_ref_name_and_seq[1]

                split_new_align = read_new_align.split('>')

                for chunk in split_new_align:
                    if len(chunk) > 0:

                        # this_seqs_list_of_values = []

                        split_new_name_and_seq = chunk.split('\n', 1)

                        new_name = split_new_name_and_seq[0]
                        new_seq = split_new_name_and_seq[1]

                        orig_name_and_seq_regex = '('+ new_name + '\s.+)\s'
                        print(orig_name_and_seq_regex)

                        orig_name_and_seq_compile = re.compile(orig_name_and_seq_regex)

                        find_orig_seq = re.findall(orig_name_and_seq_compile, orig_align_contents)
                        if find_orig_seq:

                            split_orig_name_and_seq = find_orig_seq[0].split('\n', 1)

                            orig_name_check = split_orig_name_and_seq[0]
                            orig_seq = split_orig_name_and_seq[1]

                            assert new_name == orig_name_check

                            new_name = new_name.replace('_', ' ')
                            ref_name_check = ref_name_check.replace('_', ' ')
                            try:
                                br_to_ref = br_df.at[new_name, ref_name_check]
                                # this_seqs_list_of_values.append(new_name)
                                # this_seqs_list_of_values.append(ref_name_check)
                                # this_seqs_list_of_values.append(br_to_ref)

                                current_seq_check = make_comparison(new_seq, orig_seq, ref_seq)

                                # prior_to_df_list.append(current_seq_check)

                                restructure_counts(prior_to_df_list, current_seq_check, new_name, ref_name_check, br_to_ref)


                            except KeyError:
                                print("combination not found: ", new_name + ' ' + ref_name_check)

    return prior_to_df_list


                            # assert br_to_ref != ''
                            # this_seqs_list_of_values.append(new_name)
                            # this_seqs_list_of_values.append(ref_name_check)
                            # this_seqs_list_of_values.append(br_to_ref)
                            #
                            # current_seq_check = make_comparison(new_seq, orig_seq, ref_seq, this_seqs_list_of_values)


def make_comparison(new_seq, orig_seq, ref_seq):
    """
    Make comparison between the sequences and return counts of the results
    """

    output = []

    standard_nucleotides = ['A', 'C', 'G', 'T']
    gaps = ['-']
    degen_nucleotides = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']

    # matches_orig = 0
    # matches_ref = 0
    matches_both = 0
    matches_neither_ref_nor_orig = 0
    matches_only_orig = 0
    matches_only_ref = 0
    matches_both_nucs_only = 0
    matches_orig_nucs_only = 0
    matches_ref_nucs_only = 0
    matches_neither_nucs_only = 0

    new_seq = new_seq.strip('\n')
    orig_seq = orig_seq.strip('\n')
    ref_seq = ref_seq.strip('\n')

    # print(len(new_seq))
    # print(len(orig_seq))
    # print(len(ref_seq))
    #
    # print(new_seq[:20], new_seq[-20:])
    # print(orig_seq[:20], orig_seq[-20:])
    # print(ref_seq[:20], ref_seq[-20:])

    for num, nuc in enumerate(new_seq):

        nuc = nuc.upper()
        orig_nuc = orig_seq[num].upper()
        ref_nuc = ref_seq[num].upper()

        # print(nuc, orig_nuc, ref_nuc)

        if nuc == orig_nuc and nuc == ref_nuc:
            matches_both+=1
            # matches_orig+=1
            # matches_ref+=1

            if nuc in standard_nucleotides:
                matches_both_nucs_only+=1

        elif nuc == orig_nuc and nuc != ref_nuc:
            matches_only_orig+=1

            if nuc in standard_nucleotides:
                matches_orig_nucs_only+=1

        elif nuc == ref_nuc and nuc != orig_nuc:
            matches_only_ref+=1

            if nuc in standard_nucleotides:
                matches_ref_nucs_only+=1

        elif nuc != orig_nuc and nuc != ref_nuc:
            matches_neither_ref_nor_orig+=1

            if nuc in standard_nucleotides:
                matches_neither_nucs_only+=1

        else:
            print(nuc)
            print("found something unexpected")

    output.append(matches_both)
    output.append(matches_only_orig)
    output.append(matches_only_ref)
    output.append(matches_neither_ref_nor_orig)
    output.append(matches_both_nucs_only)
    output.append(matches_orig_nucs_only)
    output.append(matches_ref_nucs_only)
    output.append(matches_neither_nucs_only)

    # print(output)

    return output


def restructure_counts(primary_list, individual_counts, taxon_name, ref_name, br):
    """
    Restructure the individual sequence counts to include taxon and ref name and branch lengths
    """

    current_list = []

    current_list.append(taxon_name)
    current_list.append(ref_name)
    current_list.append(br)

    for item in individual_counts:
        current_list.append(item)

    print(current_list)

    primary_list.append(current_list)
