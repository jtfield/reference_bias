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
        taxon_list.append(taxon1)


    br_df = build_branch_length_table(taxon_list)
    # print(br_df)

    # print(dir(taxa))
    # print(taxa.all_taxa_bitmask)

    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        # print(taxon1)
        # print(idx1)
        for taxon2 in tree.taxon_namespace:
            mrca = pdm.mrca(taxon1, taxon2)
            weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)

            br_df.at[taxon1, taxon2] = weighted_patristic_distance

            # print(taxon1)
            # print(taxon2)
            # # print(mrca)
            # print(weighted_patristic_distance)
            # print("#####")

    return br_df

def build_branch_length_table(taxon_list):
    """Input a list of taxon names and build an empty table matching each taxon
        pairwise"""

    df = pd.DataFrame(columns=taxon_list, index=taxon_list)

    return df
