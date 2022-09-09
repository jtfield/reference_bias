#!/usr/bin/env python3

import os
import argparse
import dendropy
import pathlib
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(prog='remove taxa from tree', \
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')

    return parser.parse_args()

def main():
    args = parse_args()

    dp_tree = dendropy.Tree.get(path=args.tree_file, schema='newick', preserve_underscores=True)

    get_distances_range(dp_tree)


def get_distances_range(tree):
    # pritn("waffle")

    pdm = tree.phylogenetic_distance_matrix()

    taxon_list = []

    distance_list = []

    already_checked_taxon_pairs = []

    # Populate taxon list
    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        # print(dir(taxon1))
        # print(type(taxon1.label))
        # print(type(taxon1))
        str_taxon1 = str(taxon1.label)
        # print(type(str_taxon))
        taxon_list.append(str_taxon1)

    print(taxon_list)

    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        # print(taxon1)
        # print(idx1)
        str_taxon1 = str(taxon1.label)
        for taxon2 in tree.taxon_namespace:
            str_taxon2 = str(taxon2.label)
            # mrca = pdm.mrca(taxon1, taxon2)
            weighted_patristic_distance = pdm.patristic_distance(taxon1, taxon2)
            # distance_list.append(weighted_patristic_distance)

            tax_pair = (str_taxon1, str_taxon2)
            reverse_pair = (str_taxon2, str_taxon1)
            # print(tax_pair)

            if tax_pair not in already_checked_taxon_pairs or reverse_pair not in already_checked_taxon_pairs:
                # print("new pair")
                # print(tax_pair)
                already_checked_taxon_pairs.append(tax_pair)
                already_checked_taxon_pairs.append(reverse_pair)
                distance_list.append(weighted_patristic_distance)
            # else:
                # print("already found pair")

    #         br_df.at[str_taxon1, str_taxon2] = weighted_patristic_distance

    distance_list.sort()
    # print(distance_list)

    print(len(distance_list))

    no_self_distances = []
    for item in distance_list:
        if item != 0.0:
            no_self_distances.append(item)
    print(min(no_self_distances))
    print(max(no_self_distances))

    plt.hist(no_self_distances)
    plt.xlabel('Distance')
    plt.ylabel('Pairwise Comparisons')

    plt.show()


if __name__ == '__main__':
    main()
