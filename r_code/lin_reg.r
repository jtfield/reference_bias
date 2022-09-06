#! /usr/bin/R

data = read.csv("/home/vortacs/ref_bias_work/reformat_diffs_testing/summary_files/updated_summary_diffs.csv", header = TRUE, sep=",")


# unambig_diffs = data[5]

#unambig_diffs_match_ref = data[6]

#dists = data[7]

#data["ref_names"]

colnames(data)

attach(data)

plot(data["diffs_to_true"])
