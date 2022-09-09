library(ggplot2)
library(ggpubr)

setwd("/home/vortacs/ref_bias_work/fixed_alignments/hyp_a_distance_to_errors/error_diffs/summary_files")

compare_data <- read.csv("distances_added_combined_summary.csv", header = TRUE, sep=",")
print(colnames(compare_data))

attach(compare_data)

model <- lm(unambiguous_diffs_to_true~distance)

png(filename = "all_seqs_single_plot_unambig_errors_vs_distance_linear_regression.png", width = 600, height = 600)

plt = ggplot(data=compare_data, aes(x=distance, y=unambiguous_diffs_to_true)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=0.0025)

plt = plt + labs(title = "Distance to Reference vs. Unambiguous Errors", y = "Unambiguous Errors", x = "Distance to Reference")
plt

dev.off()
