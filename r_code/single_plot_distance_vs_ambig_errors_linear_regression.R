library(ggplot2)
library(ggpubr)

setwd("/home/vortacs/ref_bias_work/fixed_alignments/hyp_a_distance_to_errors/error_diffs/summary_files")

compare_data <- read.csv("distances_added_combined_summary.csv", header = TRUE, sep=",")
print(colnames(compare_data))

attach(compare_data)

model <- lm(all_ambiguous_diffs~distance)

png(filename = "all_seqs_single_plot_ambig_errors_vs_distance_linear_regression.png", width = 600, height = 600)

plt = ggplot(data=compare_data, aes(x=distance, y=all_ambiguous_diffs)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=0.0025, label.y=100000)

plt = plt + labs(title = "Distance to Reference vs. Ambiguous Errors", y = "Ambiguous Errors", x = "Distance to Reference")
plt

dev.off()
