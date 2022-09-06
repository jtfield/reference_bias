library(ggplot2)
library(ggpubr)

# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work/manuscript/hyp_a_errors_vs_distance/2022_08_23_unambig_summaries")

data_files = lapply(Sys.glob("dists_added_diff_summary_75_query_*.csv"), read.csv)

destination = "/home/vortacs/ref_bias_work/manuscript/hyp_a_errors_vs_distance/updated_individual_taxa_plots/query_only_ambiguous_errors_scatter_plots.pdf"

pdf(file=destination)

par(mfrow = c(2,2))

for (i in 1:length(data_files)) {
  
  current_table = data_files[[i]]

  model <- lm(current_table$unambiguous_diffs_to_true~current_table$distance)

  plt = ggplot(data=current_table, aes(x=distance, y=unambiguous_diffs_to_true)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=0.0025)

  plt = plt + labs(title = paste("Distance to Reference vs. Unambiguous Errors", current_table$ref_taxon[1]), y = "Unambiguous Errors", x = "Distance to Reference")
  plt
  
}
dev.off()
