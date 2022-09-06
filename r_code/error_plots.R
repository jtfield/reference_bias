
# Analaysis of sequence comparison based on reference selection
setwd("/Users/jtfield/git_repos/reference_bias/distance_summaries")

data_files = lapply(Sys.glob("dists_added_diff_summary_75_query_*.csv"), read.csv)

#scatter.smooth(distance, diffs_to_true, main='Nucleotides differing from True seq vs. Distance to reference')

destination = "/Users/jtfield/git_repos/reference_bias/distance_summaries/query_plots/query_only_ambiguous_errors_scatter_plots.pdf"

pdf(file=destination)

par(mfrow = c(2,2))

for (i in 1:length(data_files)) {
  
  current_table = data_files[[i]]
  
  scatter.smooth(current_table$distance, current_table$all_diffs, xlab="Phylogenetic distance", ylab="Errors", main='Nucleotides errors vs. Distance to reference', sub=current_table$taxon_name)
  
  
}
dev.off()


