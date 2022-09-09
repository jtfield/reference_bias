library(ggplot2)
library(ggpubr)

# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work/fixed_alignments/hyp_a_distance_to_errors/error_diffs/query_summary_files/query_summaries_with_dists")

#data_files = lapply(Sys.glob("dists_added_diff_summary_75_query_*.csv"), read.csv)

destination = "query_only_ambiguous_errors_scatter_plots.pdf"

pdf(file=destination)

#par(mfrow = c(2,2))

data_files = Sys.glob("dists_added_diff_summary_75_query_*.csv")

for (i in 1:length(data_files)) {
  
  current_table = data_files[[i]]
  #print(current_table)
  dat = read.csv(current_table)
  #print(dat)
  
  #print(dat$taxon_name[1])
  table_title = paste("Distance to Reference vs. Unambiguous Errors: \n query", dat$taxon_name[1])
  
  model <- lm(dat$unambiguous_diffs_to_true~dat$distance)
  plt = ggplot(data=dat, aes(x=distance, y=unambiguous_diffs_to_true))
  plt = plt + geom_smooth(method="lm") 
  plt = plt + geom_point() 
  plt = plt + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=0.0025)
  plt = plt + labs(title = table_title, y = "Unambiguous Errors", x = "Distance to Reference")
  print(plt)

}
dev.off()
