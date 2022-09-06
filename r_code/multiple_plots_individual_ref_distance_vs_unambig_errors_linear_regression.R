library(ggplot2)
library(ggpubr)

# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work/manuscript/hyp_a_errors_vs_distance/distance_summaries")

#data_files = lapply(Sys.glob("dists_added_diff_summary_75_query_*.csv"), read.csv)

destination = "/home/vortacs/ref_bias_work/manuscript/hyp_a_errors_vs_distance/updated_individual_taxa_plots/ref_only_ambiguous_errors_scatter_plots.pdf"

pdf(file=destination)

par(mfrow = c(2,2))

data_files = Sys.glob("dists_added_diff_summary_75_ref_*.csv")

for (i in 1:length(data_files)) {
  
  current_table = data_files[[i]]
  #print(current_table)
  dat = read.csv(current_table)
  #print(dat)
  
  #print(dat$taxon_name[1])
  table_title = paste("Distance to Reference vs. Unambiguous Errors: \n Ref", dat$taxon_name[1])

  model <- lm(dat$diffs_to_true~dat$distance)
  plt = ggplot(data=dat, aes(x=distance, y=diffs_to_true))
  plt = plt + geom_smooth(method="lm") 
  plt = plt + geom_point() 
  plt = plt + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=0.0025)
  plt = plt + labs(title = table_title, y = "Unambiguous Errors", x = "Distance to Reference")
  print(plt)
  
}
dev.off()
