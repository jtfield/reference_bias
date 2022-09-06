#library(stargazer)

# Analaysis of sequence comparison based on reference selection
setwd("/Users/jtfield/git_repos/reference_bias/distance_summaries")

#compare_data <- read.csv("updated_summary_diffs.csv", header = TRUE, sep=",")
#colnames(compare_data)

#attach(compare_data)

# Examining any base that doesnt match the true seq
#scatter.smooth(distance, diffs_to_true, main='Nucleotides differing from True seq vs. Distance to reference')

#model <- lm(distance~diffs_to_true)

#summary(model)

data_files = lapply(Sys.glob("dists_added_diff_summary_75_query_*.csv"), read.csv)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


reg_results = matrix(ncol=4, nrow=55)

reg_results[1,1] = "taxon_name"
reg_results[1,2] = "p-value"
reg_results[1,3] = "r-squared"
reg_results[1,4] = "adj_r-squared"

#sink("single_query_taxon_ambiguous_errors_regression_summary.txt")

for (i in 1:length(data_files)) {

  current_table = data_files[[i]]
  current_row = i + 1
  
  #scatter.smooth(current_table$distance, current_table$all_diffs, xlab="Phylogenetic distance", ylab="Errors", main='Nucleotides errors vs. Distance to reference', sub=current_table$ref_name)
  
  model <- lm(current_table$diffs_to_true~current_table$distance)

  #print(current_table$taxon_name[1])
  tax_name = format(current_table$taxon_name[1])
  print(tax_name)
  #print(summary(model))
  
  print(lmp(model))
  print(summary(model)$r.squared)
  print(summary(model)$adj.r.squared)
  
  #reg_results[current_row,1] = current_table$taxon_name[1]
  reg_results[current_row,1] = tax_name
  reg_results[current_row,2] = lmp(model)
  reg_results[current_row,3] = summary(model)$r.squared
  reg_results[current_row,4] = summary(model)$adj.r.squared
  
  
  
  print("\n\n\n")

  
  
  
  #rm_underscores = gsub('_', ' ', current_table$taxon_name)
  #stargazer(model, type="text", title=rm_underscores)
  #reg_output = capture.output(stargazer(model, type="text", title=rm_underscores))
  #reg_results[[i]] = reg_output
  
}

print(reg_results)

stat_results = data.frame(reg_results)

write.csv(stat_results, "/Users/jtfield/git_repos/reference_bias/distance_summaries/query_plots/single_taxon_query_unambiguous_errors_distance_regression_results.csv", row.names = FALSE)

#sink()

#print(reg_results[1])
#lapply(reg_results, write, "single_taxon_query_regression_results.txt", append=TRUE, ncolumns=500)

#dev.off()


