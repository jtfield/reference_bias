# Analaysis of sequence comparison based on reference selection
setwd("/Users/jtfield/git_repos/reference_bias")

coverage_rate_data <- read.csv("/Users/jtfield/git_repos/full_counts_combined_unambig_coverage_counts.csv", header = TRUE, sep=",")
print(colnames(coverage_rate_data))

attach(coverage_rate_data)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#coverage_rate_data$error_rates = (coverage_rate_data$count / coverage_rate_data$total_unambig_identical)
#print(coverage_rate_data$error_rates)

subset = coverage_rate_data[1:271,,drop=F]
#print(subset)

plot(
  x=subset$coverage, 
  y=subset$error_rate, 
  xlab="Coverage", 
  ylab="Error Rate", 
  main="Coverage vs. Error Rate",
  col="grey52",
  pch=19)

#sink("all_taxa_unambig_errors_to_distance_regression_results.txt")

#model <- lm(all_diffs~distance)
model <- lm(subset$error_rate~subset$coverage)

summary(model)

summary(model)$coefficients[1]

#str(summary(model))

print(lmp(model))


abline(model, col="blue", lwd=3)
