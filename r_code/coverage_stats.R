library(ggplot2)
library(ggpubr)

# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work/fixed_alignments/hyp_c_errors_at_coverage_levels/coverage_results_and_analyses")

coverage_rate_data <- read.csv("/home/vortacs/ref_bias_work/fixed_alignments/hyp_c_errors_at_coverage_levels/coverage_results_and_analyses/combined_coverage_diffs_and_total_counts.csv", header = TRUE, sep=",")
print(colnames(coverage_rate_data))

attach(coverage_rate_data)

expanded_data = coverage_rate_data

expanded_data$rate = coverage_rate_data$diffs_count / coverage_rate_data$total_count

print(colnames(expanded_data))

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#coverage_rate_data$error_rates = (coverage_rate_data$count / coverage_rate_data$total_unambig_identical)
#print(coverage_rate_data$error_rates)

#subset = coverage_rate_data[1:271,,drop=F]
subset = expanded_data[2:315,,drop=F]
print(subset)

#plot(
#  x=expanded_data$coverage, 
#  y=expanded_data$rate, 
#  xlab="Coverage", 
#  ylab="Error Rate", 
#  main="Coverage vs. Error Rate",
#  col="grey52",
#  pch=19)

model = lm(subset$rate~subset$coverage)
print(summary(model))

#png(filename = "all_seqs_single_plot_unambig_coverage_vs_error_rate_linear_regression.png", width = 600, height = 400)

plt = ggplot(data=subset, aes(x=coverage, y=rate)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=100)
#plt = ggplot(data=expanded_data, aes(x=coverage, y=rate)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=100)


#plt = plt + scale_y_log10()
#plt = plt + scale_y_sqrt()
#breaks = c(0, 0.001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.09)

plt = plt + labs(title = "Coverage vs. Error Rate", y = "Unambiguous error Rate", x = "Coverage")
plt

#dev.off()

#sink("all_taxa_unambig_errors_to_distance_regression_results.txt")

write.csv(expanded_data, "/home/vortacs/ref_bias_work/fixed_alignments/hyp_c_errors_at_coverage_levels/coverage_results_and_analyses/unambig_coverage_vs_error_rate_table.csv")






#model <- lm(all_diffs~distance)
model <- lm(subset$error_rate~subset$coverage)

summary(model)

summary(model)$coefficients[1]

#str(summary(model))

print(lmp(model))


abline(model, col="blue", lwd=3)
