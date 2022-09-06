
# Analaysis of sequence comparison based on reference selection
setwd("/Users/jtfield/git_repos/reference_bias")

compare_data <- read.csv("all_diffs_summary_with_distances.csv", header = TRUE, sep=",")
print(colnames(compare_data))

attach(compare_data)
 
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


model <- lm(all_diffs~distance)

p_val = lmp(model)

r_sqr = summary(model)$r.squared

print(r_sqr)

#summary(model)

#print(summary(model)$r.squared)

#print(summary(model)$adj.r.squared)

# Examining any base that doesnt match the true seq
#scatter.smooth(distance, all_diffs, main='Ambiguous errors vs. Distance to reference')

plot(
  x=compare_data$distance, 
  y=compare_data$all_diffs, 
  xlab="Distance", 
  ylab="Errors", 
  main="Ambiguous errors vs. Distance to reference",
  col="grey52",
  pch=19)

#sink("all_taxa_unambig_errors_to_distance_regression_results.txt")

#model <- lm(all_diffs~distance)

summary(model)

summary(model)$coefficients[1]

str(summary(model))

print(lmp(model))


abline(model, col="blue", lwd=3)

#sink()






###### CHI SQUARED TESTING STUFF
#compare_data$expected_ref_biased_diffs = (compare_data$diffs_to_true_matching_ref * 0.25)

#k = summary(table(compare_data$distance, compare_data$diffs_to_true_matching_ref))

#k2 = chisq.test(compare_data$distance, compare_data$expected_ref_biased_diffs)

#k2
#k2$observed
#k2$expected

#compare_data$expected_ref_biased_diffs

#compare_data$unbiased_diffs = (compare_data$diffs_to_true - compare_data$diffs_to_true_matching_ref)

#compare_data$unbiased_diffs

#chisq.test(compare_data$unbiased_diffs, compare_data$diffs_to_true_matching_ref)

#chisq.test(compare_data$expected_ref_biased_diffs, compare_data$diffs_to_true_matching_ref)
