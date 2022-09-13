setwd("/home/vortacs/ref_bias_work/fixed_alignments/hyp_a_distance_to_errors/error_diffs/ref_summary_files")

compare_data <- read.csv("distances_added_combined_summary.csv", header = TRUE, sep=",")
print(colnames(compare_data))

attach(compare_data)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

model <- lm(unambiguous_diffs_to_true~distance)

print(summary(model))

print(lmp(model))
print(summary(model)$r.squared)
print(summary(model)$adj.r.squared)

model2 = lm(all_ambiguous_diffs~distance)

print(summary(model2))
