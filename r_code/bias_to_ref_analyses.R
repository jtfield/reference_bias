# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work/fixed_alignments/hyp_a_distance_to_errors/error_diffs/ref_summary_files")

compare_data <- read.csv("distances_added_combined_summary.csv", header = TRUE, sep=",")
print(colnames(compare_data))

attach(compare_data)

# Part 1 of Hypothesis B
unambig_errors_sum = sum(compare_data$unambiguous_diffs_to_true)
print(unambig_errors_sum)

unambig_errors_matching_ref_sum = sum(unambiguous_diffs_to_true_matching_ref)
print(unambig_errors_matching_ref_sum)

#ambig_errors_sum = sum(compare_data$all_diffs)
ambig_errors_sum= sum(compare_data$all_ambiguous_diffs)
print(ambig_errors_sum)

#total_bases = 2328289 * 54
#print(total_bases)

#total_bases_compared = 2916 * 2328289
#print(total_bases_compared)

prop.test(x = unambig_errors_matching_ref_sum, n = unambig_errors_sum, alternative = 'two.sided', p = 0.25)

# Part 2 of Hypothesis B
#errors_subset_df = compare_data[ , c('diffs_to_true_matching_ref', 'diffs_to_true')]
errors_subset_df = compare_data[ , c('unambiguous_diffs_to_true_matching_ref', 'unambiguous_diffs_to_true')]

errors_subset_df$unambiguous_diffs_to_true_not_matching_ref = (errors_subset_df$unambiguous_diffs_to_true - errors_subset_df$unambiguous_diffs_to_true_matching_ref)

print(errors_subset_df)

ref_errors_df = errors_subset_df[ , c('unambiguous_diffs_to_true_matching_ref', 'unambiguous_diffs_to_true_not_matching_ref')]

sum_ref_errors = sum(ref_errors_df$unambiguous_diffs_to_true_matching_ref)
avg_ref_errors = mean(ref_errors_df$unambiguous_diffs_to_true_matching_ref)
min_ref_errors = min(ref_errors_df$unambiguous_diffs_to_true_matching_ref)
max_ref_errors = max(ref_errors_df$unambiguous_diffs_to_true_matching_ref)

print(sum_ref_errors)
print(avg_ref_errors)
print(min_ref_errors)
print(max_ref_errors)

sum_non_ref_errors = sum(ref_errors_df$unambiguous_diffs_to_true_not_matching_ref)
avg_non_ref_errors = mean(ref_errors_df$unambiguous_diffs_to_true_not_matching_ref)
min_non_ref_errors = min(ref_errors_df$unambiguous_diffs_to_true_not_matching_ref)
max_non_ref_errors = max(ref_errors_df$unambiguous_diffs_to_true_not_matching_ref)
print(sum_non_ref_errors)
print(avg_non_ref_errors)
print(min_non_ref_errors)
print(max_non_ref_errors)

#ref_errors_density = hist(ref_errors_df$unambiguous_diffs_to_true_matching_ref)
#plot(ref_errors_density)

#non_ref_errors_desnity = hist(ref_errors_df$unambiguous_diffs_to_true_not_matching_ref)
#plot(non_ref_errors_desnity)

#plot(ref_errors_density, col = rgb(0,0,1,1/4), xlim = c(0,4000), ylim = c(0,2000))
#plot(non_ref_errors_desnity, col = rgb(1,0,0,1/4), xlim = c(0,4000), ylim = c(0,2000), add = T)

chisq.test(ref_errors_df)


