
# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work/birds/estandia")

compare_data <- read.csv("ref_bias_comparison.csv", header = TRUE, sep=",")
colnames(compare_data)

attach(compare_data)

# Examining any base that doesnt match the true seq
scatter.smooth(patristic_distance_to_ref, doesnt_match_orig, main='Hours studied vs. Exam Score')

boxplot(doesnt_match_orig)

model <- lm(doesnt_match_orig~patristic_distance_to_ref)

summary(model)

# Examining only called nucleotides that dont match the true seq
scatter.smooth(patristic_distance_to_ref, doesnt_match_orig_nucs_only, main='Hours studied vs. Exam Score')

boxplot(doesnt_match_orig_nucs_only)

model <- lm(doesnt_match_orig_nucs_only~patristic_distance_to_ref)

summary(model)


# Examining only nucleotides that match the true seq, no degen
scatter.smooth(patristic_distance_to_ref, matches_orig_nucs_only, main='Hours studied vs. Exam Score')

boxplot(matches_orig_nucs_only)

model <- lm(matches_orig_nucs_only~patristic_distance_to_ref)

summary(model)


# Examining bases that match the true seq
scatter.smooth(patristic_distance_to_ref, matches_orig, main='Hours studied vs. Exam Score')

boxplot(matches_orig)

model <- lm(matches_orig~patristic_distance_to_ref)

summary(model)
