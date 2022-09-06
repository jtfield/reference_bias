library(ggplot2)
library(ggpubr)

setwd("/home/vortacs/ref_bias_work/manuscript/hyp_a_errors_vs_distance")

compare_data <- read.csv("updated_summary_diffs.csv", header = TRUE, sep=",")
print(colnames(compare_data))

attach(compare_data)

model <- lm(unambiguous_diffs_to_true~distance)

plt = ggplot(data=compare_data, aes(x=distance, y=unambiguous_diffs_to_true)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=0.0025)

plt = plt + labs(title = "Distance to Reference vs. Unambiguous Errors", y = "Unambiguous Errors", x = "Distance to Reference")
plt