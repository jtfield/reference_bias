library(ggplot2)
library(ggpubr)

# Analaysis of sequence comparison based on reference selection
setwd("/home/vortacs/ref_bias_work")

tree_length_data <- read.csv("/home/vortacs/ref_bias_work/missing_data_vs_ref_tree_length_change.csv", header = TRUE, sep=",")

print(colnames(tree_length_data))

attach(tree_length_data)

model = lm(tree_length_data$tree_length_difference~tree_length_data$missing_data_count)
print(summary(model))

plt = ggplot(data=tree_length_data, aes(x=missing_data_count, y=tree_length_difference)) + geom_smooth(method="lm") + geom_point() + stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x=300000)
plt = plt + labs(title = "Missing data vs. Ref based tree change", y = "ref based tree change from true", x = "missing data")
plt
