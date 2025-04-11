# Load required libraries
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)  

# Load the data correctly
data <- read.csv("../../Heterozygosity/heterozygosity_results_vcftools.csv", header = TRUE)

# Check column names
print(colnames(data))  # Ensure 'Pop' exists

# Rename column if necessary
if (!"Pop" %in% colnames(data)) {
  stop("Error: Column 'Pop' not found in the dataset. Check column names with colnames(data).")
}

# Convert population names to factor for ordering (Re-do the factor conversion)
data <- data %>% mutate(Pop = factor(Pop, levels = c("MAN", "HAM", "LIN", "DUN", "IRE")))

# Perform Wilcoxon pairwise tests with Bonferroni correction
stat.test <- data %>% 
  wilcox_test(Obs_Het ~ Pop, p.adjust.method = "bonferroni")

# Adjust y positions for annotations
stat.test <- stat.test %>%
  add_xy_position(x = "Pop", step.increase = 0.05)

# Remove non-significant results
stat.test <- stat.test[stat.test$p.adj.signif != "ns", ]

# Generate the boxplot (Keeping ggboxplot)
p1 <- ggboxplot(data, x = "Pop", y = "Obs_Het", fill = "Pop",
           add.params = list(size = 1, alpha = 0.4)) +  # Add jitter for individual points
  #stat_summary(fun=median, geom='point', colour='black', size=3) +  # Add median points
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")) +
  labs(
    title = "",
    x = "Population",
    y = "Observed Heterozygosity",
    fill = "Population"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 16),  # Increased x-axis label size
    legend.position = "none"
  ) +
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",
    tip.length = 0.01, 
    bracket.shorten = 0.05,
    size = 8
  )

ggsave("heterozygosity_boxplot_pvalue_vcftools.png",p1, width = 14, height = 10, dpi = 300)


## Export pvalue stats
# Perform Wilcoxon pairwise tests with Bonferroni correction
stat.test <- data %>% 
  wilcox_test(Obs_Het ~ Pop, p.adjust.method = "bonferroni")

# Adjust y positions for annotations
stat.test <- stat.test %>%
  add_xy_position(x = "Pop", step.increase = 0.05)

# Convert list columns to character or numeric for export
stat.test <- stat.test %>%
  mutate(across(where(is.list), ~sapply(., paste, collapse = ", ")))

## export stats 
write.table(stat.test, file = "het_pvalue_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Compute median observed and expected heterozygosity for each population
median_het <- data %>%
  group_by(Pop) %>%
  summarise(
    Median_Obs_Het = median(Obs_Het, na.rm = TRUE),
    Median_Exp_Het = median(Exp_Het, na.rm = TRUE)
  )

# Print the median values
print(median_het)

## export median 
write.table(stat.test, file = "het_median_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Perform ANOVA tests separately for Obs_Het and Exp_Het
anova_obs <- aov(Avg_Obs_Het ~ Pop_ID, data = data)
anova_exp <- aov(Avg_Exp_Het ~ Pop_ID, data = data)
anova_pvalue_obs <- signif(summary(anova_obs)[[1]][["Pr(>F)"]][1], 4)
anova_pvalue_exp <- signif(summary(anova_exp)[[1]][["Pr(>F)"]][1], 4)