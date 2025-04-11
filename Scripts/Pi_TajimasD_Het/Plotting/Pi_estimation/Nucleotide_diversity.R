##Plotting Sliding windows analysis - pi
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(ggpubr)
library(rstatix)

#load datasets of each population 
HAMpi <- read.table("HAM_pi_100kb.windowed.pi", header=T)
MANpi <- read.table("MAN_pi_100kb.windowed.pi", header=T)
LINpi <- read.table("LIN_pi_100kb.windowed.pi", header=T)
DUNpi <- read.table("DUN_pi_100kb.windowed.pi", header=T)
IREpi <- read.table("IRE_pi_100kb.windowed.pi", header=T)


# Add a population column to each
HAMpi$Population <- "HAM"
MANpi$Population <- "MAN"
LINpi$Population <- "LIN"
DUNpi$Population <- "DUN"
IREpi$Population <- "IRE"

# Merge into a single dataframe
data <- rbind(MANpi, HAMpi, LINpi, DUNpi, IREpi)

#calculating median and 5th and 95th percentile
median(HAMpi$PI)
quantile(HAMpi$PI, probs = c(0.05,0.95))
summary(HAMpi$PI)


median(MANpi$PI)
quantile(MANpi$PI, probs = c(0.05,0.95))
summary(MANpi$PI)


median(LINpi$PI)
quantile(LINpi$PI, probs = c(0.05,0.95))
summary(LINpi$PI)


median(DUNpi$PI)
quantile(DUNpi$PI, probs = c(0.05,0.95))
summary (DUNpi$PI)

median(IREpi$PI)
quantile(IREpi$PI, probs = c(0.05,0.95))
summary(IREpi$PI)


###VIOLIN PLOTS

length_HAMpi <- length(HAMpi$PI)
length_MANpi <- length(MANpi$PI)
length_LINpi <- length(LINpi$PI)
length_DUNpi <- length(DUNpi$PI)
length_IREpi <- length(IREpi$PI)

# Customize the color palette
color_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")


# Set the custom order for the `pop` factor
data$Population <- factor(data$Population, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Calculate mean nucleotide diversity for each population
mean_rates <- data %>%
  group_by(Population) %>%
  summarize(mean_pi = mean(PI, na.rm = TRUE))

# Perform ANOVA test
anova_test <- aov(PI ~ Population, data = data)
anova_summary <- summary(anova_test)

# Extract p-value from the ANOVA test
anova_pvalue <- signif(anova_summary[[1]][["Pr(>F)"]][1], 4)  # Round p-value to 4 significant digits



bxp <- ggboxplot(
  data,  # Use filtered data
  x = "Population",
  y = "PI",
  fill = "Population",
  palette = color_palette,
  outlier.shape = NA
) +
  labs(
    title = "Nucleotide Diversity",
    x = "Population",
    y = "Nucleotide Diversity (π)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  ylim(0, 0.00015)  # Explicitly set y-axis limits


print(bxp)


# Add the exact ANOVA p-value annotation to the plot
bxp <- bxp + annotate(
  "text",
  x = 2.5,  # Adjust the position of the p-value annotation as needed
  y = 0.00015,  # Position slightly below the upper limit of the y-axis
  label = paste("p-value:", anova_pvalue),
  size = 5,
  fontface = "italic"
)

# Add mean nucleotide diversity with "µ = " format
bxp <- bxp + stat_summary(
  fun = mean,
  geom = "text",
  aes(label = paste0("µ = ", round(..y.., 5))),  # Use µ for mean, rounded to 5 decimals
  position = position_dodge(width = 0.75),
  size = 4,
  color = "black",
  vjust = -0.5  # Adjust vertical position of mean text
)

# Display the plot
print(bxp)

# Save the plot
ggsave("nucleotide_diversity_anova_plot.png", plot = bxp, width = 10, height = 8, dpi = 300)




############# Voilin plot #####################################
# Create the violin plot for Total ROH Length (Mb)
color_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")



# Perform pairwise t-tests
stat.test <- data %>% 
  wilcox_test(PI ~ Population, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test <- stat.test %>%
  add_xy_position(x = "Population", step.increase = 0.2)

# Filter out "ns" from stat.test
stat.test <- stat.test[stat.test$p.adj.signif != "ns", ]

###### Set p_value position
stat.test$y.position <- seq(0.00035, 0.00038, length.out = nrow(stat.test))

### filter just for plotting
data <- data %>%
  filter(PI >= 0 & PI <= 0.00040)

# Plot using ggviolin
p1 <- ggviolin(data, x = "Population", y = "PI", fill = "Population",
               add = "boxplot", add.params = list(fill = "white", width = 0.1, alpha = 0.4)) +
  stat_summary(fun=median, geom='point', colour='black', size=2) +  # Add median points
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")) +  # Custom colors
  theme_minimal() +
  theme(
    text = element_text(size = 22),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )  +
  labs( x = "Population",  y = "Nucleotide Diversity (π)", title = " ") +
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05,
    size = 8
  ) + 
  ylim(0, 0.00040)

ggsave("Nucleotide_diveristy_voilinplot_pvalue.png",p1, width = 14, height = 10, dpi = 300)



## Export pvalue stats
# Perform pairwise t-tests
stat.test <- data %>% 
  wilcox_test(PI ~ Population, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test <- stat.test %>%
  add_xy_position(x = "Population", step.increase = 0.2)

# Convert list columns to character or numeric for export
stat.test <- stat.test %>%
  mutate(across(where(is.list), ~sapply(., paste, collapse = ", ")))

## export stats 
write.table(stat.test, file = "100kb_Pi_pvalue_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)



