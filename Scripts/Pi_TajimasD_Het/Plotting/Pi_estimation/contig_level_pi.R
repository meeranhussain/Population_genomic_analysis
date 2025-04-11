# Load necessary libraries
library(ggpubr)
library(rstatix)
library(dplyr)

#load datasets of each population 
HAMpi <- read.table("HAM_pi_10kb.windowed.pi", header=T)
MANpi <- read.table("MAN_pi_10kb.windowed.pi", header=T)
LINpi <- read.table("LIN_pi_10kb.windowed.pi", header=T)
DUNpi <- read.table("DUN_pi_10kb.windowed.pi", header=T)
IREpi <- read.table("IRE_pi_10kb.windowed.pi", header=T)

# Add a population column to each
HAMpi$Population <- "HAM"
MANpi$Population <- "MAN"
LINpi$Population <- "LIN"
DUNpi$Population <- "DUN"
IREpi$Population <- "IRE"

# Merge into a single dataframe
data <- rbind(MANpi, HAMpi, LINpi, DUNpi, IREpi)

# Customize the color palette
color_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")

# Define the contigs of interest
selected_contigs <- c("contig_1", "contig_83", "contig_2", "contig_84", 
                      "contig_3", "contig_85", "contig_4", "contig_86")

# Filter data for selected contigs
data_filtered <- data %>%
  filter(CHROM %in% selected_contigs) 

# Set population order
data_filtered$Population <- factor(data_filtered$Population, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Perform ANOVA for each contig separately
anova_results <- data_filtered %>%
  group_by(CHROM) %>%
  summarise(p_value = signif(summary(aov(PI ~ Population, data = cur_data()))[[1]][["Pr(>F)"]][1], 4))

# Get maximum y value for annotation placement
max_y <- max(data_filtered$PI, na.rm = TRUE)

# Create Box + Violin Plot
plot_combined <- ggplot(data_filtered, aes(x = Population, y = PI, fill = Population)) +
  geom_violin(alpha = 0.4, width = 0.9) +  # Violin plot for density visualization
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Boxplot for summary stats
  facet_wrap(~CHROM, scales = "free_y") +  # Separate by contig
  labs(
    title = "Nucleotide Diversity Scaffold level",
    x = "Population",
    y = "Nucleotide Diversity (π)",
    fill = "Population"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.75, 0.1),  
    legend.justification = c(0.5, 0.5),  
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    legend.title = element_text(size = 10, face = "bold"),  # Legend title size
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),  
    legend.spacing.y = unit(0.8, "cm") 
  ) + 
  ylim(0, 0.00040)

# Add ANOVA p-values per contig
plot_combined <- plot_combined + 
  geom_text(data = anova_results, aes(x = 2.5, y = 0.00035, label = paste("pvalue:", p_value)), 
            inherit.aes = FALSE, size = 3, fontface = "italic")



# Create a mapping of contigs to scaffold names
contig_labels <- c(
  "contig_1" = "scaffold_1",
  "contig_2" = "scaffold_2",
  "contig_3" = "scaffold_3",
  "contig_4" = "scaffold_4",
  "contig_83" = "scaffold_83",
  "contig_84" = "scaffold_84",
  "contig_85" = "scaffold_85",
  "contig_86" = "scaffold_86"
)

# Update the facet labels using labeller
plot_combined <- plot_combined +
  facet_wrap(~contigs, scales = "free_y", labeller = labeller(contigs = contig_labels))




# Display the plot
print(plot_combined)

# Save the plot
ggsave("TajimasD_box_violin_contig_anova.png", plot = plot_combined, width = 12, height = 8, dpi = 300)


################################## t-Test ############################

# Filter data for selected contigs
data_filtered <- data %>%
  filter(contigs %in% selected_contigs)

# Set population order
data_filtered$name <- factor(data_filtered$name, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Perform **Wilcoxon pairwise test** for each contig
pairwise_pvals <- data_filtered %>%
  group_by(contigs) %>%
  wilcox_test(value ~ name, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "name",step.increase = 0.5)

pairwise_pvals <- pairwise_pvals[pairwise_pvals$p.adj.signif != "ns", ]


# Create Violin + Boxplot with `ggviolin`
violin_plot <- ggviolin(
  data_filtered,
  x = "name",
  y = "value",
  fill = "name",
  add = "boxplot",  # Adds boxplot inside violin
  palette = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"),
  add.params = list(fill = "white", width = 0.2)  # Boxplot customization
) +
  facet_wrap(~contigs, scales = "free_y") +  # Separate by contig
  labs(
    title = "Tajima's D Scaffold Level",
    x = "Population",
    y = "Tajima's D"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.75, 0.1),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.8, "cm")
  ) +
  stat_pvalue_manual(
    pairwise_pvals,
    label = "p.adj.signif",  # Show significance stars (***, **, *)
    tip.length = 0.01,
    bracket.shorten = 0.05
  )

# Create a mapping of contigs to scaffold names
contig_labels <- c(
  "contig_1" = "scaffold_1",
  "contig_2" = "scaffold_2",
  "contig_3" = "scaffold_3",
  "contig_4" = "scaffold_4",
  "contig_83" = "scaffold_83",
  "contig_84" = "scaffold_84",
  "contig_85" = "scaffold_85",
  "contig_86" = "scaffold_86"
)

# Update the facet labels using labeller
violin_plot <- violin_plot +
  facet_wrap(~contigs, scales = "free_y", labeller = labeller(contigs = contig_labels))

# Display the updated plot
print(violin_plot)

# Save the plot
ggsave("TajimasD_violin__contig_Wilcoxon.png", plot = violin_plot, width = 16, height = 10, dpi = 300)
