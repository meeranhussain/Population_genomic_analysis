
########################################### FROH Contig level ########################################

############### Contig level Froh ###########################
# Load Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

# Read PLINK ROH output
roh_contig <- read.table("01_parameter/ROH_results_up.hom", header = TRUE, stringsAsFactors = FALSE)
colnames(roh_contig)[2] <- "IID"  # Rename Individual column


# Customize the color palette
color_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")

# Select relevant columns and convert KB to base pairs (bp)
roh_data <- roh_contig %>%
  select(IID, CHR, KB) %>%  # Select Individual ID, Contig (CHR), ROH Length (KB)
  rename(Individual = IID, Contig = CHR, ROH_Length_kb = KB) %>%
  mutate(ROH_Length_bp = ROH_Length_kb * 1000)  # Convert KB to base pairs

# Convert Contig to Factor for Proper Ordering
roh_data$Contig <- factor(roh_data$Contig, levels = unique(roh_data$Contig))

# Assign Populations Based on Sample Naming (Modify as Needed)
roh_data$Population <- case_when(
  grepl("^DUN", roh_data$Individual) ~ "DUN",
  grepl("^HAM", roh_data$Individual) ~ "HAM",
  grepl("^IRE", roh_data$Individual) ~ "IRE",
  grepl("^LIN", roh_data$Individual) ~ "LIN",
  grepl("^MAN", roh_data$Individual) ~ "MAN",
  TRUE ~ "Other"
)

# Convert Population to Factor with Specific Order
roh_data$Population <- factor(roh_data$Population, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Contig Sizes in Base Pairs
contig_sizes <- data.frame(
  Contig = c("contig_1", "contig_83", "contig_2", "contig_84", "contig_3", 
             "contig_85", "contig_4", "contig_86"),
  Contig_Size_bp = c(29470290, 28376149, 23977460, 11926033, 
                     10466058, 9529246, 7912293, 7312752)
)

# Merge contig sizes with ROH data
roh_data <- roh_data %>%
  left_join(contig_sizes, by = "Contig")

# Calculate Froh: Proportion of Each Contig Covered by ROH
roh_data <- roh_data %>%
  group_by(Individual, Contig, Population) %>%
  summarize(
    Total_ROH_bp = sum(ROH_Length_bp, na.rm = TRUE),  # Sum ROH lengths per contig
    Contig_Size_bp = unique(Contig_Size_bp)  # Ensure correct contig size
  ) %>%
  mutate(Froh = Total_ROH_bp / Contig_Size_bp)  # Compute FROH


############ Factor contigs based on size #######################
#roh_data$Contig <- factor(roh_data$Contig, 
#                          levels = c("contig_1", "contig_83", "contig_2", "contig_84", 
#                                     "contig_3", "contig_85", "contig_4", "contig_86"))



################################## t-Test ############################


# Perform **Wilcoxon pairwise test** for each contig
pairwise_pvals <- roh_data %>%
  group_by(Contig) %>%
  wilcox_test(Froh ~ Population, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Population",step.increase = 0.1)

# Create Violin + Boxplot with `ggviolin`
violin_plot <- ggviolin(
  roh_data,
  x = "Population",
  y = "Froh",
  fill = "Population",
  add = "boxplot",  # Adds boxplot inside violin
  palette = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"),
  add.params = list(fill = "white", width = 0.2)  # Boxplot customization
) +
  facet_wrap(~Contig, scales = "free_y") +  # Separate by contig
  labs(
    title = "F(ROH) Across Scaffolds",
    x = "Population",
    y = expression(F[ROH])
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
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
  facet_wrap(~Contig, scales = "free_y", labeller = labeller(Contig = contig_labels))

# Display the updated plot
print(violin_plot)



# Save the Plot
ggsave("01_parameter/Froh_by_Contig__p_value_voilin.png", violin_plot, width = 14, height = 10, dpi = 300)



