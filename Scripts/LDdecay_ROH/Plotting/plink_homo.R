# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

# Load PLINK ROH results
roh_data <- read.table("01_parameter/ROH_results.hom_up.indiv", header = TRUE)
colnames(roh_data)[2] <- "IID"  # Rename Individual column

# Ensure population names are assigned (Modify this based on your population groups)
roh_data$Group <- factor(roh_data$FID, 
                         levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Convert total length from base pairs to Megabases
roh_data$Total_ROH_Mb <- roh_data$KB / 1000  # Convert KB to Mb

# Calculate FROH (Proportion of Genome in ROH)
# Assuming the genome size is ~130Mb (Modify if necessary)
genome_size_mb <- 130
roh_data$FROH <- roh_data$Total_ROH_Mb / genome_size_mb

# Define color palette
color_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")


######################### Box plot #########################

# Create the boxplot
p1 <- ggboxplot(
  roh_data, 
  x = "Group", 
  y = "Total_ROH_Mb", 
  fill = "Group", 
  palette = color_palette,
  outlier.shape = NA  # Hide outliers for cleaner visualization
)

stat.test <- roh_data %>% 
  t_test(Total_ROH_Mb ~ Group, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test <- stat.test %>%
  add_xy_position(x = "Group", step.increase = 0.2)

# Filter out "ns" from stat.test
stat.test <- stat.test[stat.test$p.adj.signif != "ns", ]

###### Set p_value position
stat.test$y.position <- seq(20, 25, length.out = nrow(stat.test))

# Customize the plot with p-value annotations
final_plot_p1 <- p1 + 
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05,
    size = 8
  ) +
  labs(
    title = "Total Length of ROHs (Mb)",
    x = "Population",
    y = "Total ROH Length (Mb)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),  # Increased x-axis label size
    axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave("01_parameter/Total_length_ROH.png", final_plot_p1, width = 10, height = 5, dpi = 300)


# Create the boxplot
p2 <- ggboxplot(
  roh_data, 
  x = "Group", 
  y = "FROH", 
  fill = "Group", 
  palette = color_palette,
  outlier.shape = NA  # Hide outliers for cleaner visualization
)

stat.test_froh <- roh_data %>% 
  t_test(FROH ~ Group, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test_froh <- stat.test_froh %>%
  add_xy_position(x = "Group", step.increase = 0.2)

# Filter out "ns" from stat.test
stat.test_froh <- stat.test_froh[stat.test_froh$p.adj.signif != "ns", ]

###### Set p_value position
stat.test_froh$y.position <- seq(0.16, 0.2, length.out = nrow(stat.test_froh))

# Customize the plot with p-value annotations
final_plot_p2 <- p2 + 
  stat_pvalue_manual(
    stat.test_froh, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05,
    size = 8
  ) +
  labs(
    title = "",
    x = "Population",
    y = expression(F[ROH])
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),  # Increased x-axis label size
    axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave("01_parameter/FROH_plot.png", final_plot_p2, width = 10, height = 5, dpi = 300)


############## Total ROH plot ###################
# Load PLINK ROH results
roh_data <- read.table("01_parameter/ROH_results.hom_up.indiv", header = TRUE)
colnames(roh_data)[2] <- "IID"  # Rename Individual column

# Ensure population names are assigned (Modify this based on your population groups)
roh_data$Group <- factor(roh_data$FID, 
                         levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))



# Create the boxplot
p3 <- ggboxplot(
  roh_data, 
  x = "Group", 
  y = "NSEG", 
  fill = "Group", 
  palette = color_palette,
  outlier.shape = NA  # Hide outliers for cleaner visualization
)

######## P-value ##########
stat.test_NROH <- roh_data %>% 
  t_test(NSEG ~ Group, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test_NROH <- stat.test_NROH %>%
  add_xy_position(x = "Group", step.increase = 0.2)

# Filter out "ns" from stat.test
stat.test_NROH <- stat.test_NROH[stat.test_NROH$p.adj.signif != "ns", ]

###### Set p_value position
stat.test_NROH$y.position <- seq(45, 55, length.out = nrow(stat.test_NROH))

# Customize the plot with p-value annotations
final_plot_p3 <- p3 + 
  stat_pvalue_manual(
    stat.test_NROH, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05,
    size = 8
  ) +
  labs(
    title = "",
    x = "Population",
    y = "Number of ROHs"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),  # Increased x-axis label size
    axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave("01_parameter/N0.ROH_plot.png", final_plot_p3, width = 10, height = 5, dpi = 300)

################### ROH Length distribution ##############################

library(ggplot2)
library(dplyr)
library(tidyr)

# Sample data (Replace with your actual PLINK ROH results)
roh_data <- read.table("01_parameter/ROH_results_up.hom", header = TRUE, stringsAsFactors = FALSE)
colnames(roh_data)[2] <- "IID"  # Rename Individual column

# Convert KB to MB
roh_data <- roh_data %>%
  mutate(ROH_Length_Mb = KB / 1000)

# Define ROH length categories
roh_data <- roh_data %>%
  mutate(Category = case_when(
    ROH_Length_Mb >= 0.01 & ROH_Length_Mb < 0.5 ~ "0.01–0.5 Mb",
    ROH_Length_Mb >= 0.5 & ROH_Length_Mb < 1 ~ "0.5–1 Mb",
    ROH_Length_Mb >= 1 & ROH_Length_Mb < 3 ~ "1–3 Mb",
    ROH_Length_Mb >= 3 ~ ">3 Mb",
    TRUE ~ "Other"
  ))

# Summarize ROH count per category per contig
roh_summary <- roh_data %>%
  group_by(CHR, Category, FID) %>%
  summarise(Count = n(), .groups = "drop")

# Calculate ratio within each population and contig
roh_summary <- roh_summary %>%
  group_by(FID, CHR) %>%
  mutate(Ratio = Count / sum(Count)) 

# Reshape for for just summary
roh_summary_wide <- roh_summary %>%
  pivot_wider(names_from = Category, values_from = Ratio, values_fill = 0)

roh_summary$CHR <- factor(roh_summary$CHR, 
                          levels = c("contig_1", "contig_83", "contig_2", "contig_84", 
                                     "contig_3", "contig_85", "contig_4", "contig_86"))

# Define scaffold name mapping
contig_labels <- c(
  "contig_1" = "scaffold_1",
  "contig_2" = "scaffold_3",
  "contig_3" = "scaffold_5",
  "contig_4" = "scaffold_7",
  "contig_83" = "scaffold_2",
  "contig_84" = "scaffold_4",
  "contig_85" = "scaffold_6",
  "contig_86" = "scaffold_8"
)

# Ensure population names are assigned (Modify this based on your population groups)
roh_summary$Group <- factor(roh_summary$FID, 
                         levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Plot the data faceted by population (FID)
p6 <-  ggplot(roh_summary, aes(x = CHR, y = Ratio, color = Category, shape = Category, group = Category)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~Group, scales = "free_x", ncol = 3) +  # Adjust columns to create gap
  scale_x_discrete(labels = contig_labels) +  
  theme_minimal() +
  labs(
    title = "",
    x = "Scaffolds",
    y = "Ratio",
    color = "ROH Length",
    shape = "ROH Length"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 16),  # Increased x-axis label size
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.85, 0.1),  
    legend.justification = c(0.5, 0.5),  
    legend.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    legend.title = element_text(size = 8, face = "bold"),  # Legend title size
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),  
    legend.spacing.y = unit(0.6, "cm") 
  ) + guides(color = guide_legend(override.aes = list(size = 2)),  # Increases legend key sizes
             shape = guide_legend(override.aes = list(size = 2)))

ggsave("01_parameter/ROH_length_distribution.png", p6 , width = 12, height = 10, dpi = 300)

# Arrange the two plots together
ggarrange(final_plot_p3, final_plot_p2, p6,  ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)"))

# Save the plot
ggsave("01_parameter/ROH_analysis_plots_boxplot_len_dist.png",  width = 12, height = 8, dpi = 300)

#############################################################################



############# Till above the codes are used for publication i.e. only box plot ##################





########################### Voilin plot ##############################

##################### Length of ROH ####################
# Create the violin plot for Total ROH Length (Mb)
p1_violin <- ggviolin(
  roh_data, 
  x = "Group", 
  y = "Total_ROH_Mb", 
  fill = "Group", 
  palette = color_palette,
  add = "jitter",  # Adds scatter points to show individual data
  trim = FALSE  # Ensures full distribution is shown
)

# Statistical Testing: Bonferroni-adjusted t-tests
stat.test <- roh_data %>%
  t_test(Total_ROH_Mb ~ Group, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test <- stat.test %>%
  add_xy_position(x = "Group", step.increase = 0.2)

# Customize the plot with p-value annotations
final_plot_p1_violin <- p1_violin + 
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05
  ) +
  labs(
    title = "Total Length of ROHs (Mb)",
    x = "Population",
    y = "Total ROH Length (Mb)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Display the plot
final_plot_p1_violin

#################### FROH #################################

# Create the violin plot for FROH
p2_violin <- ggviolin(
  roh_data, 
  x = "Group", 
  y = "FROH", 
  fill = "Group", 
  palette = color_palette,
  add = "jitter",  # Adds scatter points to show individual data
  trim = FALSE  # Ensures full distribution is shown
)

# Statistical Testing: Bonferroni-adjusted t-tests
stat.test_froh <- roh_data %>%
  t_test(FROH ~ Group, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test_froh <- stat.test_froh %>%
  add_xy_position(x = "Group", step.increase = 0.2)

# Customize the plot with p-value annotations
final_plot_p2_violin <- p2_violin + 
  stat_pvalue_manual(
    stat.test_froh, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05
  ) +
  labs(
    title = "F(ROH)",
    x = "Population",
    y = expression(F[ROH])  # Correctly formats FROH as a mathematical expression
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Display the plot
final_plot_p2_violin

############### Number of ROH #####################
# Create the violin plot
p3_violin <- ggviolin(
  roh_data, 
  x = "Group", 
  y = "NSEG", 
  fill = "Group", 
  palette = color_palette,
  add = "jitter",  # Adds scatter points to show individual data
  trim = FALSE  # Ensures full distribution is shown
)

# Statistical Testing: Bonferroni-adjusted t-tests
stat.test_NROH <- roh_data %>%
  t_test(NSEG ~ Group, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test_NROH <- stat.test_NROH %>%
  add_xy_position(x = "Group", step.increase = 0.2)

# Customize the plot with p-value annotations
final_plot_p3_violin <- p3_violin + 
  stat_pvalue_manual(
    stat.test_NROH, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05
  ) +
  labs(
    title = "Number of ROHs",
    x = "Population",
    y = "Number of ROHs"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) 

# Display the plot
final_plot_p3_violin



# Arrange the two plots together
ggarrange(final_plot_p1_violin, final_plot_p2_violin, final_plot_p3_violin, p6, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)"))

# Save the plot
ggsave("01_parameter/ROH_analysis_plots_voilin.png",  width = 12, height = 8, dpi = 300)

############### Contig level Froh ###########################
# Load Required Libraries
library(ggplot2)
library(dplyr)
library(tidyr)


# Read PLINK ROH output
roh_contig <- read.table("01_parameter/ROH_results_up.hom", header = TRUE, stringsAsFactors = FALSE)
colnames(roh_contig)[2] <- "IID"  # Rename Individual column

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
roh_data$Population <- factor(roh_data$Population, levels = c("DUN", "HAM", "IRE", "LIN", "MAN"))

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

roh_data$Contig <- factor(roh_data$Contig, 
                      levels = c("contig_1", "contig_83", "contig_2", "contig_84", 
                                 "contig_3", "contig_85", "contig_4", "contig_86"))

# Boxplot of Froh by Contig for Each Population
p4 <- ggplot(roh_data, aes(x = Contig, y = Froh, fill = Population)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~Population, ncol = 1, scales = "free") +  # Separate panel per population
  scale_fill_manual(values = c("DUN" = "#E41A1C", "HAM" = "#A6D854", 
                               "IRE" = "#377EB8", "LIN" = "#984EA3", 
                               "MAN" = "#FF69B4")) + # Custom colors
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Froh Across Contigs for Each Population",
    x = "Contig",
    y = expression(F[ROH])
  )

# Save the Plot
ggsave("01_parameter/Froh_by_Contig.png", p4, width = 10, height = 10, dpi = 300)


######### Density plot ###############
roh_data <- read.table("01_parameter/ROH_results.hom_up.indiv", header = TRUE, stringsAsFactors = FALSE)
colnames(roh_data)[2] <- "IID"  # Rename Individual column

# Ensure population names are assigned (Modify this based on your population groups)
roh_data$Group <- factor(roh_data$FID, 
                         levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Convert total length from base pairs to Megabases
roh_data$Total_ROH_Mb <- roh_data$KB / 1000  # Convert KB to Mb

# Calculate FROH (Proportion of Genome in ROH)
# Assuming the genome size is ~129Mb (Modify if necessary)
genome_size_mb <- 126
roh_data$FROH <- roh_data$Total_ROH_Mb / genome_size_mb
# Filter each population subset
roh_man <- roh_data %>% filter(FID == "MAN")
roh_ham <- roh_data %>% filter(FID == "HAM")
roh_lin <- roh_data %>% filter(FID == "LIN")
roh_dun <- roh_data %>% filter(FID == "DUN")
roh_ire <- roh_data %>% filter(FID == "IRE")

# Create a single long-format data frame to avoid unequal column lengths
inbreeding_coeff <- bind_rows(
  roh_man %>% mutate(Coefficient = "man_FROH"),
  roh_ham %>% mutate(Coefficient = "ham_FROH"),
  roh_lin %>% mutate(Coefficient = "lin_FROH"),
  roh_dun %>% mutate(Coefficient = "dun_FROH"),
  roh_ire %>% mutate(Coefficient = "ire_FROH")
) %>%
  select(Coefficient, FROH) %>%  # Select relevant columns
  rename(Value = FROH)           # Rename for ggplot

# Create a mapping of old facet labels to new ones
facet_labels <- c(
  "dun_FROH" = "Dunedin",
  "ham_FROH" = "Hamilton",
  "ire_FROH" = "Ireland",
  "lin_FROH" = "Lincoln",
  "man_FROH" = "Mangōnui"
)

# Plot with updated facet labels
p5 <- ggplot(inbreeding_coeff, aes(x = Value, fill = Coefficient)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Coefficient, scales = "free", labeller = labeller(Coefficient = facet_labels)) +
  theme_minimal() +
  labs(title = "F(ROH)", x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

 ggsave("01_parameter/Froh_density.png", p5 , width = 10, height = 10, dpi = 300)
 

