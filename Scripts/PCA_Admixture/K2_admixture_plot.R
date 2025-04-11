library(ggplot2)
library(reshape2)

# Define paths to your Q files and metadata
q_file <- "./Maethio_Qfiles/M_aethio_r16.2.Q"  # Adjust to the correct K2 file path
metadata_file <- "pop_map.txt"  # Update with your metadata file path

# Read metadata
metadata <- read.table(metadata_file, header = FALSE, stringsAsFactors = FALSE)
colnames(metadata) <- c("Sample", "Population")

# Read K2 Q file
q_data <- read.table(q_file, header = FALSE)

# Assign cluster names
colnames(q_data) <- c("Cluster1", "Cluster2")
q_data$Sample <- metadata$Sample
q_data$Population <- metadata$Population

# Reshape data for ggplot
admixture_melt <- melt(q_data, id.vars = c("Sample", "Population"),
                       variable.name = "Cluster", value.name = "Proportion")

# Order populations
admixture_melt$Population <- factor(admixture_melt$Population, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Plot the admixture proportions for K2
plot <- ggplot(admixture_melt, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ Population, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Paired") +  # Adjust palette if needed
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),  # Ensure Y-axis values are visible
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.05, "lines"),
    panel.spacing.y = unit(0.01, "lines")
  ) +
  labs( x = NULL,
        y = "Admixture Proportion", 
        title = NULL) +
  scale_y_continuous(labels = scales::percent_format())  # Format Y-axis as percentage

# Save the plot
ggsave(path = "./", plot = plot, filename = "Maethio_K2_admixture.png", device='png', width = 10, height = 6, dpi = 300)
