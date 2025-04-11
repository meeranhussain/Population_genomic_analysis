library(ggplot2)
library(reshape2)

# Define paths to your Q files and metadata
q_files <- list.files("./Maethio_Qfiles", full.names = TRUE)
metadata_file <- "pop_map.txt"  # Update with your metadata file path

# Read metadata
metadata <- read.table(metadata_file, header = FALSE, stringsAsFactors = FALSE)
colnames(metadata) <- c("Sample", "Population")

# Find the maximum number of clusters across all Q files
max_clusters <- max(sapply(q_files, function(file) ncol(read.table(file, header = FALSE))))

# Define a function to process each Q file
process_q_file <- function(file, metadata, max_clusters) {
  q_data <- read.table(file, header = FALSE)  # Read Q file
  n_clusters <- ncol(q_data)  # Number of clusters in this file
  
  # Pad with zeros if the number of clusters is less than max_clusters
  if (n_clusters < max_clusters) {
    q_data <- cbind(q_data, matrix(0, nrow = nrow(q_data), ncol = max_clusters - n_clusters))
  }
  
  colnames(q_data) <- paste0("Cluster", seq_len(max_clusters))  # Name clusters
  q_data$Sample <- metadata$Sample  # Add sample names
  q_data$Population <- metadata$Population  # Add population info
  q_data <- q_data[order(q_data$Population), ]  # Sort by population
  q_data$File <- gsub(".*/|\\.Q$", "", file)  # Add filename for grouping
  return(q_data)
}

# Process all Q files and combine into a single data frame
admixture_data <- do.call(rbind, lapply(q_files, process_q_file, metadata = metadata, max_clusters = max_clusters))

# Reshape data for ggplot
admixture_melt <- melt(admixture_data, id.vars = c("Sample", "Population", "File"),
                       variable.name = "Cluster", value.name = "Proportion")

admixture_melt$File[which(admixture_melt$File == "M_aethio_r16.2")] <- "K2"
admixture_melt$File[which(admixture_melt$File == "M_aethio_r42.3")] <- "K3"
admixture_melt$File[which(admixture_melt$File == "M_aethio_r16.4")] <- "K4"
admixture_melt$File[which(admixture_melt$File == "M_aethio_r44.5")] <- "K5"



admixture_melt$File <- factor(admixture_melt$File, levels = c("K2", "K3", "K4", "K5"))
admixture_melt$Population <- factor(admixture_melt$Population, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Plot the admixture proportions
plot <- ggplot(admixture_melt, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(File ~ Population, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Paired") +  # Adjust palette if needed
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
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
        y = NULL,
        title = "") +
  scale_x_discrete(
    breaks = unique(admixture_melt$Sample),  # Ensure unique sample names
    labels = unique(admixture_melt$Sample)   # Match labels with unique breaks
  )
ggsave(path = "./", plot = plot, filename = paste("Maethio", "admixture.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)
################################################################################



ggarrange(pca, plot, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))
# Save the plot
ggsave("pca_admixture.png",  width = 16, height = 8, dpi = 600)
