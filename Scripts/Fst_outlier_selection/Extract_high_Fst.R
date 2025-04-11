# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Define input file paths
fst_file <- "pixy_100kb_LD_fst.txt"

# Load data
fst_data <- read.csv(fst_file, sep = "\t", header = TRUE)

# Define the scaffold lengths in a named vector

scaffold <- c(
  "contig_1", "contig_83", "contig_2", "contig_84", 
  "contig_3", "contig_85", "contig_4", "contig_86"
)

###### Subsetting scaffolds ##############
fst_data <- subset(fst_data, chromosome %in% scaffold)

# Define scaffold lengths as a named vector
scaffold_lengths <- c(
  "contig_1" = 29470290,
  "contig_83" = 28376149,
  "contig_2" = 23977460,
  "contig_84" = 11926033,
  "contig_3" = 10466058,
  "contig_85" = 9529246,
  "contig_4" = 7912293,
  "contig_86" = 7312752
)

# Calculate cumulative start positions for each scaffold
cumulative_starts <- c(0, cumsum(scaffold_lengths[-length(scaffold_lengths)]))
names(cumulative_starts) <- names(scaffold_lengths)

# Add cumulative position to each dataset
add_cumulative_position <- function(data, scaffold_column, position_column) {
  data <- data %>%
    mutate(
      cumulative_position = window_pos_1 + cumulative_starts[as.character(!!sym(scaffold_column))]
    )
  return(data)
}

################################## Hamilton vs Mangonui ########################
############### Fetch Population ###############################
fst_ham_man <- fst_data %>% filter(pop1 == "HAM" & pop2 == "MAN")

### Add Position ###################################
fst_ham_man <- add_cumulative_position(fst_ham_man, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ham_man <- fst_ham_man %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ham_man$avg_wc_fst, 0.99)


ham_man <- ggplot(fst_ham_man, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ham_man$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ham_man$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Hamilton vs Mangōnui : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ham_man, filename = paste("HAM_MAN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ham_man %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "HAM_MAN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################################# Hamilton vs Lincoln ############################
############### Fetch Population ###############################
fst_ham_lin <- fst_data %>% filter(pop1 == "HAM" & pop2 == "LIN")

### Add Position ###################################
fst_ham_lin <- add_cumulative_position(fst_ham_lin, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ham_lin <- fst_ham_lin %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ham_lin$avg_wc_fst, 0.99)


ham_lin <- ggplot(fst_ham_lin, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ham_lin$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ham_lin$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Hamilton vs Lincoln : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ham_lin, filename = paste("HAM_LIN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ham_lin %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "HAM_LIN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##################################################################################


################################# Hamilton vs Dunedin ############################
############### Fetch Population ###############################
fst_ham_dun <- fst_data %>% filter(pop1 == "DUN" & pop2 == "HAM")

### Add Position ###################################
fst_ham_dun <- add_cumulative_position(fst_ham_dun, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ham_dun <- fst_ham_dun %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ham_dun$avg_wc_fst, 0.99)


ham_dun <- ggplot(fst_ham_dun, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ham_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ham_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Hamilton vs Dunedin : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ham_dun, filename = paste("HAM_DUN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ham_dun %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "HAM_DUN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)



################################# Hamilton vs Ireland ############################
############### Fetch Population ###############################
fst_ham_ire <- fst_data %>% filter(pop1 == "HAM" & pop2 == "IRE")

### Add Position ###################################
fst_ham_ire <- add_cumulative_position(fst_ham_ire, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ham_ire <- fst_ham_ire %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ham_ire$avg_wc_fst, 0.99)


ham_ire <- ggplot(fst_ham_ire, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ham_ire$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ham_ire$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Hamilton vs Ireland : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ham_ire, filename = paste("HAM_IRE_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ham_ire %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "HAM_IRE_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


##################### Ireland vs Mangonui ##########################################

############### Fetch Population ###############################
fst_ire_man <- fst_data %>% filter(pop1 == "IRE" & pop2 == "MAN")

### Add Position ###################################
fst_ire_man <- add_cumulative_position(fst_ire_man, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ire_man <- fst_ire_man %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ire_man$avg_wc_fst, 0.99)


ire_man <- ggplot(fst_ire_man, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ire_man$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ire_man$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Ireland vs Mangōnui : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ire_man, filename = paste("IRE_MAN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ire_man %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "IRE_MAN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################################## Ireland vs Lincoln ###############################

############### Fetch Population ###############################
fst_ire_lin <- fst_data %>% filter(pop1 == "IRE" & pop2 == "LIN")

### Add Position ###################################
fst_ire_lin <- add_cumulative_position(fst_ire_lin, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ire_lin <- fst_ire_lin %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ire_lin$avg_wc_fst, 0.99)


ire_lin <- ggplot(fst_ire_lin, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ire_lin$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ire_lin$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Ireland vs Lincoln : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ire_lin, filename = paste("IRE_LIN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ire_lin %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "IRE_LIN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################################## Ireland vs Dunedin ###############################

############### Fetch Population ###############################
fst_ire_dun <- fst_data %>% filter(pop1 == "DUN" & pop2 == "IRE")

### Add Position ###################################
fst_ire_dun <- add_cumulative_position(fst_ire_dun, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_ire_dun <- fst_ire_dun %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_ire_dun$avg_wc_fst, 0.99)


ire_dun <- ggplot(fst_ire_dun, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_ire_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_ire_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Ireland vs Dunedin : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = ire_dun, filename = paste("IRE_DUN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_ire_dun %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "IRE_DUN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################################## Mangonui vs Dunedin ###############################

############### Fetch Population ###############################
fst_man_dun <- fst_data %>% filter(pop1 == "DUN" & pop2 == "MAN")

### Add Position ###################################
fst_man_dun <- add_cumulative_position(fst_man_dun, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_man_dun <- fst_man_dun %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_man_dun$avg_wc_fst, 0.99)


man_dun <- ggplot(fst_man_dun, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_man_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_man_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Mangōnui vs Dunedin : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = man_dun, filename = paste("MAN_DUN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_man_dun %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "MAN_DUN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################################## Mangonui vs Lincoln ###############################

############### Fetch Population ###############################
fst_man_lin <- fst_data %>% filter(pop1 == "LIN" & pop2 == "MAN")

### Add Position ###################################
fst_man_lin <- add_cumulative_position(fst_man_lin, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_man_lin <- fst_man_lin %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_man_lin$avg_wc_fst, 0.99)


man_lin <- ggplot(fst_man_lin, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_man_lin$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_man_lin$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Mangōnui vs Lincoln : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = man_lin, filename = paste("MAN_LIN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_man_lin %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "MAN_LIN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)


################################## Lincoln vs Dunedin ###############################

############### Fetch Population ###############################
fst_lin_dun <- fst_data %>% filter(pop1 == "DUN" & pop2 == "LIN")

### Add Position ###################################
fst_lin_dun <- add_cumulative_position(fst_lin_dun, "chromosome", "window_pos_1")


# Filter rows with non-NA Fst values
fst_lin_dun <- fst_lin_dun %>% filter(!is.na(avg_wc_fst))

# Calculate the Fst threshold for the Top 1%
fst_threshold <- quantile(fst_lin_dun$avg_wc_fst, 0.99)


lin_dun <- ggplot(fst_lin_dun, aes(x = cumulative_position / 1e6, y = avg_wc_fst)) + 
  geom_point() + 
  geom_hline(yintercept=fst_threshold, linetype="dashed", color = "red") +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(fst_lin_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25),
    labels = seq(0, max(fst_lin_dun$cumulative_position / 1e6, na.rm = TRUE), by = 25)
  ) +
  labs(
    x = "Position (Mb)",
    y = "Fst",
    title = "Lincoln vs Dunedin : Fst Across 8 scaffolds (Threshold:Top 1%)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


ggsave(path = "./", plot = lin_dun, filename = paste("LIN_DUN_high_fst", "Pixy.png", sep = "_"), device='png', width = 10, height = 6, dpi = 300)


# Extract high-Fst windows
high_fst_windows <- fst_lin_dun %>% filter(avg_wc_fst >= fst_threshold)

# Save high-Fst windows to file
write.table(high_fst_windows, file = "LIN_DUN_high_fst_windows_tp_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

