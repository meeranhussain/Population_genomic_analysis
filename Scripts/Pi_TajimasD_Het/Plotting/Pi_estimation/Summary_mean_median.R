# Define populations
populations <- c("HAM", "MAN", "LIN", "DUN", "IRE")

# Create an empty data frame to store results
pi_summary <- data.frame(Population = character(),
                         Median = numeric(),
                         Percentile_5th = numeric(),
                         Percentile_95th = numeric(),
                         Min = numeric(),
                         Max = numeric(),
                         Mean = numeric(),
                         stringsAsFactors = FALSE)

# Loop through each population
for (pop in populations) {
  # Read the Pi estimate file (assuming it's tab-delimited)
  pi_data <- read.table(paste0(pop, "_pi_100kb.windowed.pi"), header = TRUE)
  
  # Extract Pi values (assuming column name is "PI")
  pi_values <- pi_data$PI
  
  # Calculate statistics
  med <- median(pi_values, na.rm = TRUE)
  perc_5th <- quantile(pi_values, probs = 0.05, na.rm = TRUE)
  perc_95th <- quantile(pi_values, probs = 0.95, na.rm = TRUE)
  min_val <- min(pi_values, na.rm = TRUE)
  max_val <- max(pi_values, na.rm = TRUE)
  mean_val <- mean(pi_values, na.rm = TRUE)
  
  # Store in the summary data frame
  pi_summary <- rbind(pi_summary, data.frame(Population = pop,
                                             Median = med,
                                             Percentile_5th = perc_5th,
                                             Percentile_95th = perc_95th,
                                             Min = min_val,
                                             Max = max_val,
                                             Mean = mean_val))
}

# Export to a file
write.table(pi_summary, file = "100kb_Pi_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary to console
print(pi_summary)
