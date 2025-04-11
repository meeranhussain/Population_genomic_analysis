# Define populations
populations <- c("HAM", "MAN", "LIN", "DUN", "IRE")

# Create an empty data frame to store results
summary_results <- data.frame(Population = character(),
                              Median = numeric(),
                              Percentile_5th = numeric(),
                              Percentile_95th = numeric(),
                              Min = numeric(),
                              Max = numeric(),
                              Mean = numeric(),
                              stringsAsFactors = FALSE)

# Loop through each population
for (pop in populations) {
  # Read the TajimaD file (assuming it's a tab-delimited file)
  tajima_data <- read.table(paste0(pop, "_tjd_100kb.Tajima.D"), header = TRUE)
  
  # Extract Tajima’s D values
  tajima_values <- tajima_data$TajimaD
  
  # Calculate statistics
  med <- median(tajima_values, na.rm = TRUE)
  perc_5th <- quantile(tajima_values, probs = 0.05, na.rm = TRUE)
  perc_95th <- quantile(tajima_values, probs = 0.95, na.rm = TRUE)
  min_val <- min(tajima_values, na.rm = TRUE)
  max_val <- max(tajima_values, na.rm = TRUE)
  mean_val <- mean(tajima_values, na.rm = TRUE)
  
  # Store in the summary data frame
  summary_results <- rbind(summary_results, data.frame(Population = pop,
                                                       Median = med,
                                                       Percentile_5th = perc_5th,
                                                       Percentile_95th = perc_95th,
                                                       Min = min_val,
                                                       Max = max_val,
                                                       Mean = mean_val))
}

# Export to a file
write.table(summary_results, file = "100kb_TajimaD_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary to console
print(summary_results)
