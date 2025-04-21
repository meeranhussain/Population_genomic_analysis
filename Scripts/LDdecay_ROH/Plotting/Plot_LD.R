library(ggplot2)
library(readr)
library(dplyr)

# Define file paths
files <- list(
  "MAN" = "plot_all_rep1.MANGONUI",
  "HAM" = "plot_all_rep1.HAMILTON",
  "LIN" = "plot_all_rep1.LINCOLN",
  "DUN" = "plot_all_rep1.DUNEDIN",
  "IRE" = "plot_all_rep1.IRELAND"
)

# Define color palette in specified order
color_palette <- c("MAN" = "#66c2a5", "HAM" = "#fc8d62", "LIN" = "#8da0cb", "DUN" = "#e78ac3", "IRE" = "#a6d854")

# Read and combine data into one dataframe
ld_data <- do.call(rbind, lapply(names(files), function(pop) {
  df <- read_table(files[[pop]], col_names = TRUE, col_types = cols(
    Dist = col_double(),
    `Mean_r^2` = col_double(),
    `Mean_D'` = col_skip(),  # Skip NA column
    `Sum_r^2` = col_skip(),
    `Sum_D'` = col_skip(),
    NumberPairs = col_skip()
  ))
  
  df <- df %>%
    rename(Distance = Dist, r2 = `Mean_r^2`) %>%
    mutate(Distance = Distance / 1000, Population = pop)  # Convert Distance to Kb
  
  return(df)
}))

# Set factor levels for ordered legend
ld_data$Population <- factor(ld_data$Population, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Create ggplot with custom styling
p <- ggplot(ld_data, aes(x = Distance, y = r2, color = Population)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_palette) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +  # Explicitly set breaks and limits
  labs(title = "",
       x = "Distance (Kb)",
       y = expression(r^2),
       color = "Population") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),  
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold"),
  )

# Save plot
ggsave("plot_all_rep1.pdf", p, width = 8, height = 6)
ggsave("plot_all_rep1.png", p, width = 8, height = 6, dpi = 600)
