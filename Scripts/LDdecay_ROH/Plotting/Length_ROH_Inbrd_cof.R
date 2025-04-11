# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)


# Sample data: replace with your actual PLINK ROH results
roh_data <- read.delim("01_parameter/ROH_results_up.hom", sep = " ")
colnames(roh_data)[2] <- "IID"


# Convert contig names to factor for ordered plotting
roh_data$CHR <- factor(roh_data$CHR, 
                       levels = c("contig_1", "contig_83", "contig_2", "contig_84", 
                                  "contig_3", "contig_85", "contig_4", "contig_86"))

# Summarize the number of ROH segments per contig per population
roh_summary <- roh_data %>%
  group_by(FID, CHR) %>%
  summarise(Count = n(), .groups = "drop")

roh_summary$FID <- factor(roh_summary$FID, levels = c("MAN", "HAM", "LIN", "DUN", "IRE"))

# Plot with a shared y-axis using facet_grid()
p1 <- ggplot(roh_summary, aes(x = CHR, y = Count, fill = FID)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(FID), scales = "free_x", space = "free_x") +  # Shared y-axis
  theme_minimal() +
  labs(title = "Number of ROH Segments", x = "Scaffolds", y = "Count", fill = "Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank()) +
  scale_fill_manual(values = c("MAN" = "#E41A1C", "HAM" = "#A6D854", 
                               "LIN" = "#377EB8", "DUN" = "#984EA3", 
                               "IRE" = "#FF69B4"))  # Custom colors


######### FROH Population based density plot ######################
contig_sizes <- data.frame(
  Contig = c("contig_1", "contig_83", "contig_2", "contig_84", "contig_3", 
             "contig_85", "contig_4", "contig_86"),
  contig_sizes = c(29470290, 28376149, 23977460, 11926033, 
                   10466058, 9529246, 7912293, 7312752)
)
########### MANGONUI #########################

roh_man <- roh_data %>% filter(FID == "MAN" )

# Compute FROH (Sum of ROH KB per individual / total genome size estimate)
genome_size <- sum(contig_sizes$contig_sizes) / 1000  # in KB
roh_man$FROH <- roh_man$KB / genome_size

# Example inbreeding coefficient values (replace with actual calculations)
inbreeding_coeff <- data.frame(
  FROH = roh_man$FROH
)

# Convert to long format for ggplot
inbreeding_long <- pivot_longer(inbreeding_coeff, cols = everything(), names_to = "Coefficient", values_to = "Value")

# Plot 2: Density plots for different inbreeding coefficients
p2 <- ggplot(inbreeding_long, aes(x = Value, fill = Coefficient)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Coefficient, scales = "free") +
  theme_minimal() +
  labs(title = "Mangōnui", x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


########### HAMILTON #########################################

roh_ham <- roh_data %>% filter(FID == "HAM" )

# Compute FROH (Sum of ROH KB per individual / total genome size estimate)
genome_size <- sum(contig_sizes$contig_sizes) / 1000  # in KB
roh_ham$FROH <- roh_ham$KB / genome_size

# Example inbreeding coefficient values (replace with actual calculations)
inbreeding_coeff <- data.frame(
  FROH = roh_ham$FROH
)

# Convert to long format for ggplot
inbreeding_long <- pivot_longer(inbreeding_coeff, cols = everything(), names_to = "Coefficient", values_to = "Value")

# Plot 2: Density plots for different inbreeding coefficients
p3 <- ggplot(inbreeding_long, aes(x = Value, fill = Coefficient)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Coefficient, scales = "free") +
  theme_minimal() +
  labs(title = "Hamilton", x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


########### LINCOLN #########################

roh_lin <- roh_data %>% filter(FID == "LIN" )

# Compute FROH (Sum of ROH KB per individual / total genome size estimate)
genome_size <- sum(contig_sizes$contig_sizes) / 1000  # in KB
roh_lin$FROH <- roh_lin$KB / genome_size

# Example inbreeding coefficient values (replace with actual calculations)
inbreeding_coeff <- data.frame(
  FROH = roh_lin$FROH,
)

# Convert to long format for ggplot
inbreeding_long <- pivot_longer(inbreeding_coeff, cols = everything(), names_to = "Coefficient", values_to = "Value")

# Plot 2: Density plots for different inbreeding coefficients
p4 <- ggplot(inbreeding_long, aes(x = Value, fill = Coefficient)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Coefficient, scales = "free") +
  theme_minimal() +
  labs(title = "Lincoln", x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


########### DUNEDIN #########################

roh_dun <- roh_data %>% filter(FID == "DUN" )

# Compute FROH (Sum of ROH KB per individual / total genome size estimate)
genome_size <- sum(contig_sizes$contig_sizes) / 1000  # in KB
roh_dun$FROH <- roh_dun$KB / genome_size

# Example inbreeding coefficient values (replace with actual calculations)
inbreeding_coeff <- data.frame(
  FROH = roh_dun$FROH
)

# Convert to long format for ggplot
inbreeding_long <- pivot_longer(inbreeding_coeff, cols = everything(), names_to = "Coefficient", values_to = "Value")

# Plot 2: Density plots for different inbreeding coefficients
p5 <- ggplot(inbreeding_long, aes(x = Value, fill = Coefficient)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Coefficient, scales = "free") +
  theme_minimal() +
  labs(title = "Dunedin", x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


########### IRELAND #########################

roh_ire <- roh_data %>% filter(FID == "IRE" )

# Compute FROH (Sum of ROH KB per individual / total genome size estimate)
genome_size <- sum(contig_sizes$contig_sizes) / 1000  # in KB
roh_ire$FROH <- roh_ire$KB / genome_size

# Example inbreeding coefficient values (replace with actual calculations)
inbreeding_coeff <- data.frame(
  FROH = roh_ire$FROH
)

# Convert to long format for ggplot
inbreeding_long <- pivot_longer(inbreeding_coeff, cols = everything(), names_to = "Coefficient", values_to = "Value")

# Plot 2: Density plots for different inbreeding coefficients
p6 <- ggplot(inbreeding_long, aes(x = Value, fill = Coefficient)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Coefficient, scales = "free") +
  theme_minimal() +
  labs(title = "Ireland", x = "Coefficient Value", y = "Density") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

################## SAVE PLOTS #############################
ggsave(path = "./01_parameter/", plot = p1, filename = paste("All_pop", "ROH_length.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)
ggsave(path = "./01_parameter/", plot = p2, filename = paste("MAN", "inbrd_cof.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)
ggsave(path = "./01_parameter/", plot = p3, filename = paste("HAM", "inbrd_cof.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)
ggsave(path = "./01_parameter/", plot = p4, filename = paste("LIN", "inbrd_cof.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)
ggsave(path = "./01_parameter/", plot = p5, filename = paste("DUN", "inbrd_cof.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)
ggsave(path = "./01_parameter/", plot = p6, filename = paste("IRE", "inbrd_cof.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)

############## GGRANGE ########################
p7 <- ggarrange(p2, p3, p4, p5, p6, ncol = 3, nrow = 2, align = "hv",         
          heights = c(1, 1),    # Ensure equal heights for rows
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)"))
ggsave("01_parameter/Inbreeding_coefficent_all_pop_plots.png", p7,   width = 12, height = 10, dpi = 300)
