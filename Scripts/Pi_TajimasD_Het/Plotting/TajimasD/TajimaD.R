##Plotting Sliding windows analysis - TajimaD
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(ggpubr)
library(rstatix)


#load datasets of each population 
HAMTajimaD <- read.table("../TajimasD/HAM_tjd_100kb.Tajima.D", header=T)
MANTajimaD <- read.table("../TajimasD/MAN_tjd_100kb.Tajima.D", header=T)
LINTajimaD <- read.table("../TajimasD/LIN_tjd_100kb.Tajima.D", header=T)
DUNTajimaD <- read.table("../TajimasD/DUN_tjd_100kb.Tajima.D", header=T)
IRETajimaD <- read.table("../TajimasD/IRE_tjd_100kb.Tajima.D", header=T)



# Remove rows with NA or "nan" in TajimaD column
HAMTajimaD <- HAMTajimaD[!is.na(HAMTajimaD$TajimaD) & HAMTajimaD$TajimaD != "NaN", ]
MANTajimaD <- MANTajimaD[!is.na(MANTajimaD$TajimaD) & MANTajimaD$TajimaD != "NaN", ]
LINTajimaD <- LINTajimaD[!is.na(LINTajimaD$TajimaD) & LINTajimaD$TajimaD != "NaN", ]
DUNTajimaD <- DUNTajimaD[!is.na(DUNTajimaD$TajimaD) & DUNTajimaD$TajimaD != "NaN", ]
IRETajimaD <- IRETajimaD[!is.na(IRETajimaD$TajimaD) & IRETajimaD$TajimaD != "NaN", ]



#calculating median and 5th and 95th percentile
median(HAMTajimaD$TajimaD)
quantile(HAMTajimaD$TajimaD, probs = c(0.05,0.95))
summary(HAMTajimaD$TajimaD)


median(MANTajimaD$TajimaD)
quantile(MANTajimaD$TajimaD, probs = c(0.05,0.95))
summary(MANTajimaD$TajimaD)


median(LINTajimaD$TajimaD)
quantile(LINTajimaD$TajimaD, probs = c(0.05,0.95))
summary(LINTajimaD$TajimaD)


median(DUNTajimaD$TajimaD)
quantile(DUNTajimaD$TajimaD, probs = c(0.05,0.95))
summary (DUNTajimaD$TajimaD)

median(IRETajimaD$TajimaD)
quantile(IRETajimaD$TajimaD, probs = c(0.05,0.95))
summary(IRETajimaD$TajimaD)



###VIOLIN PLOTS

length_HAMTajimaD <- length(HAMTajimaD$TajimaD)
length_MANTajimaD <- length(MANTajimaD$TajimaD)
length_LINTajimaD <- length(LINTajimaD$TajimaD)
length_DUNTajimaD <- length(DUNTajimaD$TajimaD)
length_IRETajimaD <- length(IRETajimaD$TajimaD)


#create dataset
data <- data.frame(
  name=c(rep("MAN",length_MANTajimaD), rep("HAM",length_HAMTajimaD), rep("LIN",length_LINTajimaD), rep("DUN",length_DUNTajimaD) 
         ,rep("IRE",length_IRETajimaD)),
  value=c(MANTajimaD$TajimaD,HAMTajimaD$TajimaD,LINTajimaD$TajimaD,DUNTajimaD$TajimaD,IRETajimaD$TajimaD))



############# Voilin plot #####################################
# Create the violin plot for Total ROH Length (Mb)
color_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")

# Ensure factors match
data$name <- factor(data$name, levels = c( "MAN", "HAM", "LIN", "DUN", "IRE"))

# Perform pairwise t-tests
stat.test <- data %>% 
  wilcox_test(value ~ name, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test <- stat.test %>%
  add_xy_position(x = "name", step.increase = 0.25)

# Filter out "ns" from stat.test
stat.test <- stat.test[stat.test$p.adj.signif != "ns", ]

####### Set p_value position
#stat.test$y.position <- seq(4.4, 5.8, length.out = nrow(stat.test))

###### Set p_value position manually
print(stat.test[10])
stat.test[1,10] <- 4.5
stat.test[2,10] <- 4.8
stat.test[3,10] <- 5.1
stat.test[4,10] <- 5.4
stat.test[5,10] <- 5.7



# Plot using ggviolin
p2 <- ggviolin(data, x = "name", y = "value", fill = "name",
         add = "boxplot", add.params = list(fill = "white", width = 0.1, alpha = 0.4)) +
  stat_summary(fun=median, geom='point', colour='black', size=2) +  # Add median points
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")) +  # Custom colors
  theme_minimal() +
  theme(
    text = element_text(size = 22),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  labs( x = "Population",  y = "Tajima's D", title = "") +
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",  # Show significance levels
    tip.length = 0.01, 
    bracket.shorten = 0.05,
    size = 8
  )

ggsave("TajimasD_voilinplot.png",p1, width = 14, height = 10, dpi = 300)

################################################################################

######## Export p_value #############
# Perform pairwise t-tests
stat.test <- data %>% 
  wilcox_test(value ~ name, p.adjust.method = "bonferroni")

# Adjust y positions for step-like annotations
stat.test <- stat.test %>%
  add_xy_position(x = "name", step.increase = 0.25)

# Convert list columns to character or numeric for export
stat.test <- stat.test %>%
  mutate(across(where(is.list), ~sapply(., paste, collapse = ", ")))

## export stats 
write.table(stat.test, file = "100kb_TajimasD_pvalue_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)





######################## Manhattan plot ########################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define scaffold sizes (predefined order)
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

# Compute cumulative starting positions for each scaffold
cumulative_starts <- c(0, cumsum(scaffold_lengths[-length(scaffold_lengths)]))
names(cumulative_starts) <- names(scaffold_lengths)

# Function to add cumulative genomic positions for Manhattan plot
add_cumulative_position <- function(data, scaffold_column = "CHROM", position_column = "BIN_START") {
  data %>%
    filter(!!sym(scaffold_column) %in% names(cumulative_starts)) %>%
    mutate(
      cumulative_position = !!sym(position_column) + cumulative_starts[as.character(!!sym(scaffold_column))]
    )
}

# Load and clean data for each population
HAM <- read.table("HAM_tjd_100kb.Tajima.D", header=T) %>% mutate(Population = "HAM")
MAN <- read.table("MAN_tjd_100kb.Tajima.D", header=T) %>% mutate(Population = "MAN")
LIN <- read.table("LIN_tjd_100kb.Tajima.D", header=T) %>% mutate(Population = "LIN")
DUN <- read.table("DUN_tjd_100kb.Tajima.D", header=T) %>% mutate(Population = "DUN")
IRE <- read.table("IRE_tjd_100kb.Tajima.D", header=T) %>% mutate(Population = "IRE")

# Combine data
tajima_all <- bind_rows(HAM, MAN, LIN, DUN, IRE)

# Remove NaN or NA values
tajima_all <- tajima_all %>%
  filter(!is.na(TajimaD) & TajimaD != "NaN")

# Add cumulative position for Manhattan plot
tajima_all <- add_cumulative_position(tajima_all)

# Ensure scaffolds are ordered by size in the plot
tajima_all$CHROM <- factor(tajima_all$CHROM, levels = names(scaffold_lengths))

# Manhattan plot
ggplot(tajima_all, aes(x = cumulative_position / 1e6, y = TajimaD, color = TajimaD > 0)) +
  geom_point(size=0.8, alpha=0.7) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("Negative", "Positive"),
                     name = "Tajima's D") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=0.3, color="grey20") +
  facet_wrap(~ Population, ncol=1, scales="free_y") +  # Separate by population, keeping a single x-axis
  theme_bw() +
  labs(x = "Genomic Position (Mb)", y = "Tajima's D",
       title = "Manhattan Plot of Tajima's D (10kb windows) with Ordered Contigs") +
  theme(strip.text = element_text(face="bold"),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
