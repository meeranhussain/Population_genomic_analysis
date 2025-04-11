# Load required library
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("topGO")

library(topGO)
library(dplyr)
library(ggplot2)

# Accept the argument passed to the script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide the comparison name as an argument, e.g., HAM_DUN.")
}

comparison <- args[1] # e.g., HAM_DUN
input_path <- paste0(comparison, "/02_fetch_GOs/") # Path to input and output directory
input_file <- paste0(input_path, comparison, "_trans_idssnpeff.txt") # Full path to input file
output_prefix <- paste0(input_path, "GOresults_", comparison) # Full path to output prefix

# Import tab-delimited file of all genes and their GO terms
geneID2GO <- readMappings("00_overall_vcf_ann/Formatted_GO_List_all_variant.txt")  

# Inspect the structure of geneID2GO (optional)
# str(head(geneID2GO))

geneNames <- names(geneID2GO)
# head(geneNames)

# Import transcript IDs of outlier SNPs
interesting_genes <- read.table(input_file)

myInterestingGenes <- as.vector(interesting_genes$V1) # Extract gene IDs as a vector
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
# str(geneList)

# Function to perform GO enrichment for a specific ontology
run_GO_enrichment <- function(ontology, geneList, geneID2GO, output_prefix) {
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  allRes <- GenTable(
    GOdata, 
    raw.p.value = resultFisher, 
    classicFisher = resultFisher, 
    ranksOf = "classicFisher", 
    Fis = resultFisher, 
    topNodes = length(resultFisher@score)
  )
  output_file <- paste0(output_prefix, "_", ontology, ".txt")
  write.table(allRes, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  return(allRes)
}

# Perform GO enrichment for each ontology
cat("Performing GO enrichment for Biological Process (BP)...\n")
bp_results <- run_GO_enrichment("BP", geneList, geneID2GO, output_prefix)

cat("Performing GO enrichment for Molecular Function (MF)...\n")
mf_results <- run_GO_enrichment("MF", geneList, geneID2GO, output_prefix)

cat("Performing GO enrichment for Cellular Component (CC)...\n")
cc_results <- run_GO_enrichment("CC", geneList, geneID2GO, output_prefix)

cat("GO enrichment analysis completed. Results saved to files.\n")


############### Plotting ################################
bp_results$Ontology <- "BP"
mf_results$Ontology <- "MF"
cc_results$Ontology <- "CC"

# Combine results from all three ontologies
combined_results <- bind_rows(bp_results, mf_results, cc_results)

# Filter for significant terms and select the top 5 for each ontology
top_results <- combined_results %>%
  filter(as.numeric(raw.p.value) < 0.05) %>%  # Apply p-value threshold
  group_by(Ontology) %>%
  arrange(as.numeric(`Rank in classicFisher`), .by_group = TRUE) %>%  # Use rank for sorting
  slice_head(n = 10) %>%  # Select top 5 GO terms based on rank
  ungroup() %>%
  mutate(Term = factor(Term, levels = unique(Term)))  # Preserve term order for plotting


# Plot top GO terms
GO_enrich <- ggplot(top_results, aes(x = Term, y = -log10(as.numeric(raw.p.value)), fill = Ontology)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = paste0(comparison, ":Top 10 Significant GO Terms"),
       x = "GO Term",
       y = "-log10(p-value)",
       fill = "Ontology") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face="bold"),
    plot.title = element_text(hjust = 0.5)) 

ggsave(path = input_path , plot = GO_enrich, filename = paste0(comparison, "_top_10_GO_enrich.png"), device='png', width = 12, height = 10, dpi = 300)
