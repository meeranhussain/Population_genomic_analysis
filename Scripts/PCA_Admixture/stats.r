# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))
packages <- c("vcfR","ggplot2", "dplyr", "tidyr", "reshape", "StAMPP", "adegenet", "corrplot")

# Check if packages are installed
installed_packages <- packages %in% rownames(installed.packages())

# Install packages that are not installed
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

ainstall.packages('devtools',dependencies=T)
library(devtools)
install.packages("rlang")

#install pophelper package from GitHub
install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")

devtools::install_local("pophelper-master")


library(tidyr)
library(dplyr)
library(vcfR)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(corrplot)
library(reshape)
library(LEA)
library(pophelper)
library(stringr)

###
##running PCA
snp_vcf2 = read.vcfR("M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno_LD.vcf")
pop.data2 = read.table("pop_map.txt", header = F)
gl.snp2 <- vcfR2genlight(snp_vcf2)
pop(gl.snp2) <- rep(pop.data2$V2)
snp.pca2 <- glPca(gl.snp2, nf = 10)

##write PCA scores
snp.pca.scores2 <- as.data.frame(snp.pca2$scores)
snp.pca.scores2$pop <- pop(gl.snp2)
PCA_adegen_path <- file.path("./", paste("M_aethio", "adegenetPCA.txt", sep = "_"))
write.table(snp.pca.scores2, PCA_adegen_path, sep = "\t")

##export list of eigen values and percentage variances for each PC
#eigen values for PCs
eig.val<-snp.pca2$eig
eig.val
#percentages of variance for each PC
eig.perc <- 100*snp.pca2$eig/sum(snp.pca2$eig)
eig.perc
eigen<-data.frame(eig.val,eig.perc)
eigen
#writing file with both
PCA_eigen_path <- file.path("./", paste("M_aethio", "eigen_summary.txt", sep = "_"))

write.csv(eigen,file=PCA_eigen_path,row.names=TRUE,quote=FALSE)


###Load Data
data2 <- read.delim("M_aethio_adegenetPCA.txt") #I have manually added population information to this file prior to loading it
mycol = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")

#####Grouping Samples in Order
group_order <- c('MAN', 'HAM', 'LIN', 'DUN', 'IRE')

data2$pop <- factor(data2$pop, levels = group_order)

data2 <- data2 %>% arrange(pop)
##############################################

# Plot PC1 vs PC2
ggplot(data2, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 2) + 
  scale_color_manual(values = mycol) +
  theme_classic() +
  xlab(paste("PC1 (", round(100 * snp.pca2$eig[1] / sum(snp.pca2$eig), 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(100 * snp.pca2$eig[2] / sum(snp.pca2$eig), 2), "%)", sep = ""))


########## FOR Fancy plot ###########################

library(ggplot2)

#  PCA plot
pca <- ggplot(data2, aes(x = PC1, y = PC2, color = pop, shape = pop)) +
  geom_point(size = 5, alpha = 0.8) +  # Larger points with slight transparency
  scale_color_manual(values = mycol) +  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 19, 15)) +  # Different shapes for populations
  labs(
    x = paste0("PC1 (", round(100 * snp.pca2$eig[1] / sum(snp.pca2$eig), 2), "%)"),
    y = paste0("PC2 (", round(100 * snp.pca2$eig[2] / sum(snp.pca2$eig), 2), "%)"),
    color = "Population",
    shape = "Population"
  ) +
  theme_classic(base_size = 14) + 
  labs( 
        title = "") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +  # Horizontal line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")    # Vertical line at 0

ggsave(path = "./", plot = pca, filename = paste("Maethio", "PCA.png", sep = "_"), device='png', width = 8, height = 6, dpi = 300)


# Calculate the cumulative variance explained
cum_var <- cumsum(snp.pca2$eig) / sum(snp.pca2$eig)

# Plot cumulative variance
plot(cum_var, type = "b", xlab = "Principal Component", ylab = "Cumulative Variance Explained", main = "Cumulative Variance by PCs")


########################## MDS Plot ########################## 

# Assuming you have already loaded your VCF file and created the genlight object
# Example:
gl.snp = vcfR2genlight(snp_vcf2)

# 1. Compute a distance matrix (Euclidean or Manhattan)
# Euclidean distance (default)
dist_matrix <- dist(as.matrix(gl.snp))

# 2. Perform MDS
mds_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# k = 2 means you want to reduce to 2 dimensions
# eig = TRUE tells R to return the eigenvalues

# 3. Prepare data for plotting
mds_df <- as.data.frame(mds_result$points)
pop(gl.snp) <- rep(pop.data2$V2)
mds_df$pop <- pop(gl.snp)  

#####Grouping Samples in Order
group_order <- c('MAN', 'HAM', 'LIN', 'DUN', 'IRE')

mds_df$pop <- factor(mds_df$pop, levels = group_order)

mds_df <- mds_df %>% arrange(pop)
##############################################

# Calculate variance explained by each dimension
mds_variance <- mds_result$eig / sum(mds_result$eig) * 100

# Format labels for x and y axes
x_label <- paste0("MDS1 (", round(mds_variance[1], 2), "%)")
y_label <- paste0("MDS2 (", round(mds_variance[2], 2), "%)")

# 4. Plot the MDS results
ggplot(mds_df, aes(x = V1, y = V2, color = pop)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f1c039", "#f37d21", "#51692d", "#56ba32", "#487aa1",  "#e41a1c", "#02ccf5")) +
  theme_classic() +
  xlab(x_label) +
  ylab(y_label) +
  ggtitle("MDS Plot") +
  theme(plot.title = element_text(hjust = 0.5))


######################################################### PAIRWISE Fst ###################################################################

##calculating Fst between populations
mae_fst <- stamppFst(gl.snp, nboots = 100, percent = 95, nclusters = 6)
Fst <- mae_fst$Fsts
pFst <- mae_fst$Pvalues
write.table(Fst, "Fst.txt", sep="\t")
write.table(pFst, "Fst_pvalue.txt", sep="\t")

##creating heatmap
# Melt the correlation matrix
M_Fst <- as.matrix(read.table("Fst.txt"))
MaeFs <- melt(M_Fst, na.rm = TRUE)
summary(MaeFs$value)


ggplot(data = MaeFs, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#ffd60a", high = "#001d3d", mid = "#4e9de6", 
                       midpoint = 0.056, limit = c(0.005,0.11), space = "Lab", 
                       name="Pairwise Fst") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  labs(x = NULL, y = NULL) +
  coord_fixed()

ggsave(path = "./", filename = paste("M_aethio_popgen_all", "pairwiseFstplot.png", sep = "_"), device='jpeg', dpi=1000)
##########################################################################################################################################


################################### SNMF ############################################
library(LEA)
library(pophelper)

##creating input files
vcf2geno(input.file = "M_aethio_MAN_HAM_LIN_DUN_IRE_bialSNP_MAF_geno_LD.vcf", output.file = "M_aethio.geno")

##snmf clustering
projectalpha = NULL
projectalpha = snmf("M_aethio.geno", K = 1:10, repetitions = 50, entropy = T, CPU = 8, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
pdf(file = "./cross_ent_alphadefualt.pdf")
plot(projectalpha, col = "maroon4", pch = 19, cex = 1.2)
dev.off()

best2 = which.min(cross.entropy(projectalpha, K = 2))
best2
best3 = which.min(cross.entropy(projectalpha, K = 3))
best3
best4 = which.min(cross.entropy(projectalpha, K = 4))
best4
best5 = which.min(cross.entropy(projectalpha, K = 5))
best5
best6 = which.min(cross.entropy(projectalpha, K = 6))
best6


# Create the folder to store the Q files
if (!dir.exists("./Maethio_Qfiles")) {
  dir.create("./Maethio_Qfiles")
}

# Find the best runs for each K
best2 = which.min(cross.entropy(projectalpha, K = 2))
best3 = which.min(cross.entropy(projectalpha, K = 3))
best4 = which.min(cross.entropy(projectalpha, K = 4))
best5 = which.min(cross.entropy(projectalpha, K = 5))
best6 = which.min(cross.entropy(projectalpha, K = 6))
best7 = which.min(cross.entropy(projectalpha, K = 7))
best8 = which.min(cross.entropy(projectalpha, K = 8))

# Store best runs in a list
best_runs <- list(best2, best3, best4, best5, best6, best7, best8)

# File paths to the Q files for each K
# Adjust the file pattern depending on how your Q files are named.
# Here, assuming the files are named as project_K#.Q or similar format
for (K in 2:8) {
  # Define the best run for this value of K
  best_run <- best_runs[[K-1]]
  
  # Path to the best Q file for this K
  # Corrected way to fetch the file path with pattern matching
  qfile_pattern <- file.path("M_aethio.snmf", paste0("K", K), paste0("run", best_run), "*.Q")
  qfile_path <- Sys.glob(qfile_pattern)
  
  # If the Q file exists, move it to the new folder
  if (file.exists(qfile_path)) {
    file.copy(qfile_path, "./Maethio_Qfiles", overwrite = TRUE)
  }
}

# After moving, list the files in the new folder
sfiles <- list.files(path=("./Maethio_Qfiles"), full.names=T)
slist <- readQ(files=sfiles)

# Define a palette of colors to be used for different numbers of clusters
cluster_colors <- c("#51692d", "#f1c039", "#ff6347", "#4682b4", "#8a2be2", 
                    "#ff69b4", "#32cd32", "#ffa500", "#1e90ff", "#ff4500", 
                    "#20b2aa", "#dc143c", "#00ff7f", "#ff00ff", "#d2691e")

# Loop through each element in slist to plot and save
for (i in seq_along(slist)) {
  
  # Define the number of clusters (this is the number of columns in the Q matrix for slist[[i]])
  num_clusters <- ncol(slist[[i]])
  
  
  # Dynamically select the appropriate number of colors from the color palette
  selected_colors <- cluster_colors[1:num_clusters]
  
  # Define the output PDF file for each cluster plot
  output_pdf <- paste0("./cluster_plot_K", i + 1, ".pdf")
  
  # Plot the Q matrix for the current cluster solution
  plotQ(qlist = slist[i], 
        imgtype = "pdf", 
        height = 4.5,
        width = 10,
        clustercol = selected_colors, 
        dpi = 1200, 
        exportpath = "./")
  
  # Print a message indicating completion of the plot
  print(paste("Plot for K =", i + 1, "with", num_clusters, "clusters saved to", output_pdf))
}

########## Combined grp labeled #################################
# Extract the names of the data frames from the list
list_names <- names(slist)

grps1 <- read.delim("pop_map.txt", header=F,stringsAsFactors=F)
grplab_df <- data.frame("grp" = grps1$V2)


#group_order <- c('NOR', 'HAM', 'LIN', 'INV', 'IRE')
#grplab_df$grp <- factor(grplab_df$grp, levels = group_order)
#grplab_df <- grplab_df %>% arrange(grp)
#grplab_df$grp <- as.character(grplab_df$grp)

# Extract the numbers after "r" in the names
numbers <- as.numeric(str_match(list_names, "r[0-9]+\\.([0-9]+)\\.Q")[,2])

# Reorder the list based on the extracted numbers
ordered_list <- slist[order(numbers)]


plotQ(qlist = ordered_list[1:5], grplab = grplab_df, imgoutput = "join", imgtype = "jpeg", 
      height = 4.5,
      width = 20,   
      grplabsize = 3, 
      splabsize = 10,
      divcol = "white",      # Set the group separator line color to black
      divsize = 1,           # Increase the thickness of the separator lines
      divtype = "solid",
      exportpath = getwd())
####################################################################

############################## DAPC ###############################
install.packages("rnaturalearth")
library(rnaturalearth) #for mapping
library(dplyr)

Maethio_genlight <- vcfR2genlight(snp_vcf2)
pop.data2 = read.table("pop_map.txt", header = F)

num_clust <- find.clusters(Maethio_genlight)

Maethio_dapc <- dapc(Maethio_genlight, num_clust$grp)
scatter(Maethio_dapc, posi.da="bottomleft")


#tidy the data for plotting. 
dapc_data_df <-
  # as_tibble() converts the ind.coord matrix to a special data frame called a tibble. 
  as_tibble(Maethio_dapc$ind.coord, rownames = "individual") %>%
  # mutate changes or adds columns to a data frame. Here, we're adding the population and group assignment columns to the data frame
  mutate(population = pop.data2$V2,
         group = Maethio_dapc$grp)



#### Order the group
group_order <- c('MAN', 'HAM', 'LIN', 'DUN', 'IRE')

dapc_data_df$population <- factor(dapc_data_df$population, levels = group_order)

dapc_data_df <- dapc_data_df %>% arrange(population)

# Extract variance explained by LDs
ld_variance <- Maethio_dapc$eig / sum(Maethio_dapc$eig) * 100

# Convert to formatted strings for the x and y axis labels
x_label_dapc <- paste0("LD1 (", round(ld_variance[1], 2), "%)")
y_label_dapc <- paste0("LD2 (", round(ld_variance[2], 2), "%)")


# Create the ggplot with custom axis labels
dapc_plot <- ggplot(dapc_data_df, aes(
  x = LD1,
  y = LD2,
  colour = population
)) +
  geom_point(size = 3) +
  scale_fill_viridis_d(direction = -1) + 
  scale_color_manual(values = c("#f1c039", "#f37d21", "#51692d", "#56ba32", "#487aa1",  "#e41a1c", "#02ccf5")) + 
  labs(x = x_label_dapc, y = y_label_dapc) + 
  theme_classic() +
  ggtitle("DAPC Plot") +
  theme(plot.title = element_text(hjust = 0.5))



dapc_plot

