packages <- c("vcfR","ggplot2", "dplyr", "tidyr", "reshape", "StAMPP", "adegenet", "corrplot")

# Check if packages are installed
installed_packages <- packages %in% rownames(installed.packages())

# Install packages that are not installed
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}
library(tidyr)
library(vcfR)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(corrplot)
library(reshape)
args = commandArgs(trailingOnly=TRUE)

########### Create stats directory #################
dir.create("03_Analysis/05_stats", showWarnings = T)
dir.create("03_Analysis/05_stats/01_PCA", showWarnings = T)

##running PCA
snp_vcf2 = read.vcfR(args[1])
pop.data2 = read.table(args[2], header = F)
gl.snp2 <- vcfR2genlight(snp_vcf2)
pop(gl.snp2) <- rep(pop.data2$V2)
snp.pca2 <- glPca(gl.snp2, nf = 10)

##write PCA scores
snp.pca.scores2 <- as.data.frame(snp.pca2$scores)
snp.pca.scores2$pop <- pop(gl.snp2)
PCA_adegen_path <- file.path("03_Analysis/05_stats/01_PCA", paste(args[3], "adegenetPCA.txt", sep = "_"))
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
PCA_eigen_path <- file.path("03_Analysis/05_stats/01_PCA", paste(args[3], "eigen_summary.txt", sep = "_"))
write.csv(eigen,file=PCA_eigen_path,row.names=TRUE,quote=FALSE)

##PCA plotting
data2 = read.delim(PCA_adegen_path) #I have manually added population information to this file prior to loading it
mycol = c("#f1c039","#f37d21", "#51692d", "#56ba32")
ggplot(data2, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(size = 2) + 
  scale_color_manual(values=mycol) +
  theme_classic()+
  xlab("PC1 (1.99%)") +
  ylab("PC2 (1.34%)")

PCA_plot_output <- "03_Analysis/05_stats/01_PCA"
ggsave(path = PCA_plot_output, filename = paste(args[3], "PCAplot.png", sep = "_"), device='png', dpi=1000)
