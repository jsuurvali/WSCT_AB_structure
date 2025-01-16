library("tidyverse")
library("vcfR")
library("adegenet")
library("adegraphics")

# read in data
pcadat <- read.vcfR("wsct.filtered.Job.vcf.gz") %>% vcfR2genind()

# define function to replace missing values with uninformative data points (requires center = TRUE later)
missing_to_mean <- function(x) {x[is.na(x)] <- mean(x, na.rm = TRUE); return(x)}

# calculate principal components
res.pca <- apply(pcadat$tab,2,missing_to_mean) %>% dudi.pca(center = TRUE, scannf = FALSE, nf = 3)
res.pca$eig <- sapply(res.pca$eig, function(x) round(x/sum(res.pca$eig)*100, digits = 2)) # make the eigenvectors show percentage

# write scree plots to file
pdf("wsct.filtered.Job.screeplot.pdf", height = 3, width = 3)
screeplot(res.pca, main = "all\n", xlab = "\nPrincipal Component", ylab = "Explained data (%)")
dev.off()

# pca for axes 1 and 2
pops <- c(
  rep("Job, Marvel, Job-Like", 689),
  "CasW3_17_44",
  rep("Job, Marvel, Job-Like", 8)
) %>% as.factor()


pdf("wsct.filtered.Job.pca12.pdf")
s.class(res.pca$li[,c(1,2)],
        fac = pops,
        xlab = paste0("\nPCA 1 (", res.pca$eig[1], "%)"),
        ylab = paste0("PCA 2 (", res.pca$eig[2], "%)\n"),
        main = "Job\n")
dev.off()

# define populations and plot PCA for axes 2 and 3
pops <- c(
  rep("Job", 460),
  rep("Marvel", 28),
  rep("Other Job-like", 66),
  rep("Baril 1", 32),
  rep("Other Job-like", 33),
  rep("Wilkinson", 34),
  rep("Other Job-like", 36),
  "CasW3_17_44",
  rep("Other Job-like", 8)
) %>% as.factor()


pdf("wsct.filtered.Job.pca23.pdf")
s.class(res.pca$li[,c(2,3)],
        fac = pops,
        xlab = paste0("\nPCA 2 (", res.pca$eig[2], "%)"),
        ylab = paste0("PCA 3 (", res.pca$eig[3], "%)\n"),
        col = c("limegreen", "black", "darkblue", "skyblue", "steelblue", "darkgreen"),
        main = "Job\n")
dev.off()
