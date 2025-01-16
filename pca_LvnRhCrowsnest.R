library("tidyverse")
library("vcfR")
library("adegenet")
library("adegraphics")

# read in data
pcadat <- read.vcfR("wsct.filtered.LvnRhCrowsnest.vcf.gz") %>% vcfR2genind()

# define function to replace missing values with uninformative data points (requires center = TRUE later)
missing_to_mean <- function(x) {x[is.na(x)] <- mean(x, na.rm = TRUE); return(x)}

# calculate principal components
res.pca <- apply(pcadat$tab,2,missing_to_mean) %>% dudi.pca(center = TRUE, scannf = FALSE, nf = 4)
res.pca$eig <- sapply(res.pca$eig, function(x) round(x/sum(res.pca$eig)*100, digits = 2)) # make the eigenvectors show percentage

# write scree plots to file
pdf("wsct.filtered.LvnRhCrowsnest.screeplot.pdf", height = 3, width = 3)
screeplot(res.pca, main = "all\n", xlab = "\nPrincipal Component", ylab = "Explained data (%)")
dev.off()

# define populations for plotting
pops <- c(
  rep("Livingstone: Upper", 28),
  rep("Livingstone: Ridge", 53),
  rep("Livingstone: Deep", 16),
  rep("Station", 21),
  rep("North Racehorse", 54),
  rep("Crowsnest: Gold", 25),
  rep("Crowsnest: Blairmore, Girardi, Rock", 16),
  rep("Crowsnest: Star", 26)
) %>% as.factor()

# PCA axes 1 and 4
pdf("wsct.filtered.LvnRhCrowsnest.pca14.pdf")
s.class(res.pca$li[,c(1,4)],
        fac = pops,
        xlab = paste0("\nPCA 1 (", res.pca$eig[1], "%)"),
        ylab = paste0("PCA 4 (", res.pca$eig[4], "%)\n"),
        col = c("darkslategrey", "red", "cyan", "orange", "purple", "darkgrey", "gold", "salmon"),
        main = "Livingstone, Racehorse, Station, Crowsnest\n")
dev.off()

# PCA axes 2 and 3
pdf("wsct.filtered.LvnRhCrowsnest.pca23.pdf")
s.class(res.pca$li[,c(2,3)],
        fac = pops,
        xlab = paste0("\nPCA 2 (", res.pca$eig[2], "%)"),
        ylab = paste0("PCA 3 (", res.pca$eig[3], "%)\n"),
        col = c("darkslategrey", "red", "cyan", "orange", "purple", "darkgrey", "gold", "salmon"),
        main = "Livingstone, Racehorse, Station, Crowsnest\n")
dev.off()
