library("tidyverse")
library("vcfR")
library("adegenet")
library("adegraphics")

# read in data
pcadat <- read.vcfR("wsct.filtered.Bow.vcf.gz") %>% vcfR2genind()

# define function to replace missing values with uninformative data points (requires center = TRUE later)
missing_to_mean <- function(x) {x[is.na(x)] <- mean(x, na.rm = TRUE); return(x)}

# calculate principal components
res.pca <- apply(pcadat$tab,2,missing_to_mean) %>% dudi.pca(center = TRUE, scannf = FALSE, nf = 4)
res.pca$eig <- sapply(res.pca$eig, function(x) round(x/sum(res.pca$eig)*100, digits = 2)) # make the eigenvectors show percentage

# write scree plots to file
pdf("wsct.filtered.Bow.screeplot.pdf", height = 3, width = 3)
screeplot(res.pca, main = "all\n", xlab = "\nPrincipal Component", ylab = "Explained data (%)")
dev.off()


# define populations and plot the results of PCA axes 1 and 2
pops <- c(
  rep("Marvel", 28),
  rep("Other Bow", 45),
  rep("Prairie, Trail", 18),
  rep("Other Bow", 22),
  rep("Picklejar 2", 18),
  rep("Picklejar 4", 30)
) %>% as.factor()


pdf("wsct.filtered.Bow.pca12.pdf")
s.class(res.pca$li[,c(1,2)],
        fac = pops,
        xlab = paste0("\nPCA 1 (", res.pca$eig[1], "%)"),
        ylab = paste0("PCA 2 (", res.pca$eig[2], "%)\n"),
        col = c("slateblue", "darkgoldenrod", "darkorchid", "darkorchid4", "forestgreen"),
        main = "Bow\n")
dev.off()


# define populations and plot the results of PCA axes 3 and 4
pops <- c(
  rep("Marvel", 28),
  rep("Other Bow", 45),
  rep("Prairie, Trail", 18),
  rep("Silvester", 12),
  rep("Other Bow", 2),
  rep("Highwood", 8),
  rep("Picklejar 2", 18),
  rep("Picklejar 4", 30)
) %>% as.factor()


pdf("wsct.filtered.Bow.pca34.pdf")
s.class(res.pca$li[,c(3,4)],
        fac = pops,
        xlab = paste0("\nPCA 3 (", res.pca$eig[3], "%)"),
        ylab = paste0("PCA 4 (", res.pca$eig[4], "%)\n"),
        col = c("gold", "slateblue", "darkgoldenrod", "darkorchid", "darkorchid4", "forestgreen", "magenta"),
        main = "Bow\n")
dev.off()
