library("tidyverse")
library("vcfR")
library("adegenet")
library("adegraphics")

# read in data
pcadat <- read.vcfR("wsct.filtered.Oldman.OK.vcf.gz") %>% vcfR2genind()

# define function to replace missing values with uninformative data points (requires center = TRUE later)
missing_to_mean <- function(x) {x[is.na(x)] <- mean(x, na.rm = TRUE); return(x)}

# calculate principal components
res.pca <- apply(pcadat$tab,2,missing_to_mean) %>% dudi.pca(center = TRUE, scannf = FALSE, nf = 4)
res.pca$eig <- sapply(res.pca$eig, function(x) round(x/sum(res.pca$eig)*100, digits = 2)) # make the eigenvectors show percentage

# write scree plots to file
pdf("wsct.filtered.Oldman.OK.screeplot.pdf", height = 3, width = 3)
screeplot(res.pca, main = "all\n", xlab = "\nPrincipal Component", ylab = "Explained data (%)")
dev.off()


# define populations for plotting
pops <- c(
  rep("Willow", 15),
  rep("Beaver", 5),
  rep("Callum", 41),
  rep("Uppermost", 147),
  rep("Other Oldman", 197),
  rep("Crowsnest", 16),
  rep("Star", 26),
  rep("Lynx", 39),
  rep("O'Haggen", 48),
  rep("Castle", 109),
  rep("Job-like", 45),
  rep("Castle", 77)
) %>% as.factor()

# pca axes 1 and 4
pdf("wsct.filtered.Oldman.OK.pca14.pdf")
s.class(res.pca$li[,c(1,4)],
        fac = pops,
        xlab = paste0("\nPCA 1 (", res.pca$eig[1], "%)"),
        ylab = paste0("PCA 4 (", res.pca$eig[4], "%)\n"),
        col = c("saddlebrown", "rosybrown", "thistle", "darkslategrey", "slateblue", "navy", "darkseagreen", "brown",
                "skyblue", "orange", "black"),
        main = "Oldman, distinct\n")
dev.off()

# pca axes 2 and 3
pdf("wsct.filtered.Oldman.OK.pca23.pdf")
s.class(res.pca$li[,c(2,3)],
        fac = pops,
        xlab = paste0("\nPCA 2 (", res.pca$eig[2], "%)"),
        ylab = paste0("PCA 3 (", res.pca$eig[3], "%)\n"),
        col = c("saddlebrown", "rosybrown", "thistle", "darkslategrey", "slateblue", "navy", "darkseagreen", "brown",
                "skyblue", "orange", "black"),
        main = "Oldman, distinct\n")
dev.off()
