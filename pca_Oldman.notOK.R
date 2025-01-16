library("tidyverse")
library("vcfR")
library("adegenet")
library("adegraphics")

# read in data
pcadat <- read.vcfR("wsct.filtered.Oldman.notOK.vcf.gz") %>% vcfR2genind()

# define function to replace missing values with uninformative data points (requires center = TRUE later)
missing_to_mean <- function(x) {x[is.na(x)] <- mean(x, na.rm = TRUE); return(x)}

# calculate principal components
res.pca <- apply(pcadat$tab,2,missing_to_mean) %>% dudi.pca(center = TRUE, scannf = FALSE, nf = 2)
res.pca$eig <- sapply(res.pca$eig, function(x) round(x/sum(res.pca$eig)*100, digits = 2)) # make the eigenvectors show percentage

# write scree plots to file
pdf("wsct.filtered.Oldman.notOK.screeplot.pdf", height = 3, width = 3)
screeplot(res.pca, main = "all\n", xlab = "\nPrincipal Component", ylab = "Explained data (%)")
dev.off()

# define populations for plotting
pops <- c(
  rep("Upper Oldman / Honeymoon / Hidden", 132),
  rep("Livingstone", 75),
  rep("Dutch, Fly, Ernest, Todd", 66),
  rep("Racehorse", 180),
  rep("Dutch, Fly, Ernest, Todd", 11)
) %>% as.factor()


# plot PCA axes 1 and 2

pdf("wsct.filtered.Oldman.notOK.pca12.pdf")
s.class(res.pca$li[,c(1,2)],
        fac = pops,
        xlab = paste0("\nPCA 1 (", res.pca$eig[1], "%)"),
        ylab = paste0("PCA 2 (", res.pca$eig[2], "%)\n"),
        col = c("brown", "darkgrey", "orange", "olivedrab"),
        main = "Oldman, mixed\n")
dev.off()
