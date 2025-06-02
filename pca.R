library(tidyverse)
library(vcfR)
library(adegenet)
library(adegraphics)

# read in data
pca.res.okSNPs <- read.vcfR("WSCT.okSNPs.miss25.mac10.thin1k.vcf.gz") %>% vcfR2genlight %>% glPca(nf = 8)

# define population assignments

pca.pops.okSNPs <- c(
  rep("Job", 33),
  rep("Marvel", 27),
  rep("Mockingbird", 50),
  rep("Job-like E", 27),
  rep("Trail", 44),
  rep("Silvester", 56),
  rep("Job-like H", 58),
  rep("Picklejar", 28),
  rep("Cutthroat", 45),
  rep("Westrup", 25),
  rep("Playle", 26),
  rep("U Oldman", 78),
  rep("N Racehorse", 59),
  rep("Ridge", 28),
  rep("Deep", 19),
  rep("Girardi", 36),
  rep("Gold", 50),
  rep("O'Haggen", 41),
  rep("Lost", 59),
  rep("Job-Like C", 28),
  rep("Castle S", 97)
) %>% as.factor()

# define color scheme
cols.pca <- c("thistle", "gold", "orange", "darkslategrey", "red", "lightblue", "cyan",
           "deepskyblue", "dodgerblue","royalblue", "skyblue", "tan", "yellow", "darkseagreen3",
           "darkorchid4", "rosybrown", "darkgrey", "hotpink", "limegreen", "lemonchiffon", "black")

# function to create the PCA plots
pcaplot <- function(filenamepdf, pcadata, popdata, colordata, axis1, axis2){
  pdf(filenamepdf)
  s.class(pcadata$scores[,c(axis1, axis2)],
          fac = popdata, col = colordata,
          xlab = paste0("\nPCA ", axis1, " (", round(pcadata$eig[axis1], 2), "%)"),
          ylab = paste0("PCA ", axis2, " (", round(pcadata$eig[axis2], 2), "%)\n"),
          )
  dev.off()
}

# plot the data
pcaplot("WSCT.okSNPs.pca12.pdf", pca.res.okSNPs, pca.pops.okSNPs, cols.pca, 1, 2)
pcaplot("WSCT.okSNPs.pca34.pdf", pca.res.okSNPs, pca.pops.okSNPs, cols.pca, 3, 4)
pcaplot("WSCT.okSNPs.pca56.pdf", pca.res.okSNPs, pca.pops.okSNPs, cols.pca, 5, 6)
pcaplot("WSCT.okSNPs.pca78.pdf", pca.res.okSNPs, pca.pops.okSNPs, cols.pca, 7, 8)