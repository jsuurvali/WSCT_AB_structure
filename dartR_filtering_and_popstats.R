library(dartR)
library(strataG)
library(tidyverse)

# read in the gzipped VCFs in genlight format and remove sites deviating from HWE in any population  
dat.okSNPs <- gl.read.vcf("WSCT.okSNPs.miss25.mac10.thin1k.vcf.gz", ind.metafile = "ind.metadata.okSNPs") %>%
  dartR::gl.filter.hwe(n.pop.threshold = 1, multi_comp = TRUE, multi_comp_method = "fdr") %>% gl.recalc.metrics

# extract ordered population lists
poplist.okSNPs <- dat.okSNPs$pop %>% unique %>% as.factor


# identify sites in LD with another site in any population, within 100kb of it
glmap.okSNPs <- dartR::gl.report.ld.map(dat.okSNPs, ld_max_pairwise = 100000) # none found


# calculate mean heterozygosity and Fis
gl.report.heterozygosity(dat.okSNPs)[poplist.okSNPs,] %>%
  write.table("WSCT.okSNPs.hetfis.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# calculate pairwise differentiation measures Fst and Dest
# also generates plots that can be saved directly from the graphical device
dist.okSNPs <- gl.report.fstat(dat.okSNPs)

dist.okSNPs$Stat_matrices$Fstp[poplist.okSNPs, poplist.okSNPs] %>%
  write.csv("WSCT.okSNPs.fst.csv") 

dist.okSNPs$Stat_matrices$Dest[poplist.okSNPs, poplist.okSNPs] %>% 
  write.csv("WSCT.okSNPs.Dest.csv") 

# plot FST and Dest as heatmaps
gl.report.fstat(dat.okSNPs, plot.stat = "Fstp") # WSCT.okSNPs.fst.heatmap
gl.report.fstat(dat.okSNPs, plot.stat = "Dest") # WSCT.okSNPs.Dest.heatmap

# Use LD to estimate Ne
Ne.okSNPs <- dat.okSNPs %>% genlight2gtypes %>% ldNe(maf.threshold=0.05, drop.missing = TRUE)
row.names(Ne.okSNPs) <- Ne.okSNPs$stratum
Ne.okSNPs[poplist.okSNPs,] %>% write.table("WSCT.okSNPs.Ne.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# ~~~ GENERATE PLOTS TO BE SAVED DIRECTLY FROM THE GRAPHICAL DEVICE ~~~

# identity by descent
descent.okSNPs <- gl.grm(dat.okSNPs) # WSCT.okSNPs.descent.heatmap

# genetic vs geographic distance + results of a Mantel Test
gl.ibd(dat.okSNPs, distance = 'Fst') # WSCT.okSNPs.Fst.mantel
gl.ibd(dat.okSNPs, distance = 'D') # WSCT.okSNPs.Dest.mantel
