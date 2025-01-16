library(dartR)
library(tidyverse)

# read in the population assignment file, in the same order as individuals in the VCF
pops <- read.table("sample30.txt", header = FALSE, sep = "\t", col.names = c("sample", "population"))

# read in the gzipped VCF, convert to genind format and add population names
dat1 <- gl.read.vcf("forfilter.vcf.gz")
dat1$pop <- factor(pops$population, levels = unique(pops$population))

# remove sites deviating from Hardy Weinberg equilibrium in any of the populations
dat2 <- dat1 %>% gl.filter.hwe(
  n.pop.threshold = 1,
  multi_comp = TRUE,
  multi_comp_method = "fdr"
)

# remove sites that in at least 2 population are in significant LD with another site within 100kb
glmap <- gl.report.ld.map(dat2, ld_max_pairwise = 100000)
dat3 <- dat2 %>% gl.filter.ld(ld_report = glmap,
                              threshold = 0.2,
                              pop.limit = 2
                              )

# extract the SNPs from the output and output as the file "filtered.txt"
as.character(dat3$chromosome) %>% cbind(dat3$position) %>% data.frame %>%
  write.table(file = "filtered.txt",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              sep = "\t"
              )

