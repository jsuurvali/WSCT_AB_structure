library(dartR)
library(strataG)
library(tidyverse)

# read in data
dat.annot <- read.table("annot.tsv", header = TRUE, sep = "\t")
dat <- gl.read.vcf("wsct.filtered.all.vcf.gz")

# get each site-date combo as a separate "population"
# discard sites with < 5 individuals
annot.site <- dat.annot$siteID %>%
  table %>% data.frame %>% setNames(c("pop", "Freq")) %>% filter(Freq >= 5)
annot.site <- dat.annot[dat.annot$siteID %in% annot.site$pop,]
blacklist.site <- dat.annot[! dat.annot$sample %in% annot.site$sample,1]
dat.site <- dat %>% gl.drop.ind(blacklist.site, recalc = TRUE)
dat.site$pop <- factor(annot.site$siteID, levels = unique(annot.site$siteID))

# get each population
# discard populations with < 5 individuals
annot.pops <- dat.annot$population %>%
  table %>% data.frame %>% setNames(c("pop", "Freq")) %>% filter(Freq >= 5)
annot.pops <- dat.annot[dat.annot$population %in% annot.pops$pop,]
blacklist.pops <- dat.annot[! dat.annot$sample %in% annot.pops$sample,1]
dat.pops <- dat %>% gl.drop.ind(blacklist.pops, recalc = TRUE)
dat.pops$pop <- factor(annot.pops$population, levels = unique(annot.pops$population))

# calculate mean heterozygosity and Fis
gl.report.heterozygosity(dat.pops) %>% write.table("hetfis.pops.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
gl.report.heterozygosity(dat.site) %>% write.table("hetfis.site.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# calculate pairwise Fst
FST.pops <- gl.fst.pop(dat.pops) %>% write.table("fst.pops.tsv", sep = "\t", quote = FALSE)
FST.site <- gl.fst.pop(dat.site) %>% write.table("fst.site.tsv", sep = "\t", quote = FALSE)

# impute missing genotypes for each site and use LD to estimate Ne
Ne.site <- dat.site %>% gl.impute %>% genlight2gtypes %>% ldNe(maf.threshold=0.05)
rownames(Ne.site) <- Ne.site$stratum
order.site <- annot.site$siteID %>% as.vector %>% unique
Ne.site[order.site,] %>% write.table("Ne.site.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# impute missing genotypes for each population and use LD to estimate Ne
Ne.pops <- dat.pops %>% gl.impute %>% genlight2gtypes %>% ldNe(maf.threshold=0.05)
rownames(Ne.pops) <- Ne.pops$stratum
order.pops <- annot.pops$population %>% as.vector %>% unique
Ne.pops[order.pops,] %>% write.table("Ne.pops.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
