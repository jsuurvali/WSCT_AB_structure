library(LEA)
library(tidyverse)

project.okSNPs <- load.snmfProject("WSCT.okSNPs.miss25.mac10.thin1k.snmfProject")
annot.okSNPs <- read.table("WSCT.okSNPs.annot.tsv", header = TRUE)[,c("site", "lat", "lon")]

K17.okSNPs <- which.min(cross.entropy(project.okSNPs, K = 17)) %>%
  Q(project.okSNPs, 17, .) %>% data.frame %>%
  cbind(annot.okSNPs, .) %>% group_by(site, lat, lon) %>% summarise_at(vars(V1:V17), mean)
write.table(K17.okSNPs, "WSCT.okSNPs.admixture.17.piedata.csv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)






