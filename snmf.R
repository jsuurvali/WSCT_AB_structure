library(tidyverse)
library(LEA)


samples <- scan("keeplist.samples.2", character())[,1] # read in a list of sample names

project.new = snmf("WSCT.okSNPs.miss25.mac10.thin1k.geno", K=1:25, entropy = TRUE, repetitions = 20, project = "new")
project.new <- load.snmfProject("WSCT.okSNPs.miss25.mac10.thin1k.snmfProject")

ce.new = list()
for (k in 1:25){ce.new[[k]] <- cbind(rep(k, 20), seq(1:20), cross.entropy(project.new, K = k))}
ce.new <- data.frame(do.call("rbind", ce.new)) %>% setNames(c("K", "run", "entropy"))

# plot all cross entropy results together as distributions on a single plot
pdf("WSCT.okSNPs.Xentropy.pdf", "width" = 6, "height" = 8)
ggplot(data = ce.new, aes(x = factor(K), y = entropy)) +
  geom_violin(fill = "darkgrey", draw_quantiles = 0.5) +
  xlab("Number of ancestral populations") + ylab("Cross-entropy") + theme_bw() + theme(legend.position="none")
dev.off()


# define function for plotting the minimum entropy Q matrices
Kplot <- function(dataset, K_value, col_vector){
  which.min(cross.entropy(dataset, K = K_value)) %>%
    Q(dataset, K_value, .) %>% data.frame %>%
    cbind("Sample" = samples) %>%
    gather("Ancestry", "Proportion", -Sample) %>%
    mutate(Sample = factor(Sample, levels = samples)) %>%
    ggplot(aes(x = Sample, y = Proportion, fill = Ancestry)) +
    geom_bar(position = 'fill', stat = 'identity') +
    scale_fill_manual(values = col_vector) +
    theme_bw() + ylab(paste0("K=", K_value)) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8, angle=0, vjust = 0.5),
          axis.text.x = element_text(size = 2, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),
          axis.ticks.x = element_line(linewidth = 0.2),
          axis.ticks.y = element_blank(),
          legend.position="none")
}

cols <- c("yellow", "skyblue", "limegreen", "red", "darkseagreen3", "darkorchid4",
             "darkgrey", "thistle","lemonchiffon", "darkslategrey", "black", "orange",
             "gold", "hotpink", "rosybrown", "tan", "steelblue")

pdf("WSCT.okSNPs.admixture.17.pdf", width = 30, height = 2)
	Kplot(project.new, 17, cols)
dev.off()