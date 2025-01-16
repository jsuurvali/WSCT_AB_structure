setwd("E:/Dropbox/Alberta_Trout/trout/GATK_HF.DP7.GQ30.AB_filt.AEPonly.snps.indels/reduced2")
library("LEA")
library("tidyverse")

#uncomment the below 6 rows to run snmf with 100 repetitions and K values 1 to 30

#project.all = snmf("wsct.filtered.all.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")
#project.Job = snmf("wsct.filtered.Job.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")
#project.Bow = snmf("wsct.filtered.Bow.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")
#project.Oldman.OK = snmf("wsct.filtered.Oldman.OK.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")
#project.Oldman.notOK = snmf("wsct.filtered.Oldman.notOK.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")
#project.Castle = snmf("wsct.filtered.Castle.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")
#project.LvnRhCrowsnest = snmf("wsct.filtered.LvnRhCrowsnest.geno", K=1:30, entropy = TRUE, repetitions = 100, project = "new")


# read in previously created snmf projects (from the above step)
project.all <- load.snmfProject("wsct.filtered.all.snmfProject")
project.Job <- load.snmfProject("wsct.filtered.Job.snmfProject")
project.Bow <- load.snmfProject("wsct.filtered.Bow.snmfProject")
project.Oldman.OK <- load.snmfProject("wsct.filtered.Oldman.OK.snmfProject")
project.Oldman.notOK <- load.snmfProject("wsct.filtered.Oldman.notOK.snmfProject")
project.Castle <- load.snmfProject("wsct.filtered.Castle.snmfProject")
project.LvnRhCrowsnest <- load.snmfProject("wsct.filtered.LvnRhCrowsnest.snmfProject")

# create empty lists, one for each project
ce.all = list(); ce.Job = list(); ce.Bow = list()
ce.Oldman.notOK = list(); ce.Oldman.OK = list()
ce.Castle = list(); ce.LvnRhCrowsnest = list()

# populate lists with entropy data from the project files
for (k in 1:30){
  ce.all[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.all, K = k))
  ce.Bow[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.Bow, K = k))
  ce.Job[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.Job, K = k))
  ce.Oldman.OK[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.Oldman.OK, K = k))
  ce.Oldman.notOK[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.Oldman.notOK, K = k))
  ce.Castle[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.Castle, K = k))
  ce.LvnRhCrowsnest[[k]] <- cbind(rep(k, 100), seq(1:100), cross.entropy(project.LvnRhCrowsnest, K = k))
  }

ce.all <- data.frame(do.call("rbind", ce.all)) %>% setNames(c("K", "run", "entropy"))
ce.Job <- data.frame(do.call("rbind", ce.Job)) %>% setNames(c("K", "run", "entropy"))
ce.Bow <- data.frame(do.call("rbind", ce.Bow)) %>% setNames(c("K", "run", "entropy"))
ce.Oldman.OK <- data.frame(do.call("rbind", ce.Oldman.OK)) %>% setNames(c("K", "run", "entropy"))
ce.Oldman.notOK <- data.frame(do.call("rbind", ce.Oldman.notOK)) %>% setNames(c("K", "run", "entropy"))
ce.Castle <- data.frame(do.call("rbind", ce.Castle)) %>% setNames(c("K", "run", "entropy"))
ce.LvnRhCrowsnest <- data.frame(do.call("rbind", ce.LvnRhCrowsnest)) %>% setNames(c("K", "run", "entropy"))


# plot all cross entropy results together as distributions on a single plot
pdf("wsct.Xentropy.pdf", "width" = 6, "height" = 8)
ggplot(data = ce.all, aes(x = factor(K), y = entropy)) +
  geom_violin(data = ce.Bow, fill = "tan", draw_quantiles = 0.5) +
  geom_violin(data = ce.Job, fill = "skyblue", draw_quantiles = 0.5) +
  geom_violin(data = ce.Castle, fill = "orchid", draw_quantiles = 0.5) +
  geom_violin(data = ce.LvnRhCrowsnest, fill = "darkslategrey", draw_quantiles = 0.5) +
  geom_violin(data = ce.Oldman.notOK, fill = "pink", draw_quantiles = 0.5) +
  geom_violin(data = ce.Oldman.OK, fill = "grey", draw_quantiles = 0.5) +
  geom_violin(data = ce.all, fill = "black", draw_quantiles = 0.5) +
  xlab("Number of ancestral populations") + ylab("Cross-entropy") + theme_bw() + theme(legend.position="none")
dev.off()


# define function for plotting the minimum entropy Q matrices
Kplot <- function(dataset, samples, K_value, col_vector){
  indvdat <- read.table(samples)[,1]
  which.min(cross.entropy(dataset, K = K_value)) %>%
    Q(dataset, K_value, .) %>% data.frame %>%
    cbind("Sample" = indvdat) %>%
    gather("Ancestry", "Proportion", -Sample) %>%
    mutate(Sample = factor(Sample, levels = indvdat)) %>%
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

# in the next step, the defined colors on each row are exchangable. You need one color for each ancestry present
# the order of how ggplot assigns the colors varies by run, which is why I have a defined different order for each


# define colors for the "all" plots
cols.all <- list()
cols.all[[1]] <- c("brown")
cols.all[[2]] <- c("brown", "skyblue")
cols.all[[3]] <- c("brown", "lemonchiffon", "skyblue")
cols.all[[4]] <- c("tan", "brown", "lemonchiffon", "skyblue")
cols.all[[5]] <- c("brown", "skyblue", "lemonchiffon", "thistle", "tan")
cols.all[[6]] <- c("thistle", "navy", "brown", "tan", "skyblue", "lemonchiffon")
cols.all[[7]] <- c("skyblue", "thistle", "tan", "brown", "navy", "darkorchid4", "lemonchiffon")
cols.all[[8]] <- c("brown", "thistle", "darkorchid4", "darkseagreen3", "navy", "skyblue", "tan", "lemonchiffon")
cols.all[[9]] <- c("brown", "darkseagreen3", "thistle", "navy", "skyblue", "darkorchid4", "tan", "rosybrown", "lemonchiffon")
cols.all[[10]] <- c("thistle", "tan", "brown", "darkseagreen3", "lemonchiffon", "darkorchid4", "skyblue", "rosybrown","navy", "limegreen")
cols.all[[11]] <- c("darkslategrey", "lemonchiffon", "brown", "navy", "thistle", "skyblue", "yellow", "tan","rosybrown", "darkorchid4", "darkseagreen3")
cols.all[[12]] <- c("skyblue", "tan", "rosybrown", "darkorchid4", "darkseagreen3", "lemonchiffon", "yellow", "thistle","navy", "brown", "darkslategrey", "limegreen")
cols.all[[13]] <- c("darkorchid4", "brown", "navy", "limegreen", "yellow", "tan", "lightcyan", "skyblue","darkslategrey", "rosybrown", "thistle", "lemonchiffon", "darkseagreen3")
cols.all[[14]] <- c("darkslategrey", "brown", "darkseagreen3", "tan", "navy", "rosybrown", "skyblue", "darkgrey","lemonchiffon", "darkorchid4", "steelblue", "thistle", "limegreen", "yellow")
cols.all[[15]] <- c("darkgrey", "limegreen", "thistle", "tan", "lemonchiffon", "rosybrown", "darkslategrey", "lightcyan","seashell", "skyblue", "navy", "darkseagreen3", "yellow", "brown", "darkorchid4")
cols.all[[16]] <- c("brown", "darkgrey", "navy", "seashell", "yellow", "darkslategrey", "lightcyan", "rosybrown","thistle", "darkseagreen3", "lemonchiffon", "steelblue", "tan", "limegreen", "darkorchid4", "skyblue")
cols.all[[17]] <- c("darkseagreen3", "lightcyan", "olivedrab", "rosybrown", "seashell", "navy", "limegreen", "orange","darkorchid4", "darkslategrey", "lemonchiffon", "skyblue", "brown", "darkgrey", "thistle", "yellow", "tan")
cols.all[[18]] <- c("olivedrab", "darkseagreen3", "brown", "darkgrey", "navy", "skyblue", "rosybrown", "limegreen","darkorchid4", "tan", "lightcyan", "lemonchiffon", "orange", "darkslategrey", "thistle", "steelblue", "yellow", "seashell")
cols.all[[19]] <- c("rosybrown", "navy", "darkgrey", "yellow", "palegreen", "seashell", "darkorchid4", "darkseagreen3","brown", "limegreen", "tan", "darkslategrey", "steelblue", "lightcyan", "thistle", "skyblue", "olivedrab", "black", "lemonchiffon")
cols.all[[20]] <- c("yellow", "thistle", "skyblue", "seashell", "brown", "navy", "olivedrab", "limegreen","steelblue", "darkslategrey", "red", "orange", "rosybrown", "lightcyan", "darkorchid4", "darkgrey", "darkseagreen3", "tan", "lemonchiffon", "black")
cols.all[[21]] <- c("lemonchiffon", "olivedrab", "limegreen", "rosybrown", "palegreen", "red", "tan", "lightcyan","brown", "steelblue", "darkslategrey", "thistle", "darkgrey", "orange", "yellow", "lightcoral", "darkorchid4", "navy", "seashell", "skyblue", "darkseagreen3")
cols.all[[22]] <- c("skyblue", "seashell", "darkorchid4", "olivedrab", "tan", "thistle", "darkseagreen3", "darkorchid1", "yellow", "orange", "steelblue", "brown", "limegreen", "palegreen", "rosybrown",  "black", "lightcyan", "navy", "red", "darkgrey", "darkslategrey", "lemonchiffon")
cols.all[[23]] <- c("hotpink", "darkorchid4", "lemonchiffon", "#6F4E37", "palegreen", "red", "lightcyan", "skyblue", "seashell", "darkgrey", "orange", "limegreen", "black", "navy", "thistle", "steelblue", "tan", "brown", "yellow", "olivedrab", "rosybrown", "darkslategrey", "darkseagreen3")
cols.all[[24]] <- c("rosybrown", "steelblue", "darkgoldenrod", "skyblue", "lightcyan", "lemonchiffon", "darkseagreen3", "seashell", "navy", "palegreen", "darkgrey", "limegreen", "darkorchid4", "thistle", "lightcoral", "olivedrab", "yellow", "hotpink", "tan", "red", "brown", "darkslategrey", "black", "orange")
cols.all[[25]] <- c("steelblue", "rosybrown", "brown", "saddlebrown", "lightcyan", "darkgrey", "tan", "red", "black", "skyblue", "forestgreen", "palegreen", "lemonchiffon", "yellow", "thistle", "navy", "limegreen", "hotpink", "olivedrab", "#6F4E37", "darkseagreen3", "seashell", "darkslategrey", "orange", "darkorchid4")
cols.all[[26]] <- c("palegreen", "steelblue", "darkslategrey", "darkgrey", "lemonchiffon", "darkorchid4", "brown", "yellow", "black", "lightcyan", "seashell", "tan", "olivedrab", "lightcoral", "navy", "rosybrown", "thistle", "#923582", "orange", "greenyellow", "saddlebrown", "hotpink", "limegreen", "red", "skyblue", "darkseagreen3")
cols.all[[27]] <- c("darkorchid4", "rosybrown", "hotpink", "brown", "#6F4E37", "darkorchid1", "thistle", "steelblue", "saddlebrown", "#923582", "orange", "darkgrey", "palegreen", "lightcyan", "red", "navy", "seashell", "black", "limegreen", "lightcoral", "lemonchiffon", "tan", "darkslategrey", "olivedrab", "skyblue", "yellow", "darkseagreen3")
cols.all[[28]] <- c("forestgreen", "lightcyan", "black", "rosybrown", "palegreen", "orange", "#923582", "darkorchid1", "limegreen", "brown", "saddlebrown", "skyblue", "steelblue", "darkseagreen3", "darkgrey", "lemonchiffon", "tan", "darkslategrey", "thistle", "yellow", "hotpink", "navy", "darkorchid4", "red", "olivedrab", "lightpink", "lightcoral", "seashell")
cols.all[[29]] <- c("#6F4E37", "seashell", "olivedrab", "lightcoral", "darkslategrey", "darkorchid1", "forestgreen", "darkgrey", "skyblue", "limegreen", "brown", "black", "saddlebrown", "navy", "darkorchid4", "steelblue", "yellow", "orange", "thistle", "lemonchiffon", "#923582", "tan", "hotpink", "rosybrown", "palegreen", "red", "lightcyan", "#199999", "darkseagreen3")
cols.all[[30]] <- c("darkorchid4", "rosybrown", "#199999", "lightcyan", "lightcoral", "darkslategrey", "#6F4E37", "yellow", "darkseagreen3", "#224069", "brown", "steelblue", "tan", "limegreen", "palegreen", "darkorchid1", "black", "thistle", "forestgreen", "darkgrey", "skyblue", "orange", "#923582", "saddlebrown", "lemonchiffon", "navy", "olivedrab", "red", "seashell", "darkgoldenrod")

# define colors for the Job / Marvel plots
cols.Job <- list()
cols.Job[[1]] <- c("skyblue")
cols.Job[[2]] <- c("palegreen", "skyblue")
cols.Job[[3]] <- c("skyblue", "forestgreen", "palegreen")

# define colors for the Bow plots
cols.Bow <- list()
cols.Bow[[1]] <- c("tan")
cols.Bow[[2]] <- c("darkorchid4", "tan")
cols.Bow[[3]] <- c("skyblue", "darkorchid4", "tan")
cols.Bow[[4]] <- c("skyblue", "tan", "limegreen", "darkorchid4")
cols.Bow[[5]] <- c("tan", "skyblue", "darkorchid4", "limegreen", "darkorchid1")
cols.Bow[[6]] <- c("limegreen", "skyblue", "darkorchid1", "tan", "darkorchid4", "hotpink")
cols.Bow[[7]] <- c("hotpink", "darkorchid1", "darkorchid4", "skyblue", "tan", "limegreen", "gold")
cols.Bow[[8]] <- c("darkorchid1", "limegreen", "skyblue", "beige", "gold", "darkorchid4", "tan", "hotpink")

# define colors for the "good" Oldman plots
cols.Oldman.OK <- list()
cols.Oldman.OK[[1]] <- c("brown")
cols.Oldman.OK[[2]] <- c("brown", "lemonchiffon")
cols.Oldman.OK[[3]] <- c("brown", "skyblue", "lemonchiffon")
cols.Oldman.OK[[4]] <- c("thistle", "skyblue", "lemonchiffon", "brown")
cols.Oldman.OK[[5]] <- c("navy", "brown", "skyblue", "thistle", "lemonchiffon")
cols.Oldman.OK[[6]] <- c("lemonchiffon", "navy", "thistle", "skyblue", "brown", "darkseagreen3")
cols.Oldman.OK[[7]] <- c("skyblue", "navy", "rosybrown", "thistle", "brown", "darkseagreen3", "lemonchiffon")
cols.Oldman.OK[[8]] <- c("thistle", "darkslategrey", "skyblue", "navy", "darkseagreen3", "brown", "rosybrown", "lemonchiffon")
cols.Oldman.OK[[9]] <- c("yellow", "lemonchiffon", "darkgrey", "darkslategrey", "navy", "darkseagreen3", "rosybrown", "skyblue", "thistle")
cols.Oldman.OK[[10]] <- c("darkslategrey", "rosybrown", "darkseagreen3", "darkgrey", "skyblue", "steelblue", "navy", "yellow", "lemonchiffon", "thistle")
cols.Oldman.OK[[11]] <- c("yellow", "rosybrown", "lightcyan", "lemonchiffon", "darkseagreen3", "thistle", "navy", "steelblue", "darkslategrey", "darkgrey", "skyblue")
cols.Oldman.OK[[12]] <- c("navy", "orange", "yellow", "lightcyan", "rosybrown", "darkslategrey", "steelblue", "darkgrey", "darkseagreen3", "lemonchiffon", "skyblue", "thistle")
cols.Oldman.OK[[13]] <- c("darkseagreen3", "skyblue", "yellow", "rosybrown", "lemonchiffon", "red", "orange", "thistle", "darkgrey", "steelblue", "navy", "darkslategrey", "lightcyan")
cols.Oldman.OK[[14]] <- c("orange", "darkslategrey", "navy", "darkgrey", "skyblue", "thistle", "steelblue", "darkseagreen3", "seashell", "yellow", "red", "lightcyan", "rosybrown", "lemonchiffon")
cols.Oldman.OK[[15]] <- c("darkseagreen3", "lemonchiffon", "seashell", "rosybrown", "orange", "navy", "lightcoral", "skyblue", "yellow", "red", "thistle", "lightcyan", "darkgrey", "steelblue", "darkslategrey")
cols.Oldman.OK[[16]] <- c("steelblue", "thistle", "orange", "red", "darkslategrey", "lemonchiffon", "skyblue", "lightcoral", "yellow", "darkgrey", "seashell", "navy", "rosybrown", "lightcyan", "black", "darkseagreen3")
cols.Oldman.OK[[17]] <- c("black", "thistle", "orange", "#6F4E37", "steelblue", "red", "navy", "darkseagreen3", "darkgrey", "rosybrown", "seashell", "skyblue", "lightcoral", "darkslategrey", "lightcyan", "lemonchiffon", "yellow")


# define colors for the "bad" Oldman plots
cols.Oldman.notOK <- list()
cols.Oldman.notOK[[1]] <- c("brown")
cols.Oldman.notOK[[2]] <- c("lemonchiffon", "brown")
cols.Oldman.notOK[[3]] <- c("darkgrey", "lemonchiffon", "brown")
cols.Oldman.notOK[[4]] <- c("brown", "olivedrab", "lemonchiffon", "darkgrey")
cols.Oldman.notOK[[5]] <- c("yellow", "olivedrab", "brown", "darkgrey", "lemonchiffon")
cols.Oldman.notOK[[6]] <- c("brown", "darkgrey", "lemonchiffon", "olivedrab", "yellow", "blue")

# define colors for the Castle plots
cols.Castle <- list()
cols.Castle[[1]] <- c("thistle")
cols.Castle[[2]] <- c("thistle", "skyblue")
cols.Castle[[3]] <- c("navy", "thistle", "skyblue")
cols.Castle[[4]] <- c("navy", "darkseagreen3", "thistle", "skyblue")
cols.Castle[[5]] <- c("darkseagreen3", "thistle", "navy", "skyblue", "steelblue")

# define colors for the Livingstone / N Racehorse / Crowsnest plots
cols.LvnRhCrowsnest <- list()
cols.LvnRhCrowsnest[[1]] <- c("darkgrey")
cols.LvnRhCrowsnest[[2]] <- c("darkgrey", "darkslategrey")
cols.LvnRhCrowsnest[[3]] <- c("darkslategrey", "yellow", "darkgrey")
cols.LvnRhCrowsnest[[4]] <- c("darkslategrey", "darkgrey", "yellow", "lightcyan")
cols.LvnRhCrowsnest[[5]] <- c("yellow", "lightcyan", "red", "darkslategrey", "darkgrey")
cols.LvnRhCrowsnest[[6]] <- c("lightcyan", "red", "orange", "darkslategrey", "yellow", "darkgrey")
cols.LvnRhCrowsnest[[7]] <- c("yellow", "lightcyan", "red", "seashell", "darkslategrey", "orange", "darkgrey")
cols.LvnRhCrowsnest[[8]] <- c("salmon", "seashell", "yellow", "red", "lightcyan", "darkslategrey", "orange", "darkgrey")
cols.LvnRhCrowsnest[[9]] <- c("salmon", "darkslategrey", "red", "#6F4E37", "lightcyan", "yellow", "darkgrey", "seashell", "orange")


# create the plots by using the defined function and colors
plots.all <- list()
for (i in 1:length(cols.all)){plots.all[[i]] <- Kplot(project.all, "keeplist.miss25", i, cols.all[[i]])}

plots.Job <- list()
for (i in 1:length(cols.Job)){plots.Job[[i]] <- Kplot(project.Job, "keeplist.miss25.Job", i, cols.Job[[i]])}

plots.Bow <- list()
for (i in 1:length(cols.Bow)){plots.Bow[[i]] <- Kplot(project.Bow, "keeplist.miss25.Bow", i, cols.Bow[[i]])}

plots.Oldman.OK <- list()
for (i in 1:length(cols.Oldman.OK)){plots.Oldman.OK[[i]] <- Kplot(project.Oldman.OK, "keeplist.miss25.Oldman.OK", i, cols.Oldman.OK[[i]])}

plots.Oldman.notOK <- list()
for (i in 1:length(cols.Oldman.notOK)){plots.Oldman.notOK[[i]] <- Kplot(project.Oldman.notOK, "keeplist.miss25.Oldman.notOK", i, cols.Oldman.notOK[[i]])}

plots.Castle <- list()
for (i in 1:length(cols.Castle)){plots.Castle[[i]] <- Kplot(project.Castle, "keeplist.miss25.Castle", i, cols.Castle[[i]])}

plots.LvnRhCrowsnest <- list()
for (i in 1:length(cols.LvnRhCrowsnest)){plots.LvnRhCrowsnest[[i]] <- Kplot(project.LvnRhCrowsnest, "keeplist.miss25.LvnRhCrowsnest", i, cols.LvnRhCrowsnest[[i]])}

# save plots as pdf files with 1 plot per page
pdf("wsct.snmf.all.pdf", height = 2, width = 41)
plots.all
dev.off()

pdf("wsct.snmf.Job.pdf", height = 2, width = 15)
plots.Job
dev.off()

pdf("wsct.snmf.Bow.pdf", height = 2, width = 4)
plots.Bow
dev.off()

pdf("wsct.snmf.Oldman.OK.pdf", height = 2, width = 16)
plots.Oldman.OK
dev.off()

pdf("wsct.snmf.Oldman.notOK.pdf", height = 2, width = 10)
plots.Oldman.notOK
dev.off()

pdf("wsct.snmf.Castle.pdf", height = 2, width = 7)
plots.Castle
dev.off()

pdf("wsct.snmf.LvnRhCrowsnest.pdf", height = 2, width = 5)
plots.LvnRhCrowsnest
dev.off()