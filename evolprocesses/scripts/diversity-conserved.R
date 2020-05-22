rm(list=ls())
setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCanalysis/evolprocesses/mugsy_out/')

library("ape")
library("ggtree")
require(phylobase)

BGC <- c("rif")

tr <- read.tree(file = paste('RAxML_bipartitionsBranchLabels.', BGC, 'BGC', sep = ''))
#tr <- read.tree('../DTL_analysis/RAxML_bipartitionsBranchLabels.salinicore')
data <- read.table('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/genomeID.txt', header = T, sep = "\t", comment.char = "")
taxcolors <- setNames(as.character(data$clade_color), data$clade_name)

library(phytools)
tr2 <- midpoint.root(tr)
write.tree(tr2, file = paste("../DTL_analysis/raxml-", BGC, "BGC.root.tre", sep = ''))
write.tree(tr, file = paste("../DTL_analysis/raxml-", BGC, "BGC.tre", sep = ''))
#write.tree(tr2, file = paste("../DTL_analysis/RAxML_boot.core.root.tre", sep = ''))

library(ggplot2)
#pdf(file = "coregenome-raxml.pdf", height = 12, width = 4)
pdf(file = paste(BGC, "-raxml.pdf", sep = ""), height = 6, width = 4)
t <- ggtree(tr2, layout="rectangular") 
t %<+% data + 
  geom_tippoint(aes(colour = clade_name, size = 3, shape = clade_name)) +
  scale_color_manual(values = taxcolors) +
  scale_shape_manual(values = rep(15, length(unique(data$clade_name)))) +
  theme_tree2() +
  geom_treescale(width = 0.01)
dev.off()

# if subcluster is missing altogether then skips color for some reason







############################################################
############################################################
############################################################

# compute phylo distance of genome by phylo distance of BGCs
rm(list=ls())
setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCanalysis/evolprocesses/DTL_analysis/')
library(tidyverse)
library(ape)
coretree <- read.tree("RAxML_bipartitionsBranchLabels.salinicore")

# convert to distance matrix
PatristicDistMatrix <- cophenetic(coretree)

library(reshape2)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
phylodist <- B[!is.na(B$value),]
names(phylodist) <- c("genome1", "genome2", "phylodist")


lymtree <- read.tree("RAxML_bipartitionsBranchLabels.lymBGC")
PatristicDistMatrix <- cophenetic(lymtree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
lymdist <- B[!is.na(B$value),]
names(lymdist) <- c("genome1", "genome2", "BGCdist")
lymdist$BGC <- "lym"

lomtree <- read.tree("RAxML_bipartitionsBranchLabels.lomBGC")
PatristicDistMatrix <- cophenetic(lomtree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
lomdist <- B[!is.na(B$value),]
names(lomdist) <- c("genome1", "genome2", "BGCdist")
lomdist$BGC <- "lom"

riftree <- read.tree("RAxML_bipartitionsBranchLabels.rifBGC")
PatristicDistMatrix <- cophenetic(riftree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
rifdist <- B[!is.na(B$value),]
names(rifdist) <- c("genome1", "genome2", "BGCdist")
rifdist$BGC <- "rif"

saltree <- read.tree("RAxML_bipartitionsBranchLabels.salBGC")
PatristicDistMatrix <- cophenetic(saltree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
saldist <- B[!is.na(B$value),]
names(saldist) <- c("genome1", "genome2", "BGCdist")
saldist$BGC <- "sal"

spttree <- read.tree("RAxML_bipartitionsBranchLabels.sptBGC")
PatristicDistMatrix <- cophenetic(spttree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
sptdist <- B[!is.na(B$value),]
names(sptdist) <- c("genome1", "genome2", "BGCdist")
sptdist$BGC <- "spt"

slctree <- read.tree("RAxML_bipartitionsBranchLabels.slcBGC")
PatristicDistMatrix <- cophenetic(slctree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
slcdist <- B[!is.na(B$value),]
names(slcdist) <- c("genome1", "genome2", "BGCdist")
slcdist$BGC <- "slc"

slmtree <- read.tree("RAxML_bipartitionsBranchLabels.slmBGC")
PatristicDistMatrix <- cophenetic(slmtree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
slmdist <- B[!is.na(B$value),]
names(slmdist) <- c("genome1", "genome2", "BGCdist")
slmdist$BGC <- "slm"

spotree <- read.tree("RAxML_bipartitionsBranchLabels.spoBGC")
PatristicDistMatrix <- cophenetic(spotree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
spodist <- B[!is.na(B$value),]
names(spodist) <- c("genome1", "genome2", "BGCdist")
spodist$BGC <- "spo"

statree <- read.tree("RAxML_bipartitionsBranchLabels.staBGC")
PatristicDistMatrix <- cophenetic(statree)
B <- melt(PatristicDistMatrix, varnames = c('genome1', 'genome2')) 
stadist <- B[!is.na(B$value),]
names(stadist) <- c("genome1", "genome2", "BGCdist")
stadist$BGC <- "sta"


BGCdist <- rbind(lymdist, lomdist, rifdist, stadist, spodist, slmdist, slcdist, saldist, sptdist)

totaldist <- merge(phylodist, BGCdist, by = c("genome1", "genome2"))

library(ggplot2)

pdf("bgcdistXphylodist.pdf", height = 10, width = 8)
ggplot(totaldist, aes(x = phylodist, y = BGCdist)) +
  stat_smooth(method = "lm", color = "black", size = 2, alpha = 0.7, linetype = "dashed") +
  scale_color_manual(values = c("#8463cc", "#c99044", "#d3d3d3", "#c361aa", "#4aac8b",
                                "#ca5537", "#688dcd", "#89903c", "#ca586f")) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "lightgray", linetype = "dashed", size = 2) +
  geom_point(aes(color = BGC), alpha = 0.5)

dev.off()

summary(lm(totaldist$phylodist ~ totaldist$BGCdist * totaldist$BGC))$adj.r.squared 

