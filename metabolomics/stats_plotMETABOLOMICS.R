rm(list=ls())
setwd("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/metabolomics-data/")

library(vegan)
library(ggplot2)
library(ecodist)
library(RFmarkerDetector) # https://rdrr.io/cran/RFmarkerDetector/
# RFmarkerDetector: Multivariate Analysis of Metabolomics Data using Random Forests

# read in metadata
metadata <- read.table("msIDs.txt", header = T, sep = '\t', row.names = 1, comment.char = "")
staxcolors <- setNames(as.character(metadata$clade_color), metadata$clade_name)

totalMS <- read.table(file = "allMS1data.txt", header = T, sep = '\t', row.names = 1)

totalMS[totalMS <= 1E3] <- NA
# remove all rows with NA - meaning only found in media
totalMS2 <- totalMS[rowSums(is.na(totalMS)) != ncol(totalMS), ]

my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
totalMS2[is.na(totalMS2)] <- my.min(totalMS2) /10

tflex <- t(totalMS2)
logtflex <- log(tflex, 2)
paretotflex <- paretoscale(logtflex)
scaletflex <- scale(logtflex)

dist.mat <- vegdist(paretotflex, method = "euclidean", diag = T)
clust.res <- hclust(dist.mat)
plot(clust.res)

hist(t(scaletflex))
hist(t(paretotflex))

################################################
################################################
#### Principal Component Analysis (PCoA)
################################################
################################################

# Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables).
pcatflex <- merge(metadata, paretotflex, by = 0)
rownames(pcatflex) <- pcatflex$Row.names

pcatflex$clade_name <- droplevels(pcatflex$clade_name)

# PERMANOVA stats analysis
perm.pareto <- adonis2(dist.mat ~ clade_name, data=pcatflex, permutations = 999, method = "euclidean", strata = "PLOT")
perm.pareto


# PCoA analysis for MS data
pca <- prcomp(t(paretotflex), center = F, scale = F)
pcaresults <- summary(pca)
pcaresults$importance[3,1:3] 

pcadata <- as.data.frame(pcaresults$rotation)
pcadataPC12 <- pcadata[, c(1:3)]
 
pcameta <- merge(metadata, pcadataPC12, by = 0, all = F)
staxcolors <- setNames(as.character(pcameta$clade_color), pcameta$clade_name)

pdf("allMS1data.pdf", height = 10, width = 12)
ggplot(pcameta, aes(PC1, PC2)) +
  geom_point(aes(color = clade_name), size = 4) +
  geom_text(aes(label = strainID)) +
  theme_bw() +
  scale_colour_manual(values = staxcolors)

dev.off()


################################################
################################################
###   now removed wash contaminansts
###   and substracted media components from MS
################################################
################################################

# have pseudo replication in the samples, average by strain after ordination with error bars
strainMS <- read.table(file = "nomediaMS1data.txt", header = T, sep = '\t', row.names = 1)

strainMS[strainMS <= 1E3] <- NA
# remove all rows with NA - meaning only found in media
strainMS2 <- strainMS[rowSums(is.na(strainMS)) != ncol(strainMS), ]
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
strainMS2[is.na(strainMS2)] <- my.min(strainMS2) /10

# now can normalize and scale data

tflex <- t(strainMS2)
logstrainavg <- log(tflex, 2)
plogstrainavg <- paretoscale(logstrainavg)

plogstraindist <- vegdist(plogstrainavg, method = "euclidean", diag = T)
plogstrainclust <- hclust(plogstraindist)
plot(plogstrainclust)

hist(t(plogstrainavg))

# PERMANOVA stats analysis
# need to average features by strain before we can run
# read in metadata to plot the PCoA graph
metadata <- read.table("msIDs.txt", header = T, sep = '\t', row.names = 1, comment.char = "")
metadata$Var1 <- rownames(metadata)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

library(reshape2)
strainmerge <- melt(plogstrainavg)
strainMETA <- merge(strainmerge, metadata, by = "Var1")
summPERM <- summarySE(strainMETA, measurevar = "value", groupvars = c("strainID", "Var2"))

strainPERM <- dcast(summPERM, strainID ~ Var2, value.var = "value")
rownames(strainPERM) <- strainPERM$strainID
strainPERM$strainID <- NULL

strainmeta <- read.table("../genomeID.txt", sep = '\t', comment.char = "", header = T, row.names = 1)

strainpcatflex <- merge(strainmeta, strainPERM, by = 0)
rownames(strainpcatflex) <- strainpcatflex$Row.names
strain.perm.pareto <- adonis2(strainPERM ~ clade_name, data=strainpcatflex, permutations = 999, method = "euclidean", strata = "PLOT")
strain.perm.pareto


### now plot the PCoA plot
strainpca <- prcomp(t(plogstrainavg), center = F, scale = F)
strainpcaresults <- summary(strainpca)
strainpcaresults$importance[3,1:3] 

strainpcadata <- as.data.frame(strainpcaresults$rotation)
strainpcadataPC12 <- strainpcadata[, c(1:3)]

strainpcameta <- merge(metadata, strainpcadataPC12, by = 0, all = F)
taxcolors <- setNames(as.character(strainpcameta$clade_color), strainpcameta$clade_name)

summPC1 <- summarySE(strainpcameta, measurevar = "PC1", groupvars = c("strainID", "clade_name", "clade_color"))
summPC1$PC1sd <- summPC1$sd
summPC2 <- summarySE(strainpcameta, measurevar = "PC2", groupvars = c("strainID", "clade_name", "clade_color"))
summPC2$PC2sd <- summPC2$sd

totalPC12temp <- merge(summPC1, summPC2, by = c("strainID", "clade_name", "clade_color"))
totalPC12 <- totalPC12temp[, c("clade_name", "PC1", "PC2", "PC1sd", "PC2sd", "clade_color", "strainID")]
staxcolors <- setNames(as.character(totalPC12$clade_color), totalPC12$clade_name)

library(ellipse)
centroids <- aggregate(cbind(PC1, PC2) ~ clade_name, totalPC12, mean)

df_ell <- data.frame() 
for(g in levels(totalPC12$clade_name)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(totalPC12[totalPC12$clade_name==g,], 
      ellipse(cor(PC1, PC2),scale=c(sd(PC1),sd(PC2)),centre=c(mean(PC1),mean(PC2))))),clade_name=g))} 

pdf("nomediaMS1data-avgstrain.pdf", height = 10, width = 12)
ggplot(totalPC12, aes(PC1, PC2)) +
  geom_path(data = df_ell, aes(x = x, y = y, color = clade_name)) +
  geom_errorbarh(aes(xmin = PC1 - PC1sd, xmax = PC1 + PC1sd), height = 0) +
  geom_errorbar(aes(ymin = PC2 - PC2sd, ymax = PC2 + PC2sd), width = 0) +
  geom_point(aes(color = clade_name), size = 6, shape = 15) +
  #geom_text(aes(label = strainID)) +
  theme_bw() +
  scale_colour_manual(values = staxcolors) +
  scale_fill_manual(values = staxcolors)

dev.off()


################################################
################################################
###########    Targeted metabolites  ########### 
################################################
################################################

targetfeatures <- read.table("targetedfeatures.txt", header = F, sep = '\t')
# remove redundant rifS features - also rifS had column carry over
targetfeatures <- targetfeatures[targetfeatures$V2 != "rifS", ]
rownames(targetfeatures) <- targetfeatures$V1
targetfeatures$V1 <- NULL

targetsubset <- merge(targetfeatures, strainMS2, by = 0)
rownames(targetsubset) <- paste(targetsubset$Row.names, targetsubset$V2, sep = "_")
targetsubset$Row.names <- NULL
targetsubset$V2 <- NULL

targflex <- t(targetsubset)
logtarget <- log10(targflex)

metadata <- read.table("msIDs.txt", header = T, sep = '\t', row.names = 1, comment.char = "")
metadata$Var1 <- rownames(metadata)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

library(reshape2)
targmerge <- melt(logtarget)
targMETA <- merge(targmerge, metadata, by = "Var1")
summTARG <- summarySE(targMETA, measurevar = "value", groupvars = c("strainID", "Var2"))

strainTARG <- dcast(summTARG, strainID ~ Var2, value.var = "value")
rownames(strainTARG) <- strainTARG$strainID
strainTARG$strainID <- NULL

# merge the results from the MS2 network, if strain had the molecule
# if FBMN over networked, we will remove the data point
# multiply the two matrices
targMS2 <- read.table("targetedfeatures.MS2.txt", header = T, sep = '\t')
targMSd <- dcast(targMS2, strain ~ featureID, value.var = "value")
targMSd[is.na(targMSd)] <- 0
rownames(targMSd) <- targMSd$strain
targMSd$strain <- NULL
targMS2temp <- melt(as.matrix(targMSd))
colnames(targMS2temp) <- c("strain", "featureID", "presence")

strainTARGtemp <- melt(as.matrix(strainTARG))
colnames(strainTARGtemp) <- c("strain", "featureID", "logint")
targMSTOT <- merge(strainTARGtemp, targMS2temp, by = c("strain", "featureID"))
targMSTOT$int.value <- targMSTOT$logint * targMSTOT$presence

targMSTOTDF <- dcast(targMSTOT, strain ~ featureID, value.var = "int.value")
targMSTOTDF[targMSTOTDF == 0] <- NA
targMSTOTDF[is.na(targMSTOTDF)] <- 1
rownames(targMSTOTDF) <- targMSTOTDF$strain
targMSTOTDF$strain <- NULL

# now read in the BGC presence/absence
repstrains <- read.table("repstrains-MS.txt", header = T, sep = "\t", row.names = 1)
repstrains[repstrains == 1] <- 10
bgcMS <- merge(repstrains, targMSTOTDF, by = 0)
rownames(bgcMS) <- bgcMS$Row.names
bgcMS$Row.names <- NULL

library(gplots)

my_palette <- colorRampPalette(c("white", "gray94", "gray87", 
                                 "gray77", "gray60", "gray25",  
                                 "black", "black", "black", "green"))(n = 100)
colorder <- c("spt", "feat1727_sptK",
              "lym", "feat2330_lym", "feat3548_lymol", 
              "feat280_neolymA", "feat2498_neolymolA",
              "feat317_neolymB", "feat288_neolymolB", "feat651_neolymC",
              "lom", "feat122_lomC", "feat2375_lomD", "feat2837_lomE",
              "sta", "feat169_sta", "feat1615_sta7", "feat2460_staoxo",
              "sal", "feat444_salA", "feat192_salB", "feat1777_salD", "feat475_salE",
              "feat1812_salG", "feat1203_salK",
              "slm", "feat1513_slm", 
              "spo", "feat1798_spoA", "feat1820_spoB", "slc",
              "rif", "feat1895_rifS2", "feat2446_rif8deso")

hclustfunc <- function(x) hclust(x, method = "complete")
distfunc <- function(x) dist(x, method = "euclidean")

pdf("targetmetabolomics.pdf", height = 12, width = 14)
heatmap.2(as.matrix(bgcMS[, colorder]), Colv=FALSE, 
          dendrogram="row", trace="none", 
          margin = c(8,9), col=my_palette,
          hclust = hclustfunc, distfun = distfunc)
dev.off()

# permanova with bgcs included
bgcsMStotal <- merge(strainmeta, bgcMS, by = 0)
rownames(bgcsMStotal) <- bgcsMStotal$Row.names
strain.perm.pareto <- adonis2(bgcMS ~ clade_name, data=bgcsMStotal, permutations = 999, method = "euclidean", strata = "PLOT")
strain.perm.pareto

# permanova with bgcs excluded
targMStotal <- merge(strainmeta, targMSTOTDF, by = 0)
rownames(targMStotal) <- targMStotal$Row.names
strain.perm.pareto <- adonis2(targMSTOTDF ~ clade_name, data=targMStotal, permutations = 999, method = "euclidean", strata = "PLOT")
strain.perm.pareto

