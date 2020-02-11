# network analysis for population structure

setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/antismash/bigscape/') 


library(tidyverse)

###############################################################
# call modularity of the netowkr in GEPHI and get shared BGCs
###############################################################

gephiBGC <- read.table('network-modularity.csv', sep = ',', header = T)

# remove known BGC types
gephiBGCtemp <- gephiBGC[!grepl("BGC", gephiBGC$Id),]

library(splitstackshape)
gephiBGCtemp2 <- cSplit(gephiBGCtemp, "Id", ".")
gephiBGCtemp2$Id_1 <- gsub( "_", "", as.character(gephiBGCtemp2$Id_1))

library(reshape2)

gephiABD <- dcast(gephiBGCtemp2, Id_1 ~ modularity_class)
rownames(gephiABD) <- gephiABD$Id_1
gephiABD$Id_1 <- NULL

# these should match the total BGC counts minus singletons
rowSums(gephiABD)
colSums(gephiABD)

scaled.BGCs <- as.data.frame(apply(gephiABD, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))

library(gplots)

pdf('BGC-abundance.pdf', height = 15, width = 20)

my_palette <- colorRampPalette(c("white", "darkgray", "black"))(n = 200)
heatmap.2(as.matrix(gephiABD),
          main = "PA of BGCs", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace = "row",
          tracecol = "lightgray",
          margins =c(12,9),     # widens margins around plot
          col=my_palette       # use on color palette defined earlier
)  

dev.off()


gephiPA <- gephiABD
gephiPA[gephiPA > 0] <- 1

colSums(gephiPA)

pdf('BGC-PA.pdf', height = 15, width = 20)

my_palette <- colorRampPalette(c("white", "black"))(n = 2)
heatmap.2(as.matrix(gephiPA),
          main = "PA of BGCs", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace = "row",
          tracecol = "lightgray",
          margins =c(12,9),     # widens margins around plot
          col=my_palette       # use on color palette defined earlier
)  

dev.off()


genome.info <- read.table('../bgcs/genomeID.txt', header = T, sep = '\t', row.names = 1)

total.info <- merge(genome.info, gephiABD, by = 0, all = F)
row.names(total.info) <- total.info$Row.names
total.info$clade <- factor(total.info$clade)
total.info$location <- factor(total.info$location)
total.info$clade_color <- factor(total.info$clade_color)

library(vegan)
library(ggplot2)

sol <- metaMDS(gephiABD, distance = "jaccard", k = 2, trymax = 5000)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS = data.frame(x = sol$point[,1], y = sol$point[,2], 
                  clade = as.factor(total.info$clade), 
                  location = as.factor(total.info$location),
                  color = as.factor(total.info$clade_color))

plot.new()
ord <- ordiellipse(sol, total.info$clade, display = "sites", kind ="sd", label = T)
dev.off()


veganCovEllipse <-
  function(cov, center = c(0,0), scale = 1, npoints = 100){
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    o <- attr(Q, "pivot")
    t(center + scale * t(Circle %*% Q[,o]))
  }

df_ell <- data.frame()
for(g in levels(NMDS$clade)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$clade==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center)))
                                ,clade=g))
}

NMDS.mean=aggregate(NMDS[,1:2], list(clade = NMDS$clade), mean)


# vec.sp <- envfit(sol$points, gephiABD, perm = 9999)
# vec.sp.df <- as.data.frame((vec.sp$vectors$arrows * sqrt(vec.sp$vectors$r))*5)
# vec.sp.df$species <- rownames(vec.sp.df)

pdf('mds-noscale-BGC.pdf', height = 15, width = 20)

ggplot(data = NMDS, aes(x, y)) + 
  geom_point(aes(color = color), size = 10, stroke = 3, alpha = 0.8) +
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = clade), size=1, linetype=2, alpha = 0.4) +
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y, label = NMDS.mean$clade, size = 8) + 
  #geom_segment(data = vec.sp.df, aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
  #             arrow = arrow(length = unit(0.5, "cm")), colour = "grey") + 
  #geom_text(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = species), size = 5) +
  coord_fixed() +
  theme_bw() +
  #scale_shape_manual(values = c(4, 3, 2, 5, 0)) +
  scale_colour_identity()

dev.off()

# distance matrix for BGC across
m <- as.matrix(dist(gephiPA))

my_palette <- colorRampPalette(c("black", "darkgray", "gray", "lightgray", "white", "white"))(n = 600)

pdf('shared-BGCs.pdf', width = 10, height = 10)

heatmap.2(m,
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="row")     # only draw a row dendrogram)


dev.off()


dist.mat <- vegdist(gephiPA, method = "jaccard", diag = T)
m = as.matrix(dist.mat)

library(ape)
library(dendextend) 

stree = nj(dist.gene(gephiPA))

myBoots <- boot.phylo(stree, gephiPA, function(xx) nj(dist.gene(xx)), B = 1000,  mc.cores = 6)
stree$node.label <- myBoots

write.tree(stree, file = 'BGC_PA-nj.tre')


# get PA congruent with phylo tree for iToL
genometotal <- read.table("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/total.genome.metadata.txt", sep = '\t', header = T)
rownames(genometotal) <- genometotal$antismashID

PAtotal <- merge(genometotal, gephiPA, by = 0)

write.table(PAtotal, file = "shared-BGCPA.txt", row.names = T, quote = F, sep = '\t')


#############################################################
# can subset any modular of interest for further analysis
#############################################################


# subset the salinosporamide cluster
gephiBGC[gephiBGC$Label == "Lymphostin", ]
gephiBGC[gephiBGC$modularity_class == "17", ]
gephiBGC[gephiBGC$Label == "SalinosporamideA", ]
gephiBGC[gephiBGC$Label == "Lomaiviticin", ]

gephiBGC[gephiBGC$Label == "Rifamycin", ]
gephiBGC[gephiBGC$Label == "Staurosporine", ]
gephiBGC[gephiBGC$Label == "Sioxanthin", ]
gephiBGC[gephiBGC$Label == "Salinichelins", ]
gephiBGC[gephiBGC$Label == "Salinilactam", ]

bigscape <- read.table("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/antismash/bigscape/bigscape_out/network_files/2019-05-07_12-13-24_hybrids_glocal/mix/mix_clustering_c0.30.tsv")
colnames(bigscape) <- c("Id", "GCF")

# sporolide
gephiBGC[gephiBGC$Label == "Sporolide", ]
sporolide <- gephiBGC[gephiBGC$modularity_class == 43, ]
cyanosporaside <- gephiBGC[gephiBGC$modularity_class == 37, ]

bigscape[bigscape$Id == "BGC0000150.1", ]
bigspo <- bigscape[bigscape$GCF == "5004", ]

sposubset <- unique(as.vector(as.matrix(sporolide[,c("Id")])))
sposubset2 <- unique(as.vector(as.matrix(bigspo[,c("Id")])))
sposubset3 <- unique(as.vector(as.matrix(cyanosporaside[,c("Id")])))

spo <- unique(c(sposubset, sposubset2, sposubset3))

# salinilactam
salinilactam <- gephiBGC[gephiBGC$modularity_class == 30, ]
salinilactam2 <- gephiBGC[gephiBGC$modularity_class == 87, ]
bigscape[bigscape$Id == "BGC0000142.1", ]
bigslm <- bigscape[bigscape$GCF == "3926", ]

slmsubset <- unique(as.vector(as.matrix(salinilactam[,c("Id")])))
slmsubset2 <- unique(as.vector(as.matrix(salinilactam2[,c("Id")])))
slmsubset3 <- unique(as.vector(as.matrix(bigslm[,c("Id")])))

slm <- unique(c(slmsubset, slmsubset2, slmsubset3))

#salinipostin
gephiBGC[gephiBGC$Label == "Salinipostin", ]
bigscape[bigscape$Id == "BGC0001458.1", ]

salinipostin <- gephiBGC[gephiBGC$modularity_class == 7, ]
salinipostin2 <- gephiBGC[gephiBGC$modularity_class == 24, ]
bigsap <- bigscape[bigscape$GCF == "1429", ]

sapsubset <- unique(as.vector(as.matrix(salinipostin[,c("Id")])))
sapsubset2 <- unique(as.vector(as.matrix(salinipostin2[,c("Id")])))
sapsubset3 <- unique(as.vector(as.matrix(bigsap[,c("Id")])))

sap <- unique(c(sapsubset, sapsubset2, sapsubset3))


# take all of the nodes in the large cluter in the network
Lymphostin <- gephiBGC[gephiBGC$modularity_class == 11, ]
Lymphostin2 <- gephiBGC[gephiBGC$modularity_class == 19, ]
Lymphostin3 <- gephiBGC[gephiBGC$modularity_class == 17, ]

Lomaiviticin <- gephiBGC[gephiBGC$modularity_class == 35, ]

Staurosporine <- gephiBGC[gephiBGC$modularity_class == 36, ]

Salinilactam <- gephiBGC[gephiBGC$modularity_class == 30, ]



sioxanthin <- gephiBGC[gephiBGC$modularity_class == 47, ]
sioxanthin2 <- gephiBGC[gephiBGC$modularity_class == 25, ]
sioxanthin3 <- gephiBGC[gephiBGC$modularity_class == 61, ]
sioxanthin4 <- gephiBGC[gephiBGC$modularity_class == 3, ]


Desferrioxamine <- gephiBGC[gephiBGC$modularity_class == 45, ]

rifamycin <- gephiBGC[gephiBGC$modularity_class == 8, ]
rifamycin2 <- gephiBGC[gephiBGC$modularity_class == 13, ]

lymsubset <- unique(as.vector(as.matrix(Lymphostin[,c("Id")])))
lymsubset2 <- unique(as.vector(as.matrix(Lymphostin2[,c("Id")])))
lymsubset3 <- unique(as.vector(as.matrix(Lymphostin3[,c("Id")])))

lym <- unique(c(lymsubset, lymsubset2, lymsubset3))

lomsubset <- unique(as.vector(as.matrix(Lomaiviticin[,c("Id")])))
spo <- unique(as.vector(as.matrix(sporolide[,c("Id")])))
sta <- unique(as.vector(as.matrix(Staurosporine[,c("Id")])))
slc <- unique(as.vector(as.matrix(Salinichelins[,c("Id")])))
slm <- unique(as.vector(as.matrix(Salinilactam[,c("Id")])))

siosub <- unique(as.vector(as.matrix(sioxanthin[,c("Id")])))
siosub2 <- unique(as.vector(as.matrix(sioxanthin2[,c("Id")])))
siosub3 <- unique(as.vector(as.matrix(sioxanthin3[,c("Id")])))
siosub4 <- unique(as.vector(as.matrix(sioxanthin4[,c("Id")])))

sio <- unique(c(siosub, siosub2, siosub3, siosub4))



rifsubset <- unique(as.vector(as.matrix(rifamycin[,c("Id")])))
rifsubset2 <- unique(as.vector(as.matrix(rifamycin2[,c("Id")])))
rif <- unique(c(rifsubset, rifsubset2))


write.table(spo, file = "cyanosporaside_search/cya.module.txt", row.names = F, quote = F, col.names = F)


###############################################################
# do this again for what BigScape calls a GCF
###############################################################

gephiBGC <- read.table('bigscape_c0.30_PA.tsv', sep = '\t', header = T)

gephiBGC <- gephiBGC[gephiBGC$ACC != "Salinispora arenicola CNY281",]

rownames(gephiBGC) <- gephiBGC$ACC
gephiBGC$ACC <- NULL

# these should match the total BGC counts minus singletons
rowSums(gephiBGC)
colSums(gephiBGC)



library(gplots)

pdf('bigscape_GCF-PA.pdf', height = 15, width = 20)

my_palette <- colorRampPalette(c("white", "black"))(n = 2)
heatmap.2(as.matrix(gephiBGC),
          main = "PA of BGCs", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace = "row",
          tracecol = "lightgray",
          margins =c(12,9)
)  

dev.off()


genome.info <- read.table('genomeGCF.ids.txt', header = T, sep = '\t', row.names = 1)

total.info <- merge(genome.info, gephiBGC, by = 0, all = F)
row.names(total.info) <- total.info$Row.names
total.info$clade <- factor(total.info$clade)
total.info$location <- factor(total.info$location)
total.info$clade_color <- factor(total.info$clade_color)


library(vegan)
library(ggplot2)

sol <- metaMDS(gephiBGC, distance = "jaccard", k = 2, trymax = 5000)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS = data.frame(x = sol$point[,1], y = sol$point[,2])


NMDS2 <- merge(genome.info, NMDS, by = 0, all = F)
NMDS3 <- NMDS2[,c("Row.names", "x", "y", "clade", "location", "clade_color")]                  

row.names(NMDS3) <- NMDS3$Row.names
NMDS3$Row.names <- NULL

NMDS3$clade <- factor(NMDS3$clade)
NMDS3$location <- factor(NMDS3$location)
NMDS3$clade_color <- factor(NMDS3$clade_color)

plot.new()
ord <- ordiellipse(sol, total.info$clade, display = "sites", kind ="sd", label = T)
dev.off()


veganCovEllipse <-
  function(cov, center = c(0,0), scale = 1, npoints = 100){
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    o <- attr(Q, "pivot")
    t(center + scale * t(Circle %*% Q[,o]))
  }

df_ell <- data.frame()
for(g in levels(NMDS3$clade)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS3[NMDS3$clade==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center)))
                                ,clade=g))
}

NMDS.mean=aggregate(NMDS3[,1:2], list(clade = NMDS3$clade), mean)


pdf('mds-noscale-bigscape_GCF.pdf', height = 15, width = 20)

ggplot(data = NMDS3, aes(x, y)) + 
  geom_point(aes(color = clade_color), size = 10, stroke = 3, alpha = 0.8) +
  #geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = clade), size=1, linetype=2, alpha = 0.4) +
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y, label = NMDS.mean$clade, size = 8) + 
  #geom_segment(data = vec.sp.df, aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
  #             arrow = arrow(length = unit(0.5, "cm")), colour = "grey") + 
  #geom_text(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = species), size = 5) +
  coord_fixed() +
  theme_bw() +
  #scale_shape_manual(values = c(4, 3, 2, 5, 0)) +
  scale_colour_identity()

dev.off()

# distance matrix for BGC across
m <- as.matrix(dist(gephiPA))

library(ecodist)
library(vegan)

dist.mat <- vegdist(gephiBGC, method = "jaccard", diag = T)
clust.res <- hclust(dist.mat)
plot(clust.res)
pops1 <- genome.info[,c(3,4)]
pops2 <- pops1[!rownames(pops1) %in% c("Salinispora arenicola CNY281"), ]
mgroup(dist.mat, pops2, nperm=9999)
#  nclust    mantelr       pval
#      3 0.37574005 0.00010001
#      9 0.34627085 0.00010001
#     11 0.07889497 0.00180018

my_palette <- colorRampPalette(c("black", "darkgray", "gray", "lightgray", "white", "white"))(n = 600)

pdf('shared-BGCs.pdf', width = 10, height = 10)

heatmap.2(m,
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="row")     # only draw a row dendrogram)


dev.off()
