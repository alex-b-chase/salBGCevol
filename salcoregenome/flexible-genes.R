setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/roary/roary85')

genePA <- read.csv('gene_presence_absence.csv', header = T, comment.char = "", na.strings=c("","NA"), stringsAsFactors=FALSE)

# core genes only
corePAF <- genePA[genePA$No..isolates == 118, ]

# flexible genes only
genePAF <- genePA[genePA$No..isolates != 118, ]

flexgen <- genePAF[c(1,15:132)]
rownames(flexgen) <- flexgen$Gene
flexgen$Gene <- NULL
flexgen[flexgen==""]  <- NA 
flexgen[is.na(flexgen)] <- 0

flexgen[] <- lapply(flexgen, function(x) as.integer(x!="0"))

library(vegan)
tflex <- t(flexgen)

dist.mat <- vegdist(tflex, method = "jaccard", diag = T)
clust.res <- hclust(dist.mat)
plot(clust.res)

########### ########### ########### ########### ########### ########### ########### ########### 
# stats to see if strains within populations and/or within site share more similar genes
########### ########### ########### ########### ########### ########### ########### ########### 

pops <- read.table('../../total.genome.metadata.txt', header = T, row.names = 1, sep = '\t')
pops <- pops[pops$genomefileID != "SACNY281",]

total.info <- merge(pops, tflex, by=0, all=F)
row.names(total.info) <- total.info$Row.names
total.info$Row.names <- NULL

total.info$clade <- factor(total.info$clade)
total.info$species <- factor(total.info$species)
total.info$location <- factor(total.info$location)

library(ecodist)
pops1 <- pops[,c(3,5,6)]
mgroup(dist.mat, pops1, nperm=9999)
#             nclust    mantelr  pval
# species      3 0.89132676 0.00010001 ***
# clade        9 0.90992590 0.00010001 ***
# location    11 0.06031434 0.01080108 *


vare.mds <- rda(tflex, method = "jaccard")
plot(vare.mds, display = "sites")

cmd <- cmdscale(dist.mat)
ordiplot(cmd, display = "sites", type = "text")

sol <- metaMDS(tflex, distance = "jaccard", k = 2, trymax = 500)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS = data.frame(x = sol$point[,1], y = sol$point[,2], 
                  location = as.factor(total.info$location), 
                  species = as.factor(total.info$species),
                  clade = as.factor(total.info$clade))

attach(NMDS)
sp.ano <- anosim(dist.mat, species, distance = "jaccard", permutations = 999)
summary(sp.ano)
# ANOSIM statistic R: 0.9913 
# Significance: 0.001 
# Permutation: free
# Number of permutations: 999

clade.ano <- anosim(dist.mat, clade, distance = "jaccard", permutations = 999)
summary(clade.ano)
# ANOSIM statistic R: 0.9999
# Significance: 0.001 
# Permutation: free
# Number of permutations: 9999

site.ano <- anosim(dist.mat, location, distance = "jaccard", permutations = 9999)
summary(site.ano)
# ANOSIM statistic R: 0.05592 
# Significance: 0.0454 
# Permutation: free
# Number of permutations: 9999

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

NMDS.mean=aggregate(NMDS[,1:2],list(clade=NMDS$clade),mean)

vec.sp <- envfit(sol$points, tflex, perm = 999)
vec.sp.df <- as.data.frame((vec.sp$vectors$arrows * sqrt(vec.sp$vectors$r))*5)
vec.sp.df$species <- rownames(vec.sp.df)

library(ggplot2)

pdf("../pop-specific-genes/flexible-genome-MDS.pdf", width = 12, height = 12)

ggplot(data = NMDS, aes(x, y)) + 
  #geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = subcluster), size=1, linetype=2, alpha = 0.4) +
  geom_point(aes(color = clade, shape = location), size = 10, stroke = 3) +
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y, label = NMDS.mean$subcluster, size = 8) +
  #coord_fixed() +
  theme_bw() +
  scale_colour_manual(values = c("black", "orange1", "dodgerblue", "red",
                                 "firebrick", "darkgoldenrod4", "purple1", "plum3",
                                 "orchid3", "forestgreen", "green2", "lightgreen",
                                 "darkseagreen4", "lightcoral", "lightpink")) 
  #theme(legend.position="bottom")

dev.off()


########### ########### ########### ########### ########### ########### ########### ########### 
# plot jaccard matrix and NJ tree
########### ########### ########### ########### ########### ########### ########### ########### 


dist.mat <- vegdist(tflex, method = "jaccard", diag = T)
m = as.matrix(dist.mat)

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("black", "darkgray","gray", "lightgray", "white"))(n = 299)

library(ape)
library(dendextend) 
library(phytools)



stree = nj(dist.gene(tflex))
# x <- root(stree, outgroup = "Curtobacterium_MCBA15007.bgpipe")

myBoots <- boot.phylo(stree, tflex, function(xx) nj(dist.gene(xx)), B = 100,  mc.cores = 2)
stree$node.label <- myBoots
x <- midpoint.root(stree)
plot(x)

write.tree(stree, file = 'flexgenome-nj-popspec.tree')
write.tree(x, file = 'flexgenome-nj-rooted.tree')

tip_order <- c("StCNB536_BA", "StCNB476_BA", "StCNS416_BA", "StCNY012_BA", "StCNH898_BA", "StCNT261_BA", "StCNB440_BA", 
               "StCNY678_YU", "StCNY681_YU", "StCNT250_BA", "StCNR699_BA", "StCNS197_BA", "SpCNY646_RS", "SpCNS237_PL", 
               "SpDSM45549_FJ", "SpCNY202_SC", "SpCNR942_PL", "SpCNT569_FJ", "SpCNT854_HI", "SpCNT124_FJ", "SpCNT584_FJ", 
               "SpCNT029_FJ", "SpDSM45547_FJ", "SpCNY673_MD", "SpCNY703_MD", "SpCNT403_FJ", "SpCNS996_FJ", "SpCNT045_FJ", 
               "SpCNS860_FJ", "SpDSM45543_FJ", "SpCNY666_MD", "SpCNS055_PL", "SpCNS801_FJ", "SpDSM45548_FJ", "SpCNT609_FJ", 
               "SpCNT084_FJ", "SpCNT133_FJ", "SpCNT796_HI", "SpCNT851_HI", "SpCNY330_SC", "SpCNY331_SC", "SpCNY363_SC", 
               "SpCNT855_HI", "SpDSM45544_FJ", "SpCNT003_FJ", "SpCNS103_PL", "SpCNR909_PL", "SpCNQ768_GU", "SpCNR114_GU", 
               "SpCNR510_GU", "SpCNY498_PV", "SpCNH732_RS", "SpCNY239_FJ", "SpCNR894_PL", "SpCNT131_FJ", "SpCNT001_FJ", 
               "SpCNT603_FJ", "SaCNB458_BA", "SaCNY694_YU", "SaCNB527_BA", "SaCNY685_YU", "SaCNY690_YU", "SaCNH962_SC", 
               "SaCNH963_SC", "SaCNY486_PV", "SaCNH996B_SC", "SaCNH996_SC", "SaCNH713_RS", "SaCNH905_BA", "SaCNH643_BA", 
               "SaCNH646_BA", "SaCNY679_YU", "SaCNH877_BA", "SaCNY011_BA", "SaCNH964_SC", "SaCNP105_SC", "SaCNH941_SC", 
               "SaCNP193_SC", "SaDSM45545_FJ", "SaCNS848_FJ", "SaCNY280_FJ", "SaCNT799_HI", "SaCNT850_HI", "SaCNT857_HI", 
               "SaCNT859_HI", "SaCNT849_HI", "SaCNT798_HI", "SaCNT800_HI", "SaCNX508_PM", "SaCNX891_PM", "SaCNX814_PM", 
               "SaCNX481_PM", "SaCNX482_PM", "SaCNH718_RS", "SaCNR921_PL", "SaCNR107_GU", "SaCNS051_PL", "SaCNR425_GU", 
               "SaCNQ748_GU", "SaCNS243_PL", "SaCNS205_PL", "SaCNS325_PL", "SaCNS296_PL", "SaCNS299_PL", "SaCNT005_FJ", 
               "SaCNY231_FJ", "SaCNS744_FJ", "SaCNY260_FJ", "SaCNQ884_GU", "SaCNS342_FJ", "SaCNY237_FJ", "SaCNY244_FJ", 
               "SaCNS820_FJ", "SaCNY230_FJ", "SaCNY234_FJ", "SaCNY282_FJ", "SaCNS673_FJ", "SaCNY256_FJ")

ord.mat <- m[tip_order, tip_order]
write.table(ord.mat, file = "flexgenome-jaccard.txt", row.names = T, col.names = T, quote = F, sep = '\t')

coretip_order <- c("StCNY012_BA",
                   "StCNH898_BA",
                   "StCNT261_BA",
                   "StCNB536_BA",
                   "StCNB476_BA",
                   "StCNS416_BA",
                   "StCNT250_BA",
                   "StCNR699_BA",
                   "StCNS197_BA",
                   "StCNB440_BA",
                   "StCNY678_YU",
                   "StCNY681_YU",
                   "SpCNT569_FJ",
                   "SpCNR942_PL",
                   "SpCNY202_SC",
                   "SpCNY646_RS",
                   "SpDSM45549_FJ",
                   "SpCNS237_PL",
                   "SpCNT854_HI",
                   "SpCNT124_FJ",
                   "SpCNT584_FJ",
                   "SpCNT029_FJ",
                   "SpDSM45547_FJ",
                   "SpCNY703_MD",
                   "SpCNY673_MD",
                   "SpCNS996_FJ",
                   "SpCNT045_FJ",
                   "SpCNT403_FJ",
                   "SpCNS860_FJ",
                   "SpDSM45543_FJ",
                   "SpCNY666_MD",
                   "SpCNS055_PL",
                   "SpCNS801_FJ",
                   "SpDSM45548_FJ",
                   "SpCNT609_FJ",
                   "SpCNT084_FJ",
                   "SpCNT133_FJ",
                   "SpCNT796_HI",
                   "SpCNT851_HI",
                   "SpCNY330_SC",
                   "SpCNY363_SC",
                   "SpCNY331_SC",
                   "SpCNT855_HI",
                   "SpCNY498_PV",
                   "SpCNT003_FJ",
                   "SpDSM45544_FJ",
                   "SpCNR909_PL",
                   "SpCNS103_PL",
                   "SpCNQ768_GU",
                   "SpCNR114_GU",
                   "SpCNH732_RS",
                   "SpCNR894_PL",
                   "SpCNR510_GU",
                   "SpCNY239_FJ",
                   "SpCNT001_FJ",
                   "SpCNT603_FJ",
                   "SpCNT131_FJ",
                   "SaCNY685_YU",
                   "SaCNB527_BA",
                   "SaCNY690_YU",
                   "SaCNB458_BA",
                   "SaCNY694_YU",
                   "SaCNH963_SC",
                   "SaCNH962_SC",
                   "SaCNY486_PV",
                   "SaCNH996_SC",
                   "SaCNH996B_SC",
                   "SaCNH713_RS",
                   "SaCNP105_SC",
                   "SaCNH964_SC",
                   "SaCNP193_SC",
                   "SaCNH941_SC",
                   "SaCNH877_BA",
                   "SaCNY011_BA",
                   "SaCNH643_BA",
                   "SaCNH646_BA",
                   "SaCNH905_BA",
                   "SaCNY679_YU",
                   "SaDSM45545_FJ",
                   "SaCNY280_FJ",
                   "SaCNS848_FJ",
                   "SaCNT859_HI",
                   "SaCNT857_HI",
                   "SaCNT799_HI",
                   "SaCNT850_HI",
                   "SaCNT849_HI",
                   "SaCNT800_HI",
                   "SaCNT798_HI",
                   "SaCNX482_PM",
                   "SaCNX814_PM",
                   "SaCNX508_PM",
                   "SaCNX891_PM",
                   "SaCNX481_PM",
                   "SaCNH718_RS",
                   "SaCNS299_PL",
                   "SaCNS243_PL",
                   "SaCNR921_PL",
                   "SaCNQ884_GU",
                   "SaCNR107_GU",
                   "SaCNS051_PL",
                   "SaCNS296_PL",
                   "SaCNS205_PL",
                   "SaCNS325_PL",
                   "SaCNQ748_GU",
                   "SaCNR425_GU",
                   "SaCNT005_FJ",
                   "SaCNY231_FJ",
                   "SaCNS744_FJ",
                   "SaCNS820_FJ",
                   "SaCNY260_FJ",
                   "SaCNS673_FJ",
                   "SaCNY256_FJ",
                   "SaCNY230_FJ",
                   "SaCNY237_FJ",
                   "SaCNY282_FJ",
                   "SaCNY234_FJ",
                   "SaCNS342_FJ",
                   "SaCNY244_FJ"
)

ord.mat <- m[coretip_order, coretip_order]
write.table(ord.mat, file = "flexgenome-jaccard-corephylo.txt", row.names = T, col.names = T, quote = F, sep = '\t')


pdf('flexible-genome.pdf', width = 10, height = 10)
heatmap.2(ord.mat,
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="row")     # only draw a row dendrogram)


dev.off()

########### get population specific genes

# set populations
pop1 <- c("StCNB440_BA", "StCNB476_BA", "StCNB536_BA", "StCNH898_BA", "StCNR699_BA", "StCNS197_BA", "StCNS416_BA", "StCNT250_BA", "StCNT261_BA", "StCNY012_BA", "StCNY678_YU", "StCNY681_YU")
pop2 <- c("SpCNR942_PL", "SpCNT569_FJ")
pop4 <- c("SpCNS237_PL", "SpCNY646_RS", "SpDSM45549_FJ")
pop5 <- c("SpCNS860_FJ", "SpCNS996_FJ", "SpCNT029_FJ", "SpCNT045_FJ", "SpCNT124_FJ", "SpCNT403_FJ", "SpCNT584_FJ", "SpCNT854_HI", "SpCNY673_MD", "SpCNY703_MD", "SpDSM45543_FJ", "SpDSM45547_FJ") 
pop7 <- c("SpCNS055_PL", "SpCNS801_FJ", "SpDSM45548_FJ")
pop8 <- c("SpCNH732_RS", "SpCNQ768_GU", "SpCNR114_GU", "SpCNR510_GU", "SpCNR894_PL", "SpCNR909_PL", "SpCNS103_PL", "SpCNT001_FJ", "SpCNT003_FJ", "SpCNT084_FJ", "SpCNT131_FJ", "SpCNT133_FJ", "SpCNT603_FJ", "SpCNT609_FJ", "SpCNT796_HI", "SpCNT851_HI", "SpCNT855_HI", "SpCNY239_FJ", "SpCNY330_SC", "SpCNY331_SC", "SpCNY363_SC", "SpCNY498_PV", "SpDSM45544_FJ")
pop10 <- c( "SaCNB458_BA", "SaCNB527_BA", "SaCNH643_BA", "SaCNH646_BA", "SaCNH713_RS", "SaCNH718_RS", "SaCNH877_BA", "SaCNH905_BA", "SaCNH941_SC", "SaCNH962_SC", "SaCNH963_SC", "SaCNH964_SC", "SaCNH996_SC", "SaCNH996B_SC", "SaCNP105_SC", "SaCNP193_SC", "SaCNQ748_GU", "SaCNQ884_GU", "SaCNR107_GU", "SaCNR425_GU", "SaCNR921_PL", "SaCNS051_PL", "SaCNS205_PL", "SaCNS243_PL", "SaCNS296_PL", "SaCNS299_PL", "SaCNS325_PL", "SaCNS342_FJ", "SaCNS673_FJ", "SaCNS744_FJ", "SaCNS820_FJ", "SaCNS848_FJ", "SaCNT005_FJ", "SaCNT798_HI", "SaCNT799_HI", "SaCNT800_HI", "SaCNT849_HI", "SaCNT850_HI", "SaCNT857_HI", "SaCNT859_HI", "SaCNX481_PM", "SaCNX482_PM", "SaCNX508_PM", "SaCNX814_PM", "SaCNX891_PM", "SaCNY011_BA", "SaCNY230_FJ", "SaCNY231_FJ", "SaCNY234_FJ", "SaCNY237_FJ", "SaCNY244_FJ", "SaCNY256_FJ", "SaCNY260_FJ", "SaCNY280_FJ", "SaCNY282_FJ", "SaCNY486_PV", "SaCNY679_YU", "SaCNY685_YU", "SaCNY690_YU", "SaCNY694_YU", "SaDSM45545_FJ")

tokeep <- c("Gene", "Annotation", "No..isolates", pop8)

popPA <- genePA[genePA$No..isolates <= length(pop8), ]

popPA_sub <- popPA[ , (names(popPA) %in% tokeep)]
popPA_final <- popPA_sub[complete.cases(popPA_sub), ]


write.table(popPA_final, file = "../species-specific-genes/pacifica-specific-genes.txt", row.names = F, quote = F, sep = '\t')

# run subset-popgenes.sh after you generate each population to run

##########################################
# plot gene positions
##########################################

species <- c("arenicola")
genes <- read.table(paste('../species-specific-genes/', species, '.genepos.txt', sep = ''), header = F, sep = '\t')
names(genes) <- c("contig", "Gene", "genome", "start", "end", "direction")

genes$genomecontig <- paste(genes$genome, "_", genes$contig, sep = '')

contig_sizes <- read.table('../species-specific-genes/contiglength_final.txt', header = F, sep = '\t')
names(contig_sizes) <- c("genome","contig", "size")
contig_sizes$genomecontig <- paste(contig_sizes$genome, "_", contig_sizes$contig, sep = '')

totalgenes <- merge(genes, contig_sizes, by = c("genomecontig"))

totalgenes <- totalgenes[!duplicated(totalgenes), ]


# only use genomic regions that have >=3 continuous genes
totalgenes2 <- totalgenes[totalgenes$genomecontig %in% names(which(table(totalgenes$genomecontig) > 3)), ]
totalgenes3 <- totalgenes2[totalgenes2$Gene %in% names(which(table(totalgenes2$Gene) >= 3)), ]
# make sure all genomes in species have the genes
if (length(unique(totalgenes3$genome.x)) == length(unique(totalgenes$genome.x))){
  totalgenes4 <- totalgenes3
} else {
  totalgenes4 <- NULL
}

goodcontigs <- unique(totalgenes4$genomecontig)

good_sizes <- contig_sizes[ (contig_sizes$genomecontig %in% goodcontigs), ]
table(good_sizes$genomecontig)
length(unique(totalgenes4$Gene))

# create an ordered factor level to use for the chromosomes in all the datasets
totalgenes4$genomecontig <- factor(totalgenes4$genomecontig)
totalgenes4$genome.x <- factor(totalgenes4$genome.x)

contig_order <- unique(totalgenes4$genomecontig)
contig_key <- setNames(object = as.character(seq(1, length(unique(totalgenes4$genomecontig)))), 
                       nm = contig_order)
contig_order <- factor(x = contig_order, levels = rev(contig_order))

# convert the chromosome column in each dataset to the ordered factor
good_sizes[["genomecontig"]] <- factor(x = good_sizes[["genomecontig"]], 
                                 levels = contig_order)
totalgenes4[["genomecontig"]] <- factor(x = totalgenes4[["genomecontig"]], 
                                  levels = contig_order)

totalgenes4$length <- abs(totalgenes4$start - totalgenes4$end)
totalgenes5 <- totalgenes4[totalgenes4$genome.x == "SaCNS205_PL",]
good_sizes2 <- good_sizes[good_sizes$genome == "SaCNS205_PL",]

library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

pdf(paste('../species-specific-genes/', species, '-gene-annotations.pdf', sep = ''), width = 24, height = 12)

ggplot(data = good_sizes2) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(genomecontig) - 0.2, 
                xmax = as.numeric(genomecontig) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") +
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "contig", limits = names(contig_key)) +
  # add bands for gene value
  geom_rect(data = totalgenes5, aes(xmin = as.numeric(genomecontig) - 0.2, 
                                   xmax = as.numeric(genomecontig) + 0.2, 
                                   ymax = end, ymin = start)) +
  geom_text_repel(data = totalgenes5, 
                  aes(x = genomecontig, y = start, label = Gene), 
                  color = "red", show.legend = FALSE)

dev.off()

annotation <- read.table(paste("../species-specific-genes/", species, "-specific-genes.txt", sep = ''), header = T, sep = "\t", quote = "")
annotation2 <- annotation[, c(1,2) ]

totalanno4 <- merge(totalgenes4, annotation2, by = "Gene") 
write.table(totalanno4, file = paste('../species-specific-genes/', species, '-localized-genes.txt', sep = ''), sep = '\t', row.names = F, quote = F)

setwd("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/roary/species-specific-genes/")
arenicola <- read.table("arenicola-localized-genes.txt", header = T, sep = '\t', quote = "")

SACNS205 <- arenicola[arenicola$genomecontig == "SaCNS205_PL_contig_0", ]

contig_order <- unique(SACNS205$genomecontig)
contig_key <- setNames(object = as.character(seq(1, length(unique(SACNS205$genomecontig)))), 
                       nm = contig_order)

library(ggplot2)
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

pdf("arenicola_SACNS205only.pdf", width = 12, height = 4)
ggplot(data = SACNS205) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(genomecontig) - 0.2, 
                xmax = as.numeric(genomecontig) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") +
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "red"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "contig", limits = names(contig_key)) +
  # add bands for gene value
  geom_rect(data = SACNS205, aes(xmin = as.numeric(genomecontig) - 0.2, 
                                 xmax = as.numeric(genomecontig) + 0.2, 
                                 ymax = end, ymin = start))

dev.off()



tropica <- read.table("tropica-localized-genes.txt", header = T, sep = '\t', quote = "")

CNB440 <- tropica[tropica$genomecontig == "StCNB440_BA_contig_0", ]

contig_order <- unique(CNB440$genomecontig)
contig_key <- setNames(object = as.character(seq(1, length(unique(CNB440$genomecontig)))), 
                       nm = contig_order)

library(ggplot2)
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

pdf("tropica_STCNB440only.pdf", width = 12, height = 4)
ggplot(data = CNB440) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(genomecontig) - 0.2, 
                xmax = as.numeric(genomecontig) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") +
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "red"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "contig", limits = names(contig_key)) +
  # add bands for gene value
  geom_rect(data = CNB440, aes(xmin = as.numeric(genomecontig) - 0.2, 
                                 xmax = as.numeric(genomecontig) + 0.2, 
                                 ymax = end, ymin = start))

dev.off()