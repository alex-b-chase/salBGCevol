setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/antismash/bigscape')

mixednet <- read.table('bigscape_out/network_files/2019-05-07_12-13-24_hybrids_glocal/mix/mix_c0.30.network', sep = '\t', header = T)


hist(mixednet$Raw.distance, breaks = 600)
hist(mixednet$DSS.index, breaks = 600)
hist(mixednet$Squared.similarity, breaks = 600)



finalletters <- mixednet[, c("Clustername.1", "Clustername.2", "Raw.distance") ]
finalletters$type <- "Directed"
finalletters$BGC1 <-finalletters$Clustername.1
finalletters$BGC2 <- finalletters$Clustername.2

library(splitstackshape)
letterstemp <- cSplit(finalletters, "Clustername.1", ".")
letterstemp2 <- cSplit(letterstemp, "Clustername.2", ".")
letterstemp3 <- letterstemp2[letterstemp2$Clustername.1_1 != "SA_CNY281",]
letterstemp4 <- letterstemp3[letterstemp3$Clustername.2_1 != "SA_CNY281",]

finalletters2 <- letterstemp4[, c(3,4,2,1)]
colnames(finalletters2) <- c("Source", "Target", "Type", "Weight")

write.table(finalletters2, file = "bgc_cluster_network.csv", row.names = F, sep = ',')

########################################################################
####  GET ATTRIBUTES FOR ALL NODES
########################################################################


netatt <- read.table('network.nodes.uniq.txt', sep = '\t', header = F)
colnames(netatt) <- c("cluster1")

library(splitstackshape)
netatt$genome1 <- netatt$cluster1
netatt <- cSplit(netatt, "cluster1", ".")
netatt$cluster1_1 <- gsub( "_", "", as.character(netatt$cluster1_1))
netatt2 <- netatt[netatt$cluster1_1 != "SACNY281",]
colnames(netatt2) <- c("BGC1", "genome1", "cluster1")

colors <- read.table('../bgcs/genomeID.txt', header =T, sep ='\t', stringsAsFactors = T, comment.char = '')
colors2 <- colors[colors$label != "SACNY281", ]
colors2$genome1 <- colors2$label

colors3 <- merge(netatt2, colors2, by = "genome1", all = T)
colors4 <- colors3[,c("BGC1","species", "genome1" ,"cluster1","clade")]
library(dplyr)

colors5 <- colors4 %>% distinct

bgcatt <- read.table('mibig_BGC_innetwork.reduce.txt', header = F, sep = '\t')
colnames(bgcatt) <- c("BGC1", "product")

colors6 <- merge(colors5, bgcatt, by = "BGC1", all = T)

# color known BGC by themselves to distinguish
colors6$clade2 <- ifelse(grepl("BGC", colors6$BGC1, ignore.case = T), "knownBGC", 
                         paste(colors6$clade))
colors6$size <- ifelse(grepl("BGC", colors6$BGC1, ignore.case = T), "5", 
                       "1")
colors6$cluster2 <- ifelse(grepl("BGC", colors6$BGC1, ignore.case = T), paste(colors6$product), 
                           paste(colors6$genome1, colors6$cluster1, sep = '_'))





finalatt <- colors6[, c(1,9,7,8)]
colnames(finalatt) <- c("Id","Label","clade","size")


write.table(finalatt, file = "bgc_cluster_attributes.csv", row.names = F, sep = ',')
