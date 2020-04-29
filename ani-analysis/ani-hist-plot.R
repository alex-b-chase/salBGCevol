rm(list=ls())
setwd("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCanalysis/evolprocesses/compANI")

genomeani <- read.table('total.genome.txt', sep = '\t', header = F)

genomeidtemp <- read.table("../../../genomeID.txt", sep = '\t', header = T, comment.char = "")
genomeid <- genomeidtemp[, c(1,5,7)]

colnames(genomeid) <- c("strain1", "clade1", "clade1color")
genomeid2 <- genomeid
colnames(genomeid2) <- c("strain2", "clade2", "clade2color")

colnames(genomeani) <- c("strain1", "strain2", "ani")

ani1 <- merge(genomeani, genomeid, by = "strain1")
ani2 <- merge(ani1, genomeid2, by = "strain2")

ani2$withinclade <- NA
ani2[ani2$clade1 != ani2$clade2, "withinclade"] <- "no"
ani2[ani2$clade1 == ani2$clade2, "withinclade"] <- "yes"

library(plyr)
mu <- ddply(ani2, "withinclade", summarise, grp.mean=mean(ani))

library(ggplot2)
WGani <- ggplot(ani2, aes(x = ani, color = withinclade, fill = withinclade)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 200) +
  theme_classic() +
  geom_density(alpha = 0.6) +
  scale_color_manual(values = c("lightgray", "black")) +
  scale_fill_manual(values = c("lightgray", "black")) +
  geom_vline(xintercept = 95, color = "red", linetype = "dashed") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(80,100))


# run all BGCs at once now
allBGCANI <- read.table("total.bgcANI.txt", sep = ';', header = F)

allBGCANI$V4 <- allBGCANI$V1
allBGCANI$V5 <- allBGCANI$V2
allBGCANI$V1 <- gsub("\\..*", "", allBGCANI$V1)
allBGCANI$V2 <- gsub("\\..*", "", allBGCANI$V2)

colnames(allBGCANI) <- c("strain1", "strain2", "ani", "cluster1", "cluster2")

hist(allBGCANI$ani, breaks = 200)

bgcani1 <- merge(allBGCANI, genomeid, by = "strain1")
bgcani2 <- merge(bgcani1, genomeid2, by = "strain2")

bgcani2$withinspecies <- NA
bgcani2[bgcani2$clade1 != bgcani2$clade2 & bgcani2$ani < 95, "withinspecies"] <- "no"
bgcani2[bgcani2$clade1 == bgcani2$clade2 & bgcani2$ani > 95, "withinspecies"] <- "yes"
bgcani2[bgcani2$clade1 != bgcani2$clade2 & bgcani2$ani > 95, "withinspecies"] <- "zsal"

withinani <- bgcani2[bgcani2$withinspecies == "zsal",]
withinani <- withinani[complete.cases(withinani), ]

above95 <- withinani[, c("cluster1", "cluster2")]
write.table(above95, quote = F, row.names = F, file = "BGCs95.txt")

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

bgcani5 <- completeFun(bgcani2, "withinspecies")

library(plyr)
mu <- ddply(bgcani5, "withinspecies", summarise, grp.mean=mean(ani), grp.sd=sd(ani))
mu


library(ggplot2)

allbgcani <- ggplot(bgcani5, aes(x = ani, color = withinspecies, fill = withinspecies)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 200) +
  theme_classic() +
  geom_density(alpha = 0.6) +
  scale_color_manual(values = c("lightgray", "black", "green")) +
  scale_fill_manual(values = c("lightgray", "black", "green")) +
  geom_vline(xintercept = 95, color = "red", linetype = "dashed") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(80,100))



library(gridExtra) 

pdf("bgcANIvgenomeANI.pdf", width = 12, height = 10)
grid.arrange(WGani, allbgcani, ncol=1)
dev.off()






