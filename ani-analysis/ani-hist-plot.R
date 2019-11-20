rm(list=ls())
setwd("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-genomics/ani-analysis")

ani <- read.table('total.txt', sep = '\t', header = F)

hist(ani$V3, breaks = 200)

genomeid <- read.table("../coregenes/genomeID.txt", sep = '\t', header = F)
colnames(genomeid) <- c("strain1", "clade1")
genomeid2 <- genomeid
colnames(genomeid2) <- c("strain2", "clade2")

colnames(ani) <- c("strain1", "strain2", "ani")

ani1 <- merge(ani, genomeid, by = "strain1")
ani2 <- merge(ani1, genomeid2, by = "strain2")

ani3 <- ani2[ani2$strain2 != "SaCNY281_FJ",]
ani4 <- ani3[ani3$strain1 != "SaCNY281_FJ",]



ani4$species1 <- ifelse(grepl("Sa", ani4$strain1, ignore.case = T), "arenicola", 
                      ifelse(grepl("Sp", ani4$strain1, ignore.case = T), "pacifica", 
                             ifelse(grepl("St", ani4$strain1, ignore.case = T), "tropica",
                                    "Other")))
ani4$species2 <- ifelse(grepl("Sa", ani4$strain2, ignore.case = T), "arenicola", 
                        ifelse(grepl("Sp", ani4$strain2, ignore.case = T), "pacifica", 
                               ifelse(grepl("St", ani4$strain2, ignore.case = T), "tropica",
                                      "Other")))


ani4$withinclade <- NA
ani4[ani4$clade1 != ani4$clade2 & ani4$species1 != ani4$species2, "withinclade"] <- "no"
ani4[ani4$clade1 == ani4$clade2 & ani4$species1 == ani4$species2, "withinclade"] <- "yes"
ani4[ani4$clade1 != ani4$clade2 & ani4$species1 == ani4$species2, "withinclade"] <- "intra"

library(plyr)
mu <- ddply(ani4, "withinclade", summarise, grp.mean=mean(ani))



library(ggplot2)

pdf("ani-hist.pdf", height = 8, width = 6)

ggplot(ani4, aes(x = ani, color = withinclade, fill = withinclade)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 200) +
  theme_classic() +
  geom_density(alpha = 0.6) 

dev.off()