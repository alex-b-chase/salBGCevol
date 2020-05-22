setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCanalysis/evolprocesses/conservedBGCgenes/')

######## ######## ######## ######## ######## ######## ######## 
####### calculate dNdS ratios
######## ######## ######## ######## ######## ######## ######## 

# run screen genes before to eliminate genes with "bad" reading frames (i.e. not multiple of 3bp)

library(seqinr)
library(Biostrings)
library(ggplot2)

genes <- list.files("goodalignments/")
X <- read.alignment(file = paste("goodalignments/", genes[1], sep=""), format = "fasta")
Y <- matrix(nrow = 0, ncol = (X$nb - 1))
dN <- rep(0, length(genes))
dS <- rep(0, length(genes))
vdN <- rep(0, length(genes))
vdS <- rep(0, length(genes))
z <- rep(0, length(genes))
for(i in 1:length(genes)){
  X <- read.alignment(file = paste("goodalignments/", genes[i], sep=""), format = "fasta")
  #Y <- rbind(Y, kaks(X)$ka[1:(X$nb - 1)])
  dN[i] <- kaks(X)$ka[1]
  dS[i] <- kaks(X)$ks[1]
  vdN[i] <- kaks(X)$vka[1]
  vdS[i] <- kaks(X)$vks[1]
  z[i] <- (dN[i] - dS[i])/(sqrt(vdN[i] + vdS[i]))
}
p <- 2*(1-pnorm(abs(z)))
padj <- p.adjust(p, method = "BH")
lp <- -log10(padj)
plot(1:length(p), -log10(padj))

# remove the infinities
dS[dS < 0] <- 0
dNdS <- dN/dS
dNdS[!is.finite(dNdS)] <- 0

# Make table
genelist <- sapply(strsplit((list.files("goodalignments/")), split="_"), function(x) x[[1]])

Dtable <- data.frame(genelist,  dN, dS, dNdS, z, p, padj)
colnames(Dtable) <- c("gene", "dN", "dS", "dNdS", "Z", "P", "Adjusted.P")

desc <- read.csv(file = "goodconservedBGCgenes.txt", sep = '\t', header = T)
desc$gene <- paste("total.", desc$bgcname, ".",desc$cluster, ".aln", sep = '')

z <- merge(Dtable, desc, by = "gene")

pdf("dnds-frequencies.pdf", height = 10, width = 10)

q1 <- ggplot(data=z, aes(log2(z$dNdS))) + 
  geom_histogram(breaks=seq(-8, 1, by = 0.5), 
                 col="black", 
                 fill="lightgray") + 
  labs(title="Directional Selection") +
  labs(x="dN/dS", y="Count") +
  theme_bw()

q1
dev.off()

write.csv(z, file = "dnds.allgenes.csv", qmethod='escape', quote=TRUE, row.names = F)


# ######## ######## ######## ######## ######## ######## ######## 
# ######## performs most population genomic analyses
# ######## ######## ######## ######## ######## ######## ######## 
 
library(PopGenome)

GENOME.class <- readData("goodalignments/")
get.sum.data(GENOME.class)
get.individuals(GENOME.class)

# Run necessary module
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class)
GENOME.class@n.sites
pi <- GENOME.class@Pi
td <- GENOME.class@Tajima.D
genes <- sapply(strsplit((list.files("goodalignments/")), split="_"), function(x) x[[1]])

Dtable <- data.frame(genes,  pi, td)
colnames(Dtable) <- c("gene", "pi", "Tajima.D")


genetable <- merge(Dtable, desc, by = 'gene')
aggregate(genetable$pi, by=list(genetable$bgcname), FUN=mean)


# set populations
clade1 <- c("StCNB440BA",
             "StCNB476BA",
             "StCNB536BA",
             "StCNH898BA",
             "StCNR699BA",
             "StCNS197BA",
             "StCNS416BA",
             "StCNT250BA",
             "StCNT261BA",
             "StCNY012BA",
             "StCNY678YU",
             "StCNY681YU")
clade10 <- c("SaCNB458BA",
               "SaCNB527BA",
               "SaCNH643BA",
               "SaCNH646BA",
               "SaCNH713RS",
               "SaCNH718RS",
               "SaCNH877BA",
               "SaCNH905BA",
               "SaCNH941SC",
               "SaCNH962SC",
               "SaCNH963SC",
               "SaCNH964SC",
               "SaCNH996SC",
               "SaCNH996BSC",
               "SaCNP105SC",
               "SaCNP193SC",
               "SaCNQ748GU",
               "SaCNQ884GU",
               "SaCNR107GU",
               "SaCNR425GU",
               "SaCNR921PL",
               "SaCNS051PL",
               "SaCNS205PL",
               "SaCNS243PL",
               "SaCNS296PL",
               "SaCNS299PL",
               "SaCNS325PL",
               "SaCNS342FJ",
               "SaCNS673FJ",
               "SaCNS744FJ",
               "SaCNS820FJ",
               "SaCNS848FJ",
               "SaCNT005FJ",
               "SaCNT798HI",
               "SaCNT799HI",
               "SaCNT800HI",
               "SaCNT849HI",
               "SaCNT850HI",
               "SaCNT857HI",
               "SaCNT859HI",
               "SaCNX481PM",
               "SaCNX482PM",
               "SaCNX508PM",
               "SaCNX814PM",
               "SaCNX891PM",
               "SaCNY011BA",
               "SaCNY230FJ",
               "SaCNY231FJ",
               "SaCNY234FJ",
               "SaCNY237FJ",
               "SaCNY244FJ",
               "SaCNY256FJ",
               "SaCNY260FJ",
               "SaCNY280FJ",
               "SaCNY282FJ",
               "SaCNY486PV",
               "SaCNY679YU",
               "SaCNY685YU",
               "SaCNY690YU",
               "SaCNY694YU",
               "SaDSM45545FJ")
clade2 <- c("SpCNR942PL", "SpCNT569FJ")
clade3 <- c("SpCNY202SC")
clade4 <- c("SpCNS237PL",
              "SpCNY646RS",
              "SpDSM45549FJ") 
clade5 <- c("SpCNS860FJ",
                "SpCNS996FJ",
                "SpCNT029FJ",
                "SpCNT045FJ",
                "SpCNT124FJ",
                "SpCNT403FJ",
                "SpCNT584FJ",
                "SpCNT854HI",
                "SpCNY673MD",
                "SpCNY703MD",
                "SpDSM45543FJ",
                "SpDSM45547FJ")
clade6 <- c("SpCNY666MD")
clade7 <- c("SpCNS055PL",
               "SpCNS801FJ",
               "SpDSM45548FJ")
clade8 <- c("SpCNH732RS",
              "SpCNQ768GU",
              "SpCNR114GU",
              "SpCNR510GU",
              "SpCNR894PL",
              "SpCNR909PL",
              "SpCNS103PL",
              "SpCNT001FJ",
              "SpCNT003FJ",
              "SpCNT084FJ",
              "SpCNT131FJ",
              "SpCNT133FJ",
              "SpCNT603FJ",
              "SpCNT609FJ",
              "SpCNT796HI",
              "SpCNT851HI",
              "SpCNT855HI",
              "SpCNY239FJ",
              "SpCNY330SC",
              "SpCNY331SC",
              "SpCNY363SC",
              "SpCNY498PV",
              "SpDSM45544FJ")


pops <- list(clade10, clade1, clade2, clade3, 
             clade4, clade5, clade6, clade7, clade8)

GENOME.class <- set.populations(GENOME.class, pops)

GENOME.class <- MKT(GENOME.class)

GENOME.class <- neutrality.stats(GENOME.class)
GENOME.class <- F_ST.stats(GENOME.class)

get.neutrality(GENOME.class)[[1]]
GENOME.class@nuc.diversity.within

pop.pi <- as.data.frame(GENOME.class@Pi)
pop.td <- as.data.frame(GENOME.class@Tajima.D)

pop.pi$gene <- rownames(pop.pi)
pop.td$gene <- rownames(pop.td)

colnames(pop.pi) <- c("clade10", "clade1", "clade2", "clade3", "clade4", "clade5", "clade6", "clade7", "clade8", "gene")

library(reshape2)

m.dtable <- melt(pop.pi)
gene.pop.pi <- merge(m.dtable, desc, by = 'gene')

aggregate(value ~ variable + bgcname, data = gene.pop.pi, function(x) c(mean = mean(x), sd = sd(x)))
write.table(gene.pop.pi, file = "popstats.pi.txt", row.names = F, sep = '\t', quote = F)

# combine with core genome stats
corepi <- read.table('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/genomic-data/coregenes/goodcore/popstats.pi.txt', header = T, sep = '\t')
library(reshape2)
library(ggplot2)
m.corepi <- melt(corepi)
m.corepi$bgcname <- "core"

gene.pop.pi <- read.table('popstats.pi.txt', sep = '\t', header = T)
bgcpi <- gene.pop.pi[, c(1,2,3,4)]

# compare the BGCs to core
bgcpi$bgc <- "yes"
m.corepi$bgc <- "no"

totalcomp <- rbind(bgcpi, m.corepi)

# remove clades with only 1 genome
totalcomp2 <- totalcomp[totalcomp$variable != "clade9",]
totalcomp3 <- totalcomp2[totalcomp2$variable != "clade6",]
totalcomp4 <- totalcomp3[totalcomp3$variable != "clade3",]
# need to remove BGC that do not exist in species
# ST does not have slc, rif, sta
totalcomp4$bgcclade <- paste(totalcomp4$variable, "_", totalcomp4$bgcname, sep = '')
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade1_salinichelin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade1_staurosporine",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade1_rifamycin",]
# SA does not have lom, slm,  spo
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade10_lomaiviticin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade10_sporolide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade10_salinilactam",]
# SP(clade8) does not have sta, rif
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade8_staurosporine",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade8_rifamycin",]
# SP(clade5) does not have lom, slc, spo, rif, slm
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade5_lomaiviticin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade5_salinichelin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade5_sporolide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade5_rifamycin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade5_salinilactam",]
#SP(clade2) does not have slc, spo, sta, spt, rif, slm
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade2_salinichelin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade2_sporolide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade2_staurosporine",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade2_salinipostin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade2_rifamycin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade2_salinilactam",]
#SP(clade4) does not have lom, slc, spo, rif, slm
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade4_lomaiviticin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade4_salinichelin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade4_sporolide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade4_rifamycin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade4_salinilactam",]
#SP(clade7) does not have sal, slc, spo, sta, rif
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_salinosporamide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_salinichelin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_sporolide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_staurosporine",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_rifamycin",]

# remove all BGCs with only one strain because nt == 0
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_salinipostin",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade4_salinosporamide",]
totalcomp4 <- totalcomp4[totalcomp4$bgcclade != "clade7_salinilactam",]



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

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right

################################################################################
# summarize by num of total BGCs
################################################################################
clade1BGCs <- totalcomp4[totalcomp4$variable == "clade1", ]
core1 <- mean(corepi$clade1)
core1sd <- sd(corepi$clade1)
clade1BGCs$pmean <- log2(2^(clade1BGCs$value - core1))
clade1BGCs$zscore <- (clade1BGCs$value - core1) / core1sd

clade10BGCs <- totalcomp4[totalcomp4$variable == "clade10", ]
core10 <- mean(corepi$clade10)
core10sd <- sd(corepi$clade10)
clade10BGCs$pmean <- log2(2^(clade10BGCs$value - core10))
clade10BGCs$zscore <- (clade10BGCs$value - core10) / core10sd

clade2BGCs <- totalcomp4[totalcomp4$variable == "clade2", ]
core2 <- mean(corepi$clade2)
core2sd <- sd(corepi$clade2)
clade2BGCs$pmean <- log2(2^(clade2BGCs$value - core2))
clade2BGCs$zscore <- (clade2BGCs$value - core2) / core2sd

clade4BGCs <- totalcomp4[totalcomp4$variable == "clade4", ]
core4 <- mean(corepi$clade4)
core4sd <- sd(corepi$clade4)
clade4BGCs$pmean <- log2(2^(clade4BGCs$value - core4))
clade4BGCs$zscore <- (clade4BGCs$value - core4) / core4sd

clade5BGCs <- totalcomp4[totalcomp4$variable == "clade5", ]
core5 <- mean(corepi$clade5)
core5sd <- sd(corepi$clade5)
clade5BGCs$pmean <- log2(2^(clade5BGCs$value - core5))
clade5BGCs$zscore <- (clade5BGCs$value - core5) / core5sd

clade7BGCs <- totalcomp4[totalcomp4$variable == "clade7", ]
core7 <- mean(corepi$clade7)
core7sd <- sd(corepi$clade7)
clade7BGCs$pmean <- log2(2^(clade7BGCs$value - core7))
clade7BGCs$zscore <- (clade7BGCs$value - core7) / core7sd

clade8BGCs <- totalcomp4[totalcomp4$variable == "clade8", ]
core8 <- mean(corepi$clade8)
core8sd <- sd(corepi$clade8)
clade8BGCs$pmean <- log2(2^(clade8BGCs$value - core8))
clade8BGCs$zscore <- (clade8BGCs$value - core8) / core8sd

cladestotal <- rbind(clade1BGCs, clade10BGCs, clade2BGCs, clade4BGCs, clade5BGCs, clade7BGCs, clade8BGCs)


summBGC <- summarySE(cladestotal, measurevar="pmean", groupvars=c("bgcname", "variable"))
summBGCZ <- summarySE(cladestotal, measurevar="zscore", groupvars=c("bgcname", "variable"))
library(tidyverse)
summBGC2 = summBGC %>%
  complete(bgcname, variable, fill = list(pmean = 0))

summBGC2 <- summBGC2[summBGC2$variable != "clade6",]
summBGC2 <- summBGC2[summBGC2$variable != "clade9",]
summBGC2 <- summBGC2[summBGC2$variable != "clade3",]


summBGCZ2 = summBGCZ %>%
  complete(bgcname, variable, fill = list(zscore = 0))

summBGCZ2 <- summBGCZ2[summBGCZ2$variable != "clade6",]
summBGCZ2 <- summBGCZ2[summBGCZ2$variable != "clade9",]
summBGCZ2 <- summBGCZ2[summBGCZ2$variable != "clade3",]



pdf('nt.pi-diversity-log2.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

ggplot(summBGC2, aes(x = bgcname, y = pmean, fill = variable)) +
  geom_bar(position = position_dodge(.9), colour = "black", stat = "identity", width = 0.9) +
  geom_errorbar(position = position_dodge(.9), width = .25, aes(ymin = pmean - se, ymax = pmean + se)) +
  scale_fill_manual(values = c("#FF6A6A","#00B2EE", "#ADFF2F", "#B752FF", "#00FFFF", "#006400", "#DBCC3A", "#CCCCCC", "#FFFFFF", "#FFFFFF")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "log2 Fold Change [Nucleotide Diversity]", x = "BGC") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=20))

par(op)
dev.off()

pdf('nt.pi-diversity-zscore.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

ggplot(summBGCZ2, aes(x = bgcname, y = zscore, fill = variable)) +
  geom_bar(position = position_dodge(.9), colour = "black", stat = "identity", width = 0.9) +
  geom_errorbar(position = position_dodge(.9), width = .25, aes(ymin = zscore - se, ymax = zscore + se)) +
  scale_fill_manual(values = c("#FF6A6A","#00B2EE", "#ADFF2F", "#B752FF", "#00FFFF", "#006400", "#DBCC3A", "#CCCCCC", "#FFFFFF", "#FFFFFF")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "zScore [Nucleotide Diversity]", x = "BGC") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=20))

par(op)
dev.off()




summBGCcomp <- summarySE(totalcomp4, measurevar = "value", groupvars = c("bgc", "variable"))

pdf('nt.pi-diversityBGC.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

ggplot(summBGCcomp, aes(x = variable, y = value, fill = bgc)) +
  geom_bar(position = position_dodge(.9), colour = "black", stat = "identity") +
  geom_errorbar(position = position_dodge(.9), width = .25, aes(ymin = value - se, ymax = value + se)) +
  scale_fill_manual(values = c("#FFFFFF", "#CCCCCC")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "Nucleotide Diversity", x = "BGC") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=20))

par(op)
dev.off()







