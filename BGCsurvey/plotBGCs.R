setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/micromonospora_BGCsurvey/patric_genomes')

salbcgs <- read.table("../../salinispora-genomics/percBGCs.txt", header = T, sep = '\t')
bgcsurvey <- read.table('percBGCs.txt', header = T, sep = '\t')

library(ggplot2)
# plot number of BGCs per Sal species
salbcgs2 <- salbcgs[salbcgs$clade != "clade9", ]

mycolors <- c("#ff6a6a", "#00b2ee", "#adff2f", "#8b7500", "#b752ff",
              "#00ffff", "#ffa500", "#006400", "#dbcc3a", "black")

pdf('numBGCsXspecies.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

ggplot(data = salbcgs2, aes(x = reorder(clade, numBGCs, FUN = median), y = numBGCs)) +
  geom_boxplot(aes(fill = clade)) +
  theme_bw() +
  geom_hline(yintercept = mean(salbcgs2$numBGCs), linetype = "dashed", color = "black") +
  stat_summary(fun.y = mean, colour = "darkred", geom = "point", 
               shape = 18, size = 3) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "Number of BGCs/Genome", x = "Species") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_blank())

par(op)
dev.off()

# combine dataframes
library(data.table)
bgcstotal <- rbindlist(list(bgcsurvey, salbcgs))


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

genomeorder <- c("clade7",
                 "clade6",
                 "clade8",
                 "clade5",
                 "clade3",
                 "clade4",
                 "clade2",
                 "clade1",
                 "clade10",
                 "Jishengella",
                 "Verrucosispora",
                 "Xiangella",
                 "Micromonospora",
                 "Plantactinospora",
                 "Asanoa",
                 "Actinoplanes",
                 "Couchioplanes",
                 "Mycobacterium",
                 "Amycolatopsis",
                 "Pseudonocardia",
                 "Arthrobacter",
                 "Curtobacterium",
                 "Kitasatospora",
                 "Streptomyces",
                 "Frankia",
                 "Rubrobacter",
                 "Deinococcus",
                 "Fibrobacter",
                 "Bacteroides",
                 "Prevotella",
                 "Enterococcus",
                 "Streptococcus",
                 "Lactobacillus",
                 "Bacillus",
                 "Clostridium",
                 "Gloeobacter",
                 "Prochlorococcus",
                 "Synechococcus",
                 "Microcystis",
                 "Moorea",
                 "Nostoc",
                 "Sinorhizobium",
                 "Burkholderia",
                 "Vibrio",
                 "Escherichia",
                 "Salmonella",
                 "Thermotoga",
                 "Aquifex"
)

library(ggplot2)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right

summsaltot <- summarySE(salbcgs, measurevar="numBGCs", groupvars=c("clade"))
summsal <- summarySE(salbcgs, measurevar="percBGCs", groupvars=c("clade"))
################################################################################
# summarize by num of total BGCs
################################################################################
summBGC <- summarySE(bgcstotal, measurevar="numBGCs", groupvars=c("genus"))
summBGC$genus <- factor(summBGC$genus, levels = genomeorder)
summBGC <- summBGC[complete.cases(summBGC$genus), ]


pdf('numBGCsXgenus.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

p1 <- ggplot(summBGC, aes(x = genus, y = numBGCs, fill = "lightgray")) +
  stat_summary(fun.y = "mean", geom = "bar") +
  geom_errorbar(aes(ymin = numBGCs - sd, ymax = numBGCs + sd),
                width = 0.2,                    # Width of the error bars
                position = position_dodge(.9)) +
  geom_hline(yintercept = mean(salbcgs$numBGCs), linetype = "dashed", color = "black") +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "Number of BGCs", x = "Genus") +
  scale_colour_identity() +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=20))
p1

par(op)
dev.off()


################################################################################
# summarize by percent BGCs
################################################################################
summpercBGC <- summarySE(bgcstotal, measurevar="percBGCs", groupvars=c("genus"))
# order the genera
summpercBGC$genus <- factor(summpercBGC$genus, levels = genomeorder)
summpercBGC <- summpercBGC[complete.cases(summpercBGC$genus), ]

pdf('percBGCsXgenus.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

p2 <- ggplot(summpercBGC, aes(x = genus, y = percBGCs, group=1)) +
  geom_line(linetype = "dashed") +
  geom_errorbar(aes(ymin = percBGCs - sd, ymax = percBGCs + sd),
                width = 0.2,                    # Width of the error bars
                position = position_dodge(.9), color = "lightcoral") +
  geom_point(aes(size = 2, color = "lightcoral"), shape = 15, show.legend = F) +
  geom_hline(yintercept = mean(salbcgs$percBGCs), linetype = "dashed", color = "red") +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "Percent of BGCs/Genome", x = "Genus") +
  scale_colour_identity() +
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_blank())
p2

par(op)
dev.off()

library(gtable)
library(grid)
library(gridExtra)

# extract gtable
g1<-ggplot_gtable(ggplot_build(p1))
g2<-ggplot_gtable(ggplot_build(p2))

pp<-c(subset(g1$layout,name=="panel",se=t:r))
g<-gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]],pp$t,pp$l,pp$b,pp$l)

ia<-which(g2$layout$name=="axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

pdf("totalBGCpotential.pdf", height=12, width=16)
grid.draw(g)
dev.off()



### stats for percBGC
bgcstotal$tax <- ifelse(grepl("clade", bgcstotal$genus), "Salinispora", paste(bgcstotal$genus))
model <- aov(data = bgcstotal, percBGCs ~ tax)
summary(model)
sal <- as.data.frame(TukeyHSD(model)$tax)
sal$tax <- rownames(sal)
library(data.table)
salcomp <- sal[sal$tax %like% "Salinispora", ]
