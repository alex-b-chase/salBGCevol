setwd('/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/BGCsurvey/')

salbcgs <- read.table("../BGCanalysis/percBGCs.txt", header = T, sep = '\t')
bgcsurvey <- read.table('patric_genomes/percBGCs.txt', header = T, sep = '\t')

library(ggplot2)

mycolors <- c("#00b2ee", "#8b7500", "#adff2f", "#ffa500", "#b752ff",
               "#00ffff", "#dbcc3a", "#ff6a6a", "#006400", "black")

pdf('../BGCanalysis/numBGCsXspecies.pdf', height=12, width=14)
op <- par(mar=c(5, 6, 4, 2) + 0.1)

ggplot(data = salbcgs, aes(x = reorder(clade, numBGCs, FUN = median), y = numBGCs)) +
  geom_hline(yintercept = mean(salbcgs$numBGCs), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(20,30), linetype = "dashed", color = "lightgray") +
  geom_boxplot(aes(fill = clade)) +
  theme_bw() +
  stat_summary(fun = mean, colour = "darkred", geom = "point", 
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

model <- aov(data = salbcgs, numBGCs ~ clade)
model2 <- aov(data = salbcgs, percBGCs ~ clade)
summary(model)
summary(model2)

# combine dataframes
library(data.table)
bgcstotal <- rbindlist(list(bgcsurvey, salbcgs), use.names=FALSE)


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

genomeorder <- c("Aquifex",
                 "Thermotoga",
                 "Sinorhizobium",
                 "Burkholderia",
                 "Vibrio",
                 "Escherichia",
                 "Salmonella",
                 "Fibrobacter",
                 "Prevotella",
                 "Bacteroides",
                 "Clostridium",
                 "Bacillus",
                 "Lactobacillus",
                 "Streptococcus",
                 "Enterococcus",
                 "Deinococcus",
                 "Gloeobacter",
                 "Prochlorococcus",
                 "Synechococcus",
                 "Microcystis",
                 "Nostoc",
                 "Moorea",
                 "Rubrobacter",
                 "Streptomyces",
                 "Kitasatospora",
                 "Arthrobacter",
                 "Curtobacterium",
                 "Frankia",
                 "Pseudonocardia",
                 "Mycobacterium",
                 "Amycolatopsis",
                 "Actinoplanes",
                 "Couchioplanes",
                 "Asanoa",
                 "Plantactinospora",
                 "Micromonospora",
                 "Jishengella",
                 "Xiangella",
                 "Verrucosispora",
                 "arenicola",
                 "tropica",
                 "fenicalii",
                 "oceanensis",
                 "cortesiana",
                 "mooreana",
                 "pacifica",
                 "goodfellowii",
                 "vitiensis"
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
  stat_summary(fun = "mean", geom = "bar") +
  geom_errorbar(aes(ymin = numBGCs - se, ymax = numBGCs + se),
                width = 0,                    # Width of the error bars
                position = position_dodge(.9)) +
  geom_hline(yintercept = mean(salbcgs$numBGCs), linetype = "dashed", color = "red") +
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
  geom_hline(yintercept = c(0,0.05,0.1,0.15,0.2), linetype = "dashed", color = "lightgray") +
  geom_line(linetype = "dashed") +
  geom_errorbar(aes(ymin = percBGCs - se, ymax = percBGCs + se),
                width = 0,                    # Width of the error bars
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

