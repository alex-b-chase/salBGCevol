rm(list=ls())
setwd("/Volumes/Bennett_BACKUP/Research/jensen-metabolites/salinispora-BGCevol/metabolomics-data/")

features <- read.table("featurelist.txt", header = T, row.names = 1)

sptK <- 417.24058
sptHJ <- 431.25623
sptCGI <- 445.27188
sptBF <- 459.28753
sptA <- 473.30318

salA <- 314.1159
salB <- 280.1549
salC <- 284.1053
salD <- 266.1392
salE <- 294.1705
salI <- 328.1315
salJ <- 298.121
salK <- 252.1236

slcA <- 659.34764
slcBC <- 673.36329
slcDE <- 687.37894
slcFG <- 701.39459

lomA <- 683.264
lomB <- 539.694255
lomC <- 670.2738
lomD <- 677.2816
lomD <- 1353.55535
lomE <- 684.2894

lym <- 311.11439
lymol <- 299.11439
neolymA <- 339.14569
neolymB <- 325.13004
neolymC <- 353.16154
neolymolA <- 327.14569
neolymolB <- 313.13004
neolymolC <- 341.16134

spo <- 539.0956
slm <- 470.29062

rifS <- 696.30198
rifSNa <- 718.283398
rif8 <- 640.31215

sta7 <- 483.20321
sta <- 467.20829
staspo <- 468.21557
staoxo <- 497.18247

expfeature <- slcFG

ppmmax <- abs(10 * expfeature / 1e6 + expfeature)
ppmmin <- abs(10 * expfeature / 1e6 - expfeature)

features[features$mz >= ppmmin & features$mz <= ppmmax, ]

