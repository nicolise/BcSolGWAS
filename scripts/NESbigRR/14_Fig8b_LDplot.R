#Nicole E Soltis
#081117
#plot of SNPs along gene of interest

#12_singleGeneManhattan.R
#---------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")

#take the SNPs over the threshold for each phenotype

TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}


names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,-c(1)]
#only look at chromosome 16
HEM.plotdata <- HEM.plotdata[which(HEM.plotdata$Chrom=='16'),]
HEM.topSNPs <- HEM.plotdata
#get the start position of chromosome 16
min(HEM.topSNPs$Index)
max(HEM.topSNPs$Index)
max(HEM.topSNPs$Index) - min(HEM.topSNPs$Index)

#narrow window: +- 1 kb
HEM.topSNPs$Chr16Index <- HEM.topSNPs$Index - min(HEM.topSNPs$Index) + 1
min(HEM.topSNPs$Chr16Index)
#now get target region within chromosome 16
#my feature: about 1kb
#345785 to 346542
#and I'll add 2kb on each side
HEM.topSNPsSM <- HEM.topSNPs[which(HEM.topSNPs$Pos < 347542),]
HEM.topSNPsSM <- HEM.topSNPsSM[which(HEM.topSNPsSM$Pos > 344785),]

#trying a bigger window (8kb) to find missing phenos
HEM.topSNPsSM <- HEM.topSNPs[which(HEM.topSNPs$Pos < 355042),]
HEM.topSNPsSM <- HEM.topSNPsSM[which(HEM.topSNPsSM$Pos > 337285),]

#now to draw that darn LD plot
#this requires 3 file formats, none of which I have
library("snp.plotter")


#-------------------------------------------------------------------------
#START HERE, ACTUALLY
#this might make DNAbin files, from vcf files
#let's try just a 4kb region on Chromosome 16
#goes from POS 341 to POS 664510
#try POS 344785 to POS 347542
library("vcfR")
vcf <- read.vcfR("data/genome/chr16_analysis/chr16_analysis_seg.recode.vcf", verbose = FALSE)
my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, verbose = FALSE)
#this might make haplotype files, from DNAbin files
library("pegas")
