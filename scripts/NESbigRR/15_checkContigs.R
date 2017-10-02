#Nicole E Soltis
#092917
#script to check the length of our contigs for T4 SNP set
#15_checkContigs.R
#-------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

HEM.plotdata.og <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")

HEM.plotdata <- HEM.plotdata.og
HEM.plotdata <- HEM.plotdata[,-c(1)]

names(HEM.plotdata)

HEM.plotdata$Chrom.Seg <- paste(HEM.plotdata$Chrom, HEM.plotdata$Segment, sep=".")
unique(HEM.plotdata$Chrom.Seg)
min(HEM.plotdata$Pos, )

mycontigSTART <- HEM.plotdata[ HEM.plotdata$Pos == ave(HEM.plotdata$Pos, HEM.plotdata$Chrom.Seg, FUN=min), ]
mycontigSTART <- mycontigSTART[,c(1,2,17,3,14)]
mycontigSTOP <- HEM.plotdata[ HEM.plotdata$Pos == ave(HEM.plotdata$Pos, HEM.plotdata$Chrom.Seg, FUN=max), ]
names(mycontigSTOP)[3] <- "StopPos"
mycontigSTOP <- mycontigSTOP[,c(17,3)]
mycontigs <- merge(mycontigSTART,mycontigSTOP, by="Chrom.Seg")
mycontigs$Length <- (mycontigs$StopPos - mycontigs$Pos) / 1000000
