rm(list=ls())
library(dplyr); library(ggplot2)
setwd("~/Projects/BcSolGWAS")
#read in annotation file
DomestAnt <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")

#30255 SNPs for any Domest trait

plants12 <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/12Plants_AllSNPsOver99_SegLong_trueMAF20_10NA.csv")

names(DomestAnt)
names(plants12)
DomestAnt <- DomestAnt[,c(2,3,4,11,12)]
plants12 <- plants12[,c(2,3,4,5,7)]
names(DomestAnt)[5] <- "Trait"
myblob <- rbind(plants12, DomestAnt)
myblob$Chr.Seg.Pos <- paste(myblob$Chrom, myblob$Segment, myblob$Pos, sep=".")
#myblob2 <- myblob[unique(myblob$Chr.Seg.Pos),]
myblob2 <- myblob[!duplicated(myblob[,'Chr.Seg.Pos']),]

