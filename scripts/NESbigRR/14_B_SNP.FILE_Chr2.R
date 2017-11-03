#Nicole E Soltis
#091317
#plot of SNPs along gene of interest
#and now: haplotype plots

#14_B_SNP.FILE.R
#---------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF20_20NA/SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF20_20NA/SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")

#take the SNPs over the threshold for each phenotype

#TH99pos <- HEM.thresh[3,]
#for (i in 2:ncol(TH99pos)){
#  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
#}
#TH99neg <- HEM.thresh[7,]
#for (i in 2:ncol(TH99neg)){
#  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
#}

names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,-c(1)]
#only look at chromosome 16
#chrom and segment are separate already here
HEM.plotdata <- HEM.plotdata[which(HEM.plotdata$Chrom.Seg=='2.2'),]

mySNPlist <- read.csv("data/genome/chr2_analysis/SNPlistfromPED.csv")
mySNPlist <- mySNPlist[,2]
HEM.plotdata <- HEM.plotdata[HEM.plotdata$Pos %in% mySNPlist, ]
HEM.plotdata <- HEM.plotdata[!duplicated(HEM.plotdata$Pos), ] #good

HEM.topSNPs <- HEM.plotdata
#get the start position of chromosome 2
min(HEM.topSNPs$Index)
max(HEM.topSNPs$Index)
max(HEM.topSNPs$Index) - min(HEM.topSNPs$Index)

##EDIT HEREEEE
#narrow window: +- 1 kb
HEM.topSNPs$Chr2Index <- HEM.topSNPs$Index - min(HEM.topSNPs$Index) + 1
min(HEM.topSNPs$Chr2Index)

#NOW filter to just region of interest on Chromosome 2
#823306 to 828345
#and this includes 2kb on each side
HEM.topSNPsSM <- HEM.topSNPs[which(HEM.topSNPs$Pos < 828346),]
HEM.topSNPsSM <- HEM.topSNPsSM[which(HEM.topSNPsSM$Pos > 823305),]

#now convert to SNP.FILE format
#try it just for the first plant genotype (LA1547), for now
mySNP.FILE <- HEM.topSNPsSM
mySNP.FILE$ASSOC <- ifelse(mySNP.FILE$LA1547 > 0, '+', '-')
mySNP.FILE$SNP.NAME <- paste("mysnp", mySNP.FILE$Pos, sep="_")
mySNP.FILE$LOC <- mySNP.FILE$Pos
mySNP.FILE$SS.PVAL <- mySNP.FILE$LA1547

mySNP.FILE <- mySNP.FILE[,c(21:24)]

#write it
write.table(mySNP.FILE, file="data/genome/chr2_analysis/SNP.FILE.tsv",quote=FALSE, sep='\t',col.names=TRUE,row.names=FALSE)

write.table(mySNP.FILE, file="data/genome/chr2_analysis/SNP.FILE.txt",quote=FALSE, sep=' ',col.names=TRUE,row.names=FALSE)

