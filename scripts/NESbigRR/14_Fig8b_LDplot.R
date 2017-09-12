#Nicole E Soltis
#081117
#plot of SNPs along gene of interest

#12_singleGeneManhattan.R
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
HEM.plotdata <- HEM.plotdata[which(HEM.plotdata$Chrom=='16'),]

mySNPlist <- read.csv("data/genome/chr16_analysis/SNPlistfromPED.csv")
mySNPlist <- mySNPlist[,2]
HEM.plotdata <- HEM.plotdata[HEM.plotdata$Pos %in% mySNPlist, ]
#there is some duplication here due to Segment ambiguity. I am going to assume that Segment 11 is the worst contig and on down, and as such only keep the first mention of each Chr.Pos
HEM.plotdata <- HEM.plotdata[with(HEM.plotdata, order(Segment)), ]
HEM.plotdata <- HEM.plotdata[!duplicated(HEM.plotdata$Pos), ]

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

#now convert to SNP.FILE format
#try it just for the first plant genotype (LA1547), for now
mySNP.FILE <- HEM.topSNPsSM
mySNP.FILE$ASSOC <- ifelse(mySNP.FILE$LA1547 > 0, '+', '-')
mySNP.FILE$SNP.NAME <- paste("mysnp", mySNP.FILE$Pos, sep="_")
mySNP.FILE$LOC <- mySNP.FILE$Pos
mySNP.FILE$SS.PVAL <- mySNP.FILE$LA1547

mySNP.FILE <- mySNP.FILE[,c(21:24)]

#-------------------------------------------------------------------------
#START HERE, ACTUALLY
#now to draw that darn LD plot
#this requires 3 file formats, none of which I have
library("snp.plotter")

#this might make DNAbin files, from vcf files
#let's try just a 4kb region on Chromosome 16
#goes from POS 341 to POS 664510
#try POS 344785 to POS 347542
library("vcfR")
vcf <- read.vcfR("data/genome/chr16_analysis/chr16_analysis_seg.recode.vcf", verbose = FALSE)
my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, verbose = FALSE)
my_dnabin1
#plot it by site
ape::image.DNAbin(my_dnabin1[,ape::seg.sites(my_dnabin1)])
#save it as a fasta
write.dna( my_dnabin1, file = 'data/genome/chr16_analysis/chr16_analysis_seg.fasta', format = 'fasta' )
#this might make haplotype files, from DNAbin files
#pegas fails on ubuntu -- run on windows PC
library("pegas")
my_haplo <- haplotype(my_dnabin1)
plot(sort(my_haplo))
attr(my_haplo, "index")[[2]]
print(my_haplo)
my_haplo_dist <- dist.haplotype.loci(my_haplo)
?haplotype()

subset(my_haplo, minfreq=20)
as.data.frame(my_haplo)
