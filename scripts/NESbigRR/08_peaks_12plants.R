#Nicole E Soltis 
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input File: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv and .Thresh.csv
#Output File: results/Domestication_TopSNPs_SegLong.csv, results/Domestication_TopSNPs_SegWide.csv
#this goes into gene annotation and then venn diagrams
#Plots: NONE
############################################################################

#Load packages
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.Thresh.csv")

#take the top 50 over the threshold for each phenotype
TH95pos <- HEM.thresh[1,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}
TH95neg <- HEM.thresh[5,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}
TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}
TH999pos <- HEM.thresh[4,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}
TH999neg <- HEM.thresh[8,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,-c(1)]
#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in 4:15){
assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}

subset(HEM.plotdata, HEM.plotdata[,16] < -0.0002432714)

#for wide format plot need fewer SNPs -- take top 1000 for each? but all of DmWoD
#top 1000
#4:15
#basic version, no loop:
#HEMpos.Domesticated <- head(arrange(HEMpos.Domesticated, desc(Domesticated)), n=500)
#design the loop:
#assign(paste("HEMpos.", names(HEM.plotdata[4]), sep=""), head(arrange(HEMpos.Domesticated, desc(Domesticated)), n=500))
#HEMpos.Domesticated == get(paste("HEMpos.", names(HEM.plotdata[4]), sep=""))
#Domesticated == get(paste("HEMpos.", names(HEM.plotdata[4]), sep=""))[,5]

#test script
#assign(paste("HEMpos.", names(HEM.plotdata[4]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[4]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[4]), sep=""))[,5])), n=500))

#for top 500 only
for (i in 4:15){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
}
for (i in 4:15){
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), tail(arrange(get(paste("HEMneg.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
}

#combine pos and neg by group
for (i in 4:15){
assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
}

#then combine
assign(dfname, rename(dfname, c("Trait" = "Effect")))

mydfname <- paste("HEM.", names(HEM.plotdata[5]), sep="")
assign(mydfname, rename(get(mydfname), c("Trait" = "Effect")))

#Trait? names(get(mydfname))[5]
HEM.Domesticated <- rename(HEM.Domesticated, c("Domesticated" ="Effect"))
HEM.Domesticated$Trait <- "Domesticated"
HEM.Wild <- rename(HEM.Wild, c("Wild" ="Effect"))
HEM.Wild$Trait <- "Wild"
HEM.DmWoD <- rename(HEM.DmWoD, c("DmWoD" ="Effect"))
HEM.DmWoD$Trait <- "DmWoD"

HEM.topSNPs <- rbind(HEM.Domesticated, HEM.Wild, HEM.DmWoD)

library(ggplot2)
plot1 <- ggplot(HEM.plotdata, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Trait)))+
  theme_bw()
  
#make it wide format
#currently long format : Chrom, Segment, Pos, Index, Effect, Trait
write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")
write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1000SNPs_SegLong_trueMAF20_10NA.csv")

TopSNP.wide.DM <- reshape(HEM.topSNPs, 
                         timevar = "Trait",
                         idvar = c("Chrom","Segment","Pos","Index"),
                         direction = "wide")

write.csv(TopSNP.wide.DM, "results/Domestication_Top1000SNPs_SegWide_trueMAF20_10NA.csv")
