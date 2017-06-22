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
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.Thresh.csv")

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

HEM.plotdata <- HEM.plotdata[,-c(1)]

names(HEM.plotdata)

#conditionally replace values < threshold with zero
HEM.plotdata$Domesticated[HEM.plotdata$Domesticated < TH99pos_Domesticated & HEM.plotdata$Domesticated > 0] <- 0
HEM.plotdata$Wild[HEM.plotdata$Wild < TH99pos_Wild & HEM.plotdata$Wild > 0] <- 0
HEM.plotdata$DmWoD[HEM.plotdata$DmWoD < TH99pos_DmWoD & HEM.plotdata$DmWoD > 0] <- 0
HEM.plotdata$Domesticated[HEM.plotdata$Domesticated > TH99neg_Domesticated & HEM.plotdata$Domesticated < 0] <- 0
HEM.plotdata$Wild[HEM.plotdata$Wild > TH99neg_Wild & HEM.plotdata$Wild < 0] <- 0
HEM.plotdata$DmWoD[HEM.plotdata$DmWoD > TH99neg_DmWoD & HEM.plotdata$DmWoD < 0] <- 0
#remove rows if all 3 = 0 
HEM.plotdata <- HEM.plotdata[!(HEM.plotdata$Domesticated==0 & HEM.plotdata$Wild==0 & HEM.plotdata$DmWoD==0),]
#now add counting variable
HEM.plotdata$TotTraits <- ifelse(abs(HEM.plotdata$Domesticated) >0 & abs(HEM.plotdata$Wild) >0 & abs(HEM.plotdata$DmWoD) >0, "ALL", 
                              ifelse(abs(HEM.plotdata$Domesticated) >0 & abs(HEM.plotdata$Wild) >0, "DW",
                                     ifelse(abs(HEM.plotdata$DmWoD) >0 & abs(HEM.plotdata$Wild) >0, "WS",
                                            ifelse(abs(HEM.plotdata$Domesticated) >0 & abs(HEM.plotdata$DmWoD) >0, "DS",
                                                   ifelse(abs(HEM.plotdata$Domesticated) >0, "D",
                                                          ifelse(abs(HEM.plotdata$Wild) >0, "W", "S"))))))

table(HEM.plotdata$TotTraits)
#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in 4:6){
assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}

#for wide format plot need fewer SNPs -- take top 1000 for each? but all of DmWoD
#top 1000
HEMpos.Domesticated <- head(arrange(HEMpos.Domesticated, desc(Domesticated)), n=500)
HEMneg.Domesticated <- tail(arrange(HEMneg.Domesticated, desc(Domesticated)), n=500)
HEMpos.Wild <- head(arrange(HEMpos.Wild, desc(Wild)), n=500)
HEMneg.Wild <- tail(arrange(HEMneg.Wild, desc(Wild)), n=500)
#combine pos and neg by group
HEM.DmWoD <- rbind(HEMpos.DmWoD, HEMneg.DmWoD)
HEM.Domesticated <- rbind(HEMpos.Domesticated, HEMneg.Domesticated)
HEM.Wild <- rbind(HEMpos.Wild, HEMneg.Wild)


#then combine
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
write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1kSNPs_SegLong_trueMAF20_10NA.csv")

TopSNP.wide.DM <- reshape(HEM.topSNPs, 
                         timevar = "Trait",
                         idvar = c("Chrom","Segment","Pos","Index"),
                         direction = "wide")

write.csv(TopSNP.wide.DM, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1kSNPs_SegWide_trueMAF20_10NA.csv")
