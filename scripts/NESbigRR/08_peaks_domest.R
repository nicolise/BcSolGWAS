#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/domestication/Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/domestication/Sl_DomesticationLS_MAF20.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]

#take the top 50 over the threshold for each phenotype
library(plyr)

names(HEM.plotdata)

TH999 <- HEM.thresh[4,]
for (i in 2:ncol(TH999)){
  assign(paste("TH999_", names(TH999[i]), sep=""),as.numeric(TH999[i]))
}

names(HEM.plotdata)

#All groups
assign(paste("HEM.", names(HEM.plotdata[4]), sep=""), subset(HEM.plotdata, HEM.plotdata[4] > get(paste("TH999_", names(HEM.plotdata[4]), sep="")),
                                                             select=c(Chrom,Segment,Pos,4)))

#or each by hand
# HEM.Domesticated <- subset(HEM.plotdata, Domesticated > 6.285562e-04, 
#                                 select=c(Chrom,Pos,LA0410))
# HEM.Domesticated <- rename(HEM.LA0410, c("LA0410" = "Effect"))
# HEM.LA0410$Plant <- "LA0410"
# HEM.LA0410 <- head(arrange(HEM.LA0410,desc(Effect)), n = 50)

#then combine
HEM.Domesticated <- rename(HEM.Domesticated, c("Domesticated" ="Effect"))
HEM.Domesticated$Trait <- "Domesticated"
HEM.Wild <- rename(HEM.Wild, c("Wild" ="Effect"))
HEM.Wild$Trait <- "Wild"
HEM.DmWoD <- rename(HEM.DmWoD, c("DmWoD" ="Effect"))
HEM.DmWoD$Trait <- "DmWoD"

Top50SNP <- rbind(HEM.Domesticated, HEM.Wild, HEM.DmWoD)
write.csv(Top50SNP, "data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest.csv")

#max pos is 1001108

Top50SNP$Plot <- (Top50SNP$Chrom*1000000 + Top50SNP$Pos)

Top50SNP$Chrom <- gsub("Chromosome", "", Top50SNP$Chrom)
Top50SNP$Chrom <- as.numeric(as.character(Top50SNP$Chrom))
Top50SNP$Pos <- as.numeric(as.character(Top50SNP$Pos))

#sort dataframe rows in order of Chrom, then Pos
Top50SNP <- Top50SNP[with(Top50SNP, order(Chrom, Pos)), ]

#Make plotting variables
Top50SNP$Index = NA
ticks = NULL
lastbase = 0

#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(Top50SNP$Chrom)) {
  print(i)
  if (i==1) {
    Top50SNP[Top50SNP$Chrom==i, ]$Index=Top50SNP[Top50SNP$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=+lastbase+max(subset(Top50SNP,Top50SNP$Chrom==i-1)$Pos, 1)
    Top50SNP[Top50SNP$Chrom==i, ]$Index=Top50SNP[Top50SNP$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, Top50SNP[Top50SNP$Chrom==i, ]$Index[floor(length(Top50SNP[Top50SNP$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(Top50SNP$Index),max(Top50SNP$Index))

library(ggplot2)
plot1 <- ggplot(Top50SNP, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Plant)))+
  theme_bw()
 
  write.csv(Top50SNP, "results/Domestication_Top50SNP.csv")
  