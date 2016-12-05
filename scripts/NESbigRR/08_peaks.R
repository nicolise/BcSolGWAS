#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/NewModel0711/Sl_LesionSize_MAF20.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/NewModel0711/Sl_LesionSize_MAF20.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]

#take the top 50 over the threshold for each phenotype
library(plyr)

TH99 <- HEM.thresh[3,]
for (i in 2:ncol(TH99)){
  assign(paste("TH99_", names(TH99[i]), sep=""),as.numeric(TH99[i]))
}

TH999 <- HEM.thresh[4,]
for (i in 2:ncol(TH999)){
  assign(paste("TH999_", names(TH999[i]), sep=""),as.numeric(TH999[i]))
}

TH95 <- HEM.thresh[1,]
for (i in 2:ncol(TH95)){
  assign(paste("TH95_", names(TH95[i]), sep=""),as.numeric(TH95[i]))
}

names(HEM.plotdata)

#LA0410: 0.999
HEM.LA0410 <- subset(HEM.plotdata, LA410 > get(paste("TH999_", "LA410", sep="")), 
                                select=c(Chrom,Segment, Pos,LA410))
HEM.LA0410 <- rename(HEM.LA0410, c("LA410" = "Effect"))
HEM.LA0410$Plant <- "LA0410"
HEM.LA0410 <- head(arrange(HEM.LA0410,desc(Effect)), n = 50)

#LA0480: 0.999
HEM.LA0480 <- subset(HEM.plotdata, LA480 > get(paste("TH999_", "LA480", sep="")), 
                     select=c(Chrom,Segment,Pos,LA480))
HEM.LA0480 <- rename(HEM.LA0480, c("LA480" = "Effect"))
HEM.LA0480$Plant <- "LA0480"
HEM.LA0480 <- head(arrange(HEM.LA0480,desc(Effect)), n = 50)

#LA1547 has none > 0.95

#LA1589: 0.999
HEM.LA1589 <- subset(HEM.plotdata, LA1589 > get(paste("TH999_", "LA1589", sep="")), 
                     select=c(Chrom,Segment,Pos,LA1589))
HEM.LA1589 <- rename(HEM.LA1589, c("LA1589" = "Effect"))
HEM.LA1589$Plant <- "LA1589"
HEM.LA1589 <- head(arrange(HEM.LA1589,desc(Effect)), n = 50)

#LA1684: 0.999
HEM.LA1684 <- subset(HEM.plotdata, LA1684 > get(paste("TH999_", "LA1684", sep="")), 
                     select=c(Chrom,Segment,Pos,LA1684))
HEM.LA1684 <- rename(HEM.LA1684, c("LA1684" = "Effect"))
HEM.LA1684$Plant <- "LA1684"
HEM.LA1684 <- head(arrange(HEM.LA1684,desc(Effect)), n = 50)

#LA2093: 0.999
HEM.LA2093 <- subset(HEM.plotdata, LA2093 > get(paste("TH999_", "LA2093", sep="")), 
                     select=c(Chrom,Segment,Pos,LA2093))
HEM.LA2093 <- rename(HEM.LA2093, c("LA2093" = "Effect"))
HEM.LA2093$Plant <- "LA2093"
HEM.LA2093 <- head(arrange(HEM.LA2093,desc(Effect)), n = 50)

#LA2176: 0.999
HEM.LA2176 <- subset(HEM.plotdata, LA2176 > get(paste("TH999_", "LA2176", sep="")), 
                     select=c(Chrom,Segment,Pos,LA2176))
HEM.LA2176 <- rename(HEM.LA2176, c("LA2176" = "Effect"))
HEM.LA2176$Plant <- "LA2176"
HEM.LA2176 <- head(arrange(HEM.LA2176,desc(Effect)), n = 50)

#LA2706: 0.999
HEM.LA2706 <- subset(HEM.plotdata, LA2706 > get(paste("TH999_", "LA2706", sep="")), 
                     select=c(Chrom,Segment,Pos,LA2706))
HEM.LA2706 <- rename(HEM.LA2706, c("LA2706" = "Effect"))
HEM.LA2706$Plant <- "LA2706"
HEM.LA2706 <- head(arrange(HEM.LA2706,desc(Effect)), n = 50)

#LA3008: 0.99 !! or 0.95 to actually get 50 (only 15 > 0.99)
HEM.LA3008 <- subset(HEM.plotdata, LA3008 > get(paste("TH95_", "LA3008", sep="")), 
                     select=c(Chrom,Segment,Pos,LA3008))
HEM.LA3008 <- rename(HEM.LA3008, c("LA3008" = "Effect"))
HEM.LA3008$Plant <- "LA3008"
HEM.LA3008 <- head(arrange(HEM.LA3008,desc(Effect)), n = 50)

#LA3475: 0.95
HEM.LA3475 <- subset(HEM.plotdata, LA3475 > get(paste("TH95_", "LA3475", sep="")), 
                     select=c(Chrom,Segment,Pos,LA3475))
HEM.LA3475 <- rename(HEM.LA3475, c("LA3475" = "Effect"))
HEM.LA3475$Plant <- "LA3475"
HEM.LA3475 <- head(arrange(HEM.LA3475,desc(Effect)), n = 50)

#LA4345: 0.999
HEM.LA4345 <- subset(HEM.plotdata, LA4345 > get(paste("TH999_", "LA4345", sep="")), 
                     select=c(Chrom,Segment,Pos,LA4345))
HEM.LA4345 <- rename(HEM.LA4345, c("LA4345" = "Effect"))
HEM.LA4345$Plant <- "LA4345"
HEM.LA4345 <- head(arrange(HEM.LA4345,desc(Effect)), n = 50)

#LA4355: 0.99 !! or 0.95 to actually get 50. Only 1 SNP > 0.99
HEM.LA4355 <- subset(HEM.plotdata, LA4355 > get(paste("TH95_", "LA4355", sep="")), 
                     select=c(Chrom,Segment,Pos,LA4355))
HEM.LA4355 <- rename(HEM.LA4355, c("LA4355" = "Effect"))
HEM.LA4355$Plant <- "LA4355"
HEM.LA4355 <- head(arrange(HEM.LA4355,desc(Effect)), n = 50)

Top50SNP <- rbind(HEM.LA0410, HEM.LA0480, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)

#max pos is 1001108

#Top50SNP$Plot <- (Top50SNP$Chrom*1000000 + Top50SNP$Pos)

Top50SNP$Chrom <- gsub("Chromosome", "", Top50SNP$Chrom)
Top50SNP$Chrom <- as.numeric(as.character(Top50SNP$Chrom))
Top50SNP$Segment <- as.numeric(as.character(Top50SNP$Segment))
Top50SNP$Pos <- as.numeric(as.character(Top50SNP$Pos))

#sort dataframe rows in order of Chrom, then Pos
Top50SNP <- Top50SNP[with(Top50SNP, order(Chrom, Segment, Pos)), ]

#Make plotting variables
Top50SNP$Index = NA
ticks = NULL
lastbase = 0
Top50SNPhold <- Top50SNP
#subsample DF to work on script
Top50SNPpr <- Top50SNP[sample(nrow(Top50SNP), 20),]
#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
#go through each Chromosome sequentially (i from 1 to 16)
for (i in unique(Top50SNPpr$Chrom)) {
  print(i)
#special rules for Chromosome 1 because counting from zero
  if (i==1) {
    #takes the dataframe (but only the rows that have Chrom = i) and sets Index for those rows to equal...
    #The Pos value for those same rows.
    Top50SNPpr[Top50SNPpr$Chrom==i, ]$Index=Top50SNPpr[Top50SNPpr$Chrom==i, ]$Pos
  }	else {
    print(i,"done")
    }
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(Top50SNP,Top50SNP$Chrom==i-1)$Pos, 1)
    Top50SNP[Top50SNP$Chrom==i, ]$Index=Top50SNP[Top50SNP$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, Top50SNP[Top50SNP$Chrom==i, ]$Index[floor(length(Top50SNP[Top50SNP$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(Top50SNP$Index),max(Top50SNP$Index))

library(ggplot2)
plot1 <- ggplot(Top50SNP, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Plant)))+
  theme_bw()
 
  write.csv(Top50SNP, "data/GWAS_files/04_bigRRoutput/NewModel0711/TopSNPs_alltraits.csv")
  