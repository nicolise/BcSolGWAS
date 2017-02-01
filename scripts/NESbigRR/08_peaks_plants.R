#Nicole E Soltis
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input file: Sl_LesionSize_MAF20.HEM.PlotFormat.csv AND Sl_LesionSize_MAF20.HEM.Thresh.csv
#Output file: TopSNPs_alltraits.csv
#Plots: NONE

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

# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Segment <- as.numeric(as.character(HEM.plotdata$Segment))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Segment, Pos)), ]

#now make segments line up consecutively
HEM.plotdata$Chrom.Seg <- paste(HEM.plotdata$Chrom, HEM.plotdata$Segment, sep=".")
HEM.plotdata$Chrom.Seg <- as.numeric(HEM.plotdata$Chrom.Seg)

#let's try making the chrom.seg integers so that R isn't confused
unique(HEM.plotdata$Chrom.Seg)
library(plyr)
HEM.plotdata$Chrom.Seg.F <- as.factor(HEM.plotdata$Chrom.Seg)
unique(HEM.plotdata$Chrom.Seg.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6,  6.1, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 13.2, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

HEM.plotdata$Chrom.Seg.Int <- recode.vars$newvals[match(HEM.plotdata$Chrom.Seg.F, recode.vars$OGvals)]
unique(HEM.plotdata$Chrom.Seg.Int)

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing

for (i in unique(HEM.plotdata$Chrom.Seg.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom.Seg.Int==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Pos+lastbase
  }
  #set ticks to be a list of existing ticks, plus the current Index
  #floor rounds it down to the nearest whole number
  # ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
  
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

library(ggplot2)
plot1 <- ggplot(HEM.plotdata, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Plant)))+
  theme_bw()
 
  write.csv(HEM.plotdata, "data/GWAS_files/04_bigRRoutput/NewModel0711/TopSNPs_alltraits_Segments.csv")
  
names(HEM.plotdata)
HEM.plotdata
#now long format : Chrom, Segment, Pos, Index, Effect, Trait
HEM.plotdata
  