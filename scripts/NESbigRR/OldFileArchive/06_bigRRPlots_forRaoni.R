#Nicole E Soltis
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcSolGWAS/")

#Plots: Basic and greyscale Manhattan plots
###########################################################################
#Plotting the HEM results

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

#make plots for each phenotype

#create a custom greyscale color scale for the chromosomes
library(RColorBrewer)
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)
#when plotting, add +colScale + as a line if you want it in greyscale

#without the loop 
#find [6] and replace each with the number for the phenotype you want to plot
#95% and 99% version of the plot
  jpeg(paste("plots/Gm_LesionSize_MAF20_lowTR_bw_", names(HEM.plotdata[6]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
  ggplot(HEM.plotdata, aes(x=Index, y=abs(HEM.plotdata[6])))+
    theme_bw()+
    geom_point(aes(color = factor(Chrom)))+
    colScale+
    labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[6]))))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    geom_hline(yintercept=get(paste("TH95_", names(HEM.plotdata[6]), sep="")), colour = "blue") +
    geom_text(aes(0,get(paste("TH95_", names(HEM.plotdata[6]), sep="")), label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")+
    geom_hline(yintercept=get(paste("TH99_", names(HEM.plotdata[6]), sep=""))) +
    geom_text(aes(0,get(paste("TH99_", names(HEM.plotdata[6]), sep="")), label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black")+
    scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    expand_limits(y=0)
  dev.off()
  
#this is the 99.9% threshold plot
  jpeg(paste("plots/MultiPlot/domest/Sl_LesionSize_MAF20_highTR_", names(HEM.plotdata[6]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
  ggplot(HEM.plotdata, aes(x=Index, y=abs(HEM.plotdata[6])))+
    theme_bw()+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[6]))))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    geom_hline(yintercept=get(paste("TH999_", names(HEM.plotdata[6]), sep=""))) +
    geom_text(aes(0,get(paste("TH999_", names(HEM.plotdata[6]), sep="")), label =  ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")+
    scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    expand_limits(y=-0.001)
  dev.off()