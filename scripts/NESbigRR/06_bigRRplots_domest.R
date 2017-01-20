#Nicole E Soltis
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
#setwd("~/Projects/BcSolGWAS/")

#Input file: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv AND Sl_DomesticationLS_MAF20.HEM.Thresh.csv
#Output file: NONE
#Plots: Basic and greyscale domestication Manhattan plots
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
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata2 <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Pos)), ]

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=+lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#make plots for each phenotype
#it isn't working with the loop for some reason
3:ncol(HEM.plotdata)
for (y in names(HEM.plotdata[,4:6])){
print(y)
}

#create a custom color scale
library(RColorBrewer)
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#when plotting, add +colScale + as a line

#without the loop [6]
for (y in 4:6){
jpeg(paste("plots/MultiPlot/domest/Sl_LesionSize_MAF20_lowTR_bw_", names(HEM.plotdata[6]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
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
   expand_limits(y=0)
dev.off()

jpeg(paste("plots/MultiPlot/domest/Sl_LesionSize_MAF20_highTR_", names(HEM.plotdata[6]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index, y=abs(HEM.plotdata[6])))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[6]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH999_", names(HEM.plotdata[6]), sep=""))) +
  geom_text(aes(0,get(paste("TH999_", names(HEM.plotdata[6]), sep="")), label =
  ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")+
  expand_limits(y=-0.001)
dev.off()
#}
#stop [6]

#how many SNPs are above a certain threshhold?
topSNP <- sum(HEM.plotdata[6] >= 0.001) #41
highSNP <- sum(HEM.plotdata[6] >= get(paste("TH999_", names(HEM.plotdata[6]), sep=""))) #spits out number
totalSNP <- sum(HEM.plotdata[6] >= 0)
highSNP/totalSNP*100