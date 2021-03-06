#Nicole E Soltis 
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input File: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv and .Thresh.csv
#Output File: results/Domestication_TopSNPs_SegLong.csv, results/Domestication_TopSNPs_SegWide.csv
#Plots: NONE
############################################################################

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

#Reformat Chromosomes and Positions
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

names(HEM.plotdata)
#All groups (4:6)
#keep only: SNPs over 99.9% Threshold
#assign(paste("HEM.", names(HEM.plotdata[6]), sep=""), subset(HEM.plotdata, abs(HEM.plotdata[6]) > get(paste("TH999_", names(HEM.plotdata[6]), sep="")), select=c(Chrom,Segment,Pos,Index,6)))

#keeping all the SNPs
assign(paste("HEM.", names(HEM.plotdata[4]), sep=""), subset(HEM.plotdata, abs(HEM.plotdata[4]) > 0, select=c(Chrom,Segment,Pos,Index,4)))

#then combine
HEM.Domesticated <- rename(HEM.Domesticated, c("Domesticated" ="Effect"))
HEM.Domesticated$Trait <- "Domesticated"
HEM.Wild <- rename(HEM.Wild, c("Wild" ="Effect"))
HEM.Wild$Trait <- "Wild"
HEM.DmWoD <- rename(HEM.DmWoD, c("DmWoD" ="Effect"))
HEM.DmWoD$Trait <- "DmWoD"

HEM.plotdata <- rbind(HEM.Domesticated, HEM.Wild, HEM.DmWoD)

#now narrow it down to plot
HEM.plotdata.6.1 <- HEM.plotdata[which(HEM.plotdata$Chrom=='6' &
                                         HEM.plotdata$Segment=='1'),]
HEM.plotdata.2.3 <- HEM.plotdata[which(HEM.plotdata$Chrom=='2' &
                                         HEM.plotdata$Segment=='3'),]

HEM.plotdata.2 <- HEM.plotdata[which(HEM.plotdata$Index<= 6069302 &
                                     HEM.plotdata$Index>= 6048940
),]

HEM.plotdata.6 <- HEM.plotdata[which(HEM.plotdata$Index<= 17827395 &
                                       HEM.plotdata$Index>= 17812613
),]

#domesticated is blue, wild is orange, sensitivity is black
myColors <- c("#050505", "#1C86EE", "#EE7600")

#names(myColors) <- levels(ModDat$Species)
colScale <- scale_colour_manual(name = "Species",values = myColors)



library(ggplot2)
#plot 6
plot1 <- ggplot(HEM.plotdata.6.1, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Trait)))+
  colScale+
  scale_fill_discrete(labels=c("Domestication Sensitivity", "Domesticated", "Wild"))+
  labs(list(y="SNP Effect Estimate",  title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title=element_blank()))+
  scale_x_continuous(name="Chromosome Position", breaks=c(17812000, 17816000, 17820000, 17824000, 17828000), labels=c("0kb","4kb","8kb", "12kb", "16kb"))+
  theme_bw()

#plot 2
plot1 <- ggplot(HEM.plotdata.2, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Trait)))+
  colScale+
  scale_fill_discrete(labels=c("Domestication Sensitivity", "Domesticated", "Wild"))+
  labs(list(y="SNP Effect Estimate",  title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title=element_blank()))+
  scale_x_continuous(name="Chromosome Position", breaks=c(6045000, 6050000, 6055000, 6060000, 6065000, 6070000), labels=c("0kb","5kb","10kb", "15kb", "20kb", "25kb"))+
  theme_bw()
  
#make it wide format
names(HEM.plotdata)
#currentlylong format : Chrom, Segment, Pos, Index, Effect, Trait
write.csv(HEM.plotdata, "results/Domestication_TopSNPs_SegLong.csv")

Top50SNP.wide.DM <- reshape(HEM.plotdata, 
                         timevar = "Trait",
                         idvar = c("Chrom","Segment","Pos","Index"),
                         direction = "wide")

write.csv(Top50SNP.wide.DM, "results/Domestication_TopSNPs_SegWide.csv")
