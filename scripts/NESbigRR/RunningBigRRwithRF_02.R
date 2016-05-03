#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/Sl_LesionSize.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/Sl_LesionSize.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]
TH99 <- HEM.thresh[3,]
TH99_LA0410 <- as.numeric(TH99[2])
TH99_LA0480 <- as.numeric(TH99[3])
TH99_LA1547 <- as.numeric(TH99[4])
TH99_LA1589 <- as.numeric(TH99[5])
TH99_LA1684 <- as.numeric(TH99[6])
TH99_LA2093 <- as.numeric(TH99[7])
TH99_LA2176 <- as.numeric(TH99[8])
TH99_LA2706 <- as.numeric(TH99[9])
TH99_LA3008 <- as.numeric(TH99[10])
TH99_LA3475 <- as.numeric(TH99[11])
TH99_LA4345 <- as.numeric(TH99[12])
TH99_LA4355 <- as.numeric(TH99[13])

TH999 <- HEM.thresh[4,]
TH999_LA0410 <- as.numeric(TH999[2])
TH999_LA0480 <- as.numeric(TH999[3])
TH999_LA1547 <- as.numeric(TH999[4])
TH999_LA1589 <- as.numeric(TH999[5])
TH999_LA1684 <- as.numeric(TH999[6])
TH999_LA2093 <- as.numeric(TH999[7])
TH999_LA2176 <- as.numeric(TH999[8])
TH999_LA2706 <- as.numeric(TH999[9])
TH999_LA3008 <- as.numeric(TH999[10])
TH999_LA3475 <- as.numeric(TH999[11])
TH999_LA4345 <- as.numeric(TH999[12])
TH999_LA4355 <- as.numeric(TH999[13])

TH95 <- HEM.thresh[1,]
TH95_LA0410 <- as.numeric(TH95[2])
TH95_LA0480 <- as.numeric(TH95[3])
TH95_LA1547 <- as.numeric(TH95[4])
TH95_LA1589 <- as.numeric(TH95[5])
TH95_LA1684 <- as.numeric(TH95[6])
TH95_LA2093 <- as.numeric(TH95[7])
TH95_LA2176 <- as.numeric(TH95[8])
TH95_LA2706 <- as.numeric(TH95[9])
TH95_LA3008 <- as.numeric(TH95[10])
TH95_LA3475 <- as.numeric(TH95[11])
TH95_LA4345 <- as.numeric(TH95[12])
TH95_LA4355 <- as.numeric(TH95[13])


# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata2 <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Pos)), ]
range(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Pos)
min(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Pos)
max(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Pos)
range(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Pos)
min(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Pos)
max(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Pos)
#range = min to max
#head and tail not equivalent

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

##NES failing here
#Redo the positions to make them sequencial		-- accurate position indexing
##RFF code isn't working... replaced "Pos" with "pos"
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#make plots for each phenotype
jpeg("plots/Sl_LesionSize_LA0410.ManhattanPlot.jpg")
qplot(Index,abs(LA0410), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA0410", colour=factor(Chrom)) +
 geom_hline(yintercept=TH99_LA0410) +
  geom_text(aes(0,TH99_LA0410, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black") +
geom_hline(yintercept=TH95_LA0410, colour = "blue") +
  geom_text(aes(0,TH95_LA0410, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

jpeg("plots/Sl_LesionSize_LA0410.ManhattanPlot_999.jpg")
qplot(Index,abs(LA0410), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA0410", colour=factor(Chrom)) +
  geom_hline(yintercept=TH999_LA0410) +
  geom_text(aes(0,TH999_LA0410, label = ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")
dev.off()

#how many SNPs are above a certain threshhold?
sum(HEM.plotdata$LA0410 >= TH999_LA0410) #spits out number
#491

#see top 20 SNPs for each
HEM.plotdata$LA0410
require(data.table)
d <- data.table(HEM.plotdata, key="LA0410")
dprint <- d[, head(.SD, 20), by=LA0410]

library(plyr)

maxLA0410 <- HEM.plotdata[order(HEM.plotdata$LA0410,decreasing=T)[1:20],c("Chrom","Pos","LA0410")]
maxLA0410 <- rename(maxLA0410, c("LA0410" = "Effects"))
maxLA0410$Geno <- "LA0410"
maxLA0480 <- HEM.plotdata[order(HEM.plotdata$LA0480,decreasing=T)[1:20],c("Chrom","Pos","LA0480")]
maxLA0480 <- rename(maxLA0480, c("LA0480" = "Effects"))
maxLA0480$Geno <- "LA0480"
maxLA1547 <- HEM.plotdata[order(HEM.plotdata$LA1547,decreasing=T)[1:20],c("Chrom","Pos","LA1547")]
maxLA1547 <- rename(maxLA1547, c("LA1547" = "Effects"))
maxLA1547$Geno <- "LA1547"
maxLA1589 <- HEM.plotdata[order(HEM.plotdata$LA1589,decreasing=T)[1:20],c("Chrom","Pos","LA1589")]
maxLA1589 <- rename(maxLA1589, c("LA1589" = "Effects"))
maxLA1589$Geno <- "LA1589"
maxLA1684 <- HEM.plotdata[order(HEM.plotdata$LA1684,decreasing=T)[1:20],c("Chrom","Pos","LA1684")]
maxLA1684 <- rename(maxLA1684, c("LA1684" = "Effects"))
maxLA1684$Geno <- "LA1684"
maxLA2093 <- HEM.plotdata[order(HEM.plotdata$LA2093,decreasing=T)[1:20],c("Chrom","Pos","LA2093")]
maxLA2093 <- rename(maxLA2093, c("LA2093" = "Effects"))
maxLA2093$Geno <- "LA2093"
maxLA2176 <- HEM.plotdata[order(HEM.plotdata$LA2176,decreasing=T)[1:20],c("Chrom","Pos","LA2176")]
maxLA2176 <- rename(maxLA2176, c("LA2176" = "Effects"))
maxLA2176$Geno <- "LA2176"
maxLA2706 <- HEM.plotdata[order(HEM.plotdata$LA2706,decreasing=T)[1:20],c("Chrom","Pos","LA2706")]
maxLA2706 <- rename(maxLA2706, c("LA2706" = "Effects"))
maxLA2706$Geno <- "LA2706"
maxLA3008 <- HEM.plotdata[order(HEM.plotdata$LA3008,decreasing=T)[1:20],c("Chrom","Pos","LA3008")]
maxLA3008 <- rename(maxLA3008, c("LA3008" = "Effects"))
maxLA3008$Geno <- "LA3008"
maxLA3475 <- HEM.plotdata[order(HEM.plotdata$LA3475,decreasing=T)[1:20],c("Chrom","Pos","LA3475")]
maxLA3475<- rename(maxLA3475, c("LA3475" = "Effects"))
maxLA3475$Geno <- "LA3475"
maxLA4345 <- HEM.plotdata[order(HEM.plotdata$LA4345,decreasing=T)[1:20],c("Chrom","Pos","LA4345")]
maxLA4345 <- rename(maxLA4345, c("LA4345" = "Effects"))
maxLA4345$Geno <- "LA4345"
maxLA4355 <- HEM.plotdata[order(HEM.plotdata$LA4355,decreasing=T)[1:20],c("Chrom","Pos","LA4355")]
maxLA4355 <- rename(maxLA4355, c("LA4355" = "Effects"))
maxLA4355$Geno <- "LA4355"
topSNPs <- rbind(maxLA0410,maxLA0480,maxLA1547,maxLA1589,maxLA1684,maxLA2093,maxLA2176,maxLA2706,maxLA3008,maxLA3475,maxLA4345,maxLA4355)
topSNPs$Num <- rep(c(1:20),12)
topSNPs$Chrom <- as.numeric(topSNPs$Chrom)
topSNPs$Pos <- as.numeric(topSNPs$Pos)
topSNPs$Chr.Pos <- as.numeric(paste(topSNPs$Chrom, topSNPs$Pos, sep="."))

table(unique(topSNPs$Chr.Pos))

hist(topSNPs$Chr.Pos)

library(reshape2)
attach(topSNPs)
wideSNPs <- (dcast(topSNPs, Num~Geno, value.var=c('Chr.Pos'), fun=mean))
pairs(wideSNPs)

#rest of plots
jpeg("plots/Sl_LesionSize_LA1547.ManhattanPlot_95.jpg")
qplot(Index,abs(LA1547), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA1547", colour=factor(Chrom)) +
geom_hline(yintercept=TH95_LA1547, colour = "blue") +
  geom_text(aes(0,TH95_LA1547, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

#view list of Pos indices per Chrom
m <- ggplot(HEM.plotdata, aes(x=Index, fill=as.factor(Chrom)))
m + geom_bar() 

###############################
###Make a BiPlot for two phenotypes (LA0410 and LA0480)
#This code has several dependent objects from the single HEM plotting code above

#set my colors
mycols=rep(c("gray10","gray60"),max(HEM.plotdata$Chrom))

for(k in c(3:7)) {
  #Make temp data for figure
  print(names(thresh.HEM[[1]])[k-2])
  
  CamLes.plot.data<-cbind("Cam",HEM.plotdata[,k],HEM.plotdata[,c(1,13)])
  colnames(CamLes.plot.data)[1:2] <- c("Pheno","Effect")
  CamLes.plot.data <- rbind(CamLes.plot.data,cbind(Pheno="Lesion",Effect=-HEM.plotdata[,k+5],HEM.plotdata[,c(1,13)]))
  CamLes.plot.data <- cbind(CamLes.plot.data,PhenoChrom=paste(CamLes.plot.data[,1],CamLes.plot.data[,3]))
  
  #Scale the effects and threshold
  temp.scaled <- scale(c(thresh.HEM[[3]][k-2],CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Cam"]))
  Cam.Scaled.Thresh <- abs(temp.scaled[1])
  CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Cam"] <- abs(temp.scaled[-1])
  
  temp.scaled <- scale(c(thresh.HEM[[3]][k+3],CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Lesion"]))
  Les.Scaled.Thresh <- -abs(temp.scaled[1])
  CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Lesion"] <- -abs(temp.scaled[-1])
  
  rm(temp.scaled)
  
  
  #Plot Data
  
  mycols2=c("royalblue","royalblue4","royalblue","royalblue4","royalblue","olivedrab3","olivedrab","olivedrab3","olivedrab","olivedrab3")
  
  
  if(grepl("Ctrl",colnames(HEM.plotdata)[k])){
    plot2=qplot(pos,Effect,data=CamLes.plot.data[CamLes.plot.data$Pheno=="Cam",], ylab="|Scaled SNP Effect Estimate|" , colour=factor(PhenoChrom))
    
    ##All the thresholds
    #thresholds <- c(thresh.HEM[[1]][k-2],thresh.HEM[[2]][k-2],thresh.HEM[[3]][k-2],thresh.HEM[[4]][k-2])
    #thresholds <- c(thresh.HEM[[3]][k-2])
    
    plot2=plot2 + geom_hline(yintercept=Cam.Scaled.Thresh, colour="red")
    
  } else {
    plot2=qplot(pos,Effect,data=CamLes.plot.data, ylab="|Scaled SNP Effect Estimate|" , colour=factor(PhenoChrom))
    ##All the thresholds
    #thresholds <- c(thresh.HEM[[1]][c(k-2,k+3)],thresh.HEM[[2]][c(k-2,k+3)],thresh.HEM[[3]][c(k-2,k+3)],thresh.HEM[[4]][c(k-2,k+3)])
    #thresholds <- c(thresh.HEM[[3]][c(k-2,k+3)])
    #thresholds[grepl("Lesion",names(thresholds))] <- thresholds[grepl("Lesion",names(thresholds))] * -1
    plot2=plot2 + geom_hline(yintercept=c(Cam.Scaled.Thresh,Les.Scaled.Thresh), colour="red")
    
    
  }
  #plot2=qplot(pos,Effect,data=CamLes.plot.data, ylab="SNP Effect Estimate" , colour=factor(PhenoChrom))
  plot2=plot2+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(HEM.plotdata$Chrom)))
  #plot2=plot2+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
  plot2=plot2+scale_colour_manual(values=mycols2) + theme(legend.position="none", axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"))
  #plot2=plot2 + geom_hline(yintercept=c(as.numeric(thresh.HEM[k-2]),-as.numeric(thresh.HEM[k+3])), colour="red")
  
  #write to file.  May want something in a vector graphic format
  jpeg(filename = paste(strsplit(names(thresh.HEM[[1]][k-2]),"_")[[1]][2],".jpeg",sep=""), width = 1500, height = 800)
  print(plot2)
  dev.off()
}



###############################
###Make a HEM plot of a region (Cam and Lesion)
#This code has several dependent objects from the single HEM plotting code above

#PRLM3 (Chr4 Pos 9560135:9565454) in this example.  Just change the Chrom and positions to find another gene.
Reg.min <- 9560135 - 50000
Reg.max <- Reg.min + 1000000


Reg.min <- 13750000 - 50000
Reg.max <- Reg.min + 100000


Region2Plot <- HEM.plotdata[HEM.plotdata$Chrom == 4 & HEM.plotdata$Pos %in% c(Reg.min:Reg.max),]
Region2Plot <- Region2Plot[,-match("pos",colnames(Region2Plot))]
Region2Plot <- melt(Region2Plot, id=c("Chrom","Pos"), variable.name = "Isolate", value.name = "Effect")
Region2Plot$Pheno <- NA
Region2Plot$Pheno[grepl("Cam",Region2Plot$Isolate)] <- "Camalexin"
Region2Plot$Pheno[grepl("Lesion",Region2Plot$Isolate)] <- "Lesion"


Region2Plot$Isolate <- gsub("Camalexin.ng.LesPer._|Lesion.Size.mm2_|.HEM","",Region2Plot$Isolate)

plot.reg <- qplot(Pos,abs(as.numeric(Effect)),data=Region2Plot[Region2Plot$Pheno == "Camalexin",], ylab="|Scaled SNP Effect Estimate|" , colour=factor(Isolate))

#Missing a gene at 13755649-13759740; At4g27550 (wouldn't print)
gene.start = c(13759841,13763464,13765661,13767763,13771185,13772819)
gene.stop = c(13761559,13765064,13766597,13769961,13772594,13777290)

plot.reg <- plot.reg + geom_segment(aes(x = gene.start, y = -0.001, xend = gene.stop, yend = -.001), colour = "red")

#####
###For some reason, I had to use the following code at the terminal to get this to print right.
#xwd > PolyMorphicRegOnChr4.xwd
#convert PolyMorphicRegOnChr4.xwd PolyMorphicRegOnChr4.jpg

jpeg(filename = "PolyMorphicRegionOnChr4.jpg", width = 1500, height = 800)
print(plot2)
dev.off()