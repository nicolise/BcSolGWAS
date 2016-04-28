#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files/bigRRinput")
setwd("~/Projects/BcSolGWAS/data/bigRRout")
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data from file modified in script "ReformatBigRRouts.R"
HEM.plotdata <- read.csv("Sl_LesionSizePoisson.HEM.PlotFormat.csv")
HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation

#RF-Extract thresholds - 4 rows of 4 threshhold values, we just want the result of nonpermuted data to plot
HEM.Perm.Threshold <- as.vector(HEM.plotdata[1,])
HEM.plotdata <- HEM.plotdata[-1:4,] #already did this 

#get threshhold values from HEM file
HEM.thresh <- read.csv("Sl_LesionSizePoisson.HEM.PlotFormat.csv")
TH99 <- HEM.thresh[3,]
TH99_Apple517 <- as.numeric(TH99[3])
TH99_B05.10 <- as.numeric(TH99[4])
TH99_Supersteak <- as.numeric(TH99[5])
TH99_UKRazz <- as.numeric(TH99[6])

TH999 <- HEM.thresh[4,]
TH999_Apple517 <- as.numeric(TH999[3])
TH999_B05.10 <- as.numeric(TH999[4])
TH999_Supersteak <- as.numeric(TH999[5])
TH999_UKRazz <- as.numeric(TH999[6])

TH95 <- HEM.thresh[1,]
TH95_Apple517 <- as.numeric(TH95[3])
TH95_B05.10 <- as.numeric(TH95[4])
TH95_Supersteak <- as.numeric(TH95[5])
TH95_UKRazz <- as.numeric(TH95[6])



#Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("IV", 4, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("V", 5, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("III", 3, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("II", 2, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("I", 1, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))


#Give column names to the thresholds from the HEM list
#for(j in 1:length(thresh.HEM)) {
#  names(thresh.HEM[[j]]) <- colnames(dat)
#}
#RF - we already ran this command above

#Make plotting variables
HEM.plotdata$Index=NA
ticks = NULL
lastbase=0

#Redo the positions to make them sequencial		
##RFF code isn't working... replaced "Pos" with "pos"
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    lastbase=lastbase+tail(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
   # previous = lastbase+tail(subset(HEM.plotdata,Chrom==i-1)$Pos, 1) 
    # lastbase=lastbase+tail(subset(HEM.plotdata,Chrom==previous)$Pos, 1)
    #lastbase=lastbase+tail(subset(HEM.plotdata,Chrom==i-1)$Pos, 1)
   # HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos = HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#make plots for each phenotype
jpeg( "LesionEccentricityPoisson_Apple517.95.ManhattanPlot.jpg")

qplot(Index,abs(Lesion.0.m.eccentricity_Apple517), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionEccentricityPoisson_Apple517", colour=factor(Chrom)) +
 geom_hline(yintercept=TH99_Apple517) +
  geom_text(aes(0,TH99_Apple517, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black") +
geom_hline(yintercept=TH95_Apple517, colour = "blue") +
  geom_text(aes(0,TH95_Apple517, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")

dev.off()
#how many SNPs are above a certain threshhold?
sum(HEM.plotdata$Lesion.0.m.eccentricity_Apple517 >= TH99_Apple517) #spits out number
#Eccentricity = 597, Area = 1709, Grn = 1573

#try ggplot? - must add "environment" argument - code below doesn't work yet
# ggplot(data=HEM.plotdata, aes(HEM.plotdata$Index, abs(Lesion.0.m.eccentricity_Apple517)), ylab="SNP Effect Estimate" , 
#        main = "Lesion.0.m.eccentricity_Apple517", colour=factor(Chrom)) +
#  geom_hline(aes(yintercept=TH99))


 
sum(HEM.plotdata$Lesion.0.m.eccentricity_B05.10 >= TH99_B05.10) #spits out number



 


jpeg(file = "LesionEccentricityPoisson_UKRazz.95.ManhattanPlot.jpg")
qplot(Index,abs(Lesion.0.m.eccentricity_UKRazz), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionEccentricityPoisson_UKRazz", colour=factor(Chrom)) +
  geom_hline(yintercept=TH99_UKRazz) +
  geom_text(aes(0,TH99_UKRazz, label = ".99 Threshhold", vjust = 1.5, hjust = .05), col = "black") +
geom_hline(yintercept=TH95_UKRazz, colour = "blue") +
  geom_text(aes(0, TH95_UKRazz, label = ".95 Threshhold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

#make plot with all phenos?



#set my colors
mycols=rep(c("gray10","gray60"),max(HEM.plotdata$Chrom))

#Make plot (single)
# plot=qplot(Index,abs(Lesion.Size_Apple517.HEM),data=HEM.plotdata, ylab="SNP Effect Est" , colour=factor(Chrom))
# plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(HEM.plotdata$Chrom)))
# plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
# plot=plot+scale_colour_manual(values=mycols)
# plot=plot + geom_hline(yintercept=HEM.Perm.Threshold$Lesion.Size_Apple517.HEM, colour="red")


###############################
###Make a BiPlot for two phenotypes (Cam and Lesion)
#This code has several dependent objects from the single HEM plotting code above


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








