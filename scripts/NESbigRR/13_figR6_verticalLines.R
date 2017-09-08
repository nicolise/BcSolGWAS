#Nicole E Soltis 
#082317
#short script to find overlap between 2 GWAS visualizations
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#just Manhattan plots for individual plant genotypes
############################################################################
###Plotting the HEM results

#NEED TO CHECK CONTIGS FOR THIS

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata.og <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")
SNPlist  <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")
Thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.Thresh.final.csv")

names(SNPlist)
SNPlist <- SNPlist[,c(2:14,16,17,15)]
names(Thresh)
Thresh <- Thresh[,c(2:14)]

#positive: over 99% positive thresh
Threshval <- Thresh[7,]
names(Threshval)
SNPbin <- SNPlist
for (i in 2:13){
  assign(paste("Thresh_",colnames(SNPbin[i+2]),sep=''), Threshval[,i])
}
#4 to 15
for (i in 4:15){
  fxcol = SNPbin[,paste(colnames(SNPbin[i]),sep='')]
  mythresh = get(paste("Thresh_",colnames(SNPbin[i]),sep=''))
  SNPbin[,paste(colnames(SNPbin[i]),"pos",sep='')] <- ifelse(fxcol > mythresh, 1, 0)
}

#negatives
#but if the positive > 0, force negative to zero.
Threshval <- Thresh[3,]
for (i in 2:13){
  assign(paste("Thresh_",colnames(SNPbin[i+2]),sep=''), Threshval[,i])
}
#4 to 15
for (i in 4:15){
  fxcol = SNPbin[,paste(colnames(SNPbin[i]),sep='')]
  mythresh = get(paste("Thresh_",colnames(SNPbin[i]),sep=''))
  poscol = SNPbin[,paste(colnames(SNPbin[i]),"pos",sep='')]
  SNPbin[,paste(colnames(SNPbin[i]),"neg",sep='')] <- ifelse( poscol > 0, 0, ifelse(fxcol < mythresh, 1, 0))
}

#now create a summation column
names(SNPbin)
SNPbin$SUMMpos <- rowSums(SNPbin[,c(17:28)])
SNPbin$SUMMneg <- rowSums(SNPbin[,c(29:40)])
SNPbin$SUMMneg
SNPbin$SUMM <- SNPbin$SUMMneg + SNPbin$SUMMpos

table(SNPbin$SUMM)

#high overlap SNP list for annotation
HOSNP <- SNPbin[SNPbin$SUMM > 6,]

HOSNP$v <- paste(HOSNP$Chrom,HOSNP$Segment,HOSNP$Pos, sep=".")
HOSNP.v <- HOSNP$v
HOSNP.v <- HOSNP$Index

#and now for single plot
HEM.plotdata <- HEM.plotdata.og
HEM.plotdata <- HEM.plotdata[,-c(1)]

blah <- HEM.plotdata[order(abs(HEM.plotdata$LA2093)),]
blah <- tail(blah,100)
blah.v <- paste(blah$Chrom, blah$Segment, blah$Pos, sep=".")
blah.v <- blah$Index

intersect(HOSNP.v, blah.v)

#ones I kept:

+
  geom_vline(xintercept=1099438, lty=2)+
  geom_vline(xintercept=1836245, lty=2)+
  geom_vline(xintercept=31154483, lty=2)+
  geom_vline(xintercept=33853054, lty=2)+
  geom_vline(xintercept=36555407, lty=2)+
  geom_vline(xintercept=41350409, lty=2)+
  #geom_vline(xintercept=5703231, lty=2)+ #no apparent SNP > 6 phenos
  geom_vline(xintercept=6589294, lty=2)+
  geom_vline(xintercept=7955289, lty=2)+
  geom_vline(xintercept=11188054, lty=2)+
  geom_vline(xintercept=22332692, lty=2)+
  geom_vline(xintercept=24530790, lty=2)
)

#full list:
#line 1: 1099438  1099555  1099716  1099797  1100254  1100490  1100913  1100964  1100985  1101262  1101421 1101487  1105689  
#line 2: 1836245 
#line 3: 31154483 
#line 4: 33853054 33853095 33853132 33853211 33853358 33853411 
#line 5: 36555407
#line 6: 41350409  
#line 7: 5703231  5703285  5706355  5706501  5707812  5708041  5708421  5708914  5709576  5709578 5710912  5711683  5713025  5713606  5713873  5714065  5714544  5714555  5714601  5714751  5714910  5761086  5775334  
#line 8: 6589294  6603140  6603417  6605325  6607536  6608006  6612290  6614855  6615095 6615710  6615923  6616279  6619063  6619749  6619976  6624262  6624491  6628516  6628531  6628543 6628558  6628623  
#line 9: 7955289 
#line 10: 11188054 
#line 11: 22332692 
#line 12: 24530790 24530862 24531371 24531377 24531505 24531506 24531527 24532002 24532500 24533334 24534056 24539133 24539645 24539972 24540981 24542254 24545685



