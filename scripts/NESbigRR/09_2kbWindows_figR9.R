#Nicole E Soltis
#040317
#plot of SNPs/ kb along genome

#---------------------------------------------
rm(list=ls())
library(plyr); library(tidyr); library(ggplot2); library(reshape)
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
#read in MAF20 10NA SNPs file
#can also read in the MAF20 10NA bigRR output file. Retains all SNPs and adds effect estimates
#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")
HEM.plotdata <- HEM.plotdata[,c(2:14,16,17,15)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.Thresh.final.csv")
HEM.thresh <- HEM.thresh[,-1]

#goals:
#plot line 1: % of SNPs significant (any phenotype) within 2kb chunks
#plot line 2: # of traits significant (0-12) within 2kb chunks

#for both: convert phenotypes to 0,1 based on effect > 99% Thr

#positives
Threshval <- HEM.thresh[7,]
names(Threshval)
SNPbin <- HEM.plotdata
for (i in 4:15){
  assign(paste("Thresh_",colnames(SNPbin[i]),sep=''), Threshval[,i-2])
}
#4 to 13
for (i in 4:15){
  fxcol = SNPbin[,paste(colnames(SNPbin[i]),sep='')]
  mythresh = get(paste("Thresh_",colnames(SNPbin[i]),sep=''))
  SNPbin[,paste(colnames(SNPbin[i]),"pos",sep='')] <- ifelse(fxcol > mythresh, 1, 0)
}

#negatives
#but if the positive > 0, force negative to zero.
Threshval <- HEM.thresh[3,]
for (i in 4:15){
  assign(paste("Thresh_",colnames(SNPbin[i]),sep=''), Threshval[,i-2])
}
#4 to 13
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
#for plot 2: Add a column that is sum of phenotype columns
SNPbin$SUMM <- SNPbin$SUMMneg + SNPbin$SUMMpos
#for plot 1: Add a column that is 0 for never significant, 1 for ever significant
SNPbin$Sig <- ifelse(SNPbin$SUMM == 0, 0, 1)

#only keep SUMM and Sig variables for now. Plus Chrom, Seg, Pos, Index.
mydat <- SNPbin[,c(1,2,3,16,43,44)]

#new dataframe: 
#for SNPs/2kb: break index into 2 kb increments, 1 observation per row = number of rows in old df per increment
#Index ranges from 1615 to 41458474 (~40,000 kb)

#build it up
mystart <- min(mydat$Index)
mystop <- max(mydat$Index)
#wanted to try 2kb windows. going to try larger windows
mylist <- seq(mystart, mystop, 20000)
mydfsplit <- split(mydat, cut(mydat$Index, mylist, include.lowest=T))
plotdf <- as.data.frame(seq(1,length(mydfsplit),1))

plotmatrix <- matrix(ncol=6, nrow=length(mydfsplit))

for (i in 1:length(mydfsplit)){
  plotmatrix[i,1] <- nrow(mydfsplit[[i]])
  plotmatrix[i,2] <- mean(mydfsplit[[i]]$Index)
  plotmatrix[i,3] <- mean(mydfsplit[[i]]$Chrom)
  plotmatrix[i,4] <- max(mydfsplit[[i]][,5])
  plotmatrix[i,5] <- mean(mydfsplit[[i]][,5])
  plotmatrix[i,6] <- mean(mydfsplit[[i]][,6]*100)
  #and can extract the phenotypes of interest
  #could be 4:15 for effects sizes, now just 43, 44 for summary plot
  # for (y in 5:6){
  #   plotmatrix[i,y] <- mean(mydfsplit[[i]][,y])
  # }
}
plotdf <- data.frame(plotmatrix)
colnames(plotdf)[1] <- "SNPnum"
colnames(plotdf)[2] <- "Index"
colnames(plotdf)[3] <- "Chrom"
colnames(plotdf)[4] <- "MaxSUMM"
for (y in 5:6){
  colnames(plotdf)[y]= paste(names(mydat[y]))
}

#and remove rows with NaN for index
plotdf <- plotdf[!(is.na(plotdf$Index)),]

#now plot it!!
#jpeg(paste("plots/paper/SlBc_trueMAF20_10NA_lowTR_manhattan",names(HEM.plotdata[i]),".jpg", sep=""), width=8, height=5, units='in', res=600)
plotdf$ScaleSig <- plotdf$Sig / 8.3333333

((max(plotdf[ which(plotdf$Chrom=='1'),]$Index) - min(plotdf[ which(plotdf$Chrom=='1'),]$Index))/2+min(plotdf[ which(plotdf$Chrom=='1'),]$Index))

jpeg(paste("paper/plots/ActualPaper/SharedPhenos", ".GWPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
ggplot(plotdf, aes(x=Index))+
  theme_bw()+
  #colScale+
  geom_point(aes(y=MaxSUMM), colour="aquamarine2")+
        #geom_point(aes(color = factor(Chrom)))+
        labs(y="Phenotypes with Significant SNP")+
        #theme(legend.position="none")+
        scale_x_continuous(name="Chromosome", breaks = c(1671217,5242822, 8992216, 11061288, 13572024, 17204187, 20011250, 22370383, 24391089, 26758962, 28562747, 30110453, 31874774, 33982299, 35769894, 38871001), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16")) +
        expand_limits(y=0)
dev.off()

jpeg(paste("paper/plots/ActualPaper/SNPdistribution", ".GWPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
ggplot(plotdf, aes(x=Index))+
  theme_bw()+
  geom_point(aes(y=Sig), colour="darkorchid4")+
  #geom_point(aes(color = factor(Chrom)))+
  labs(y="Percent of SNPs with effect > 99% threshold")+
  #theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1671217,5242822, 8992216, 11061288, 13572024, 17204187, 20011250, 22370383, 24391089, 26758962, 28562747, 30110453, 31874774, 33982299, 35769894, 38871001), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16")) +
  expand_limits(y=0)
dev.off()


myColors <- c("#999999", "#292929","#684800" ,"#CBA22A", "#63B2D3", "#1FA69D", "#57B761", "#DAD94C","#2B869D", "#EE82EE", "#D2652D", "#CC79A7")
names(myColors) <- levels(Top50SNP$Trait)
colScale <- scale_colour_manual(name = "Plant",values = myColors)




#-------------------------------------------------------
#calculate minor allele frequency
names(SolSNPs)
SolSNPs$Count <- rowSums(SolSNPs[4:96])
table(SolSNPs$Count)
hist(SolSNPs$Count)
SolSNPs$Freq <- (SolSNPs$Count)/91
hist(SolSNPs$Freq)
SolSNPs$MAF <- ifelse(SolSNPs$Freq > 0.5, 1-SolSNPs$Freq, SolSNPs$Freq) 
hist(SolSNPs$MAF)