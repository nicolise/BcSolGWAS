#Nicole E Soltis
#081117
#plot of SNPs along gene of interest

#12_singleGeneManhattan.R
#---------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")

#take the SNPs over the threshold for each phenotype

TH99pos <- HEM.thresh[7,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[3,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}


names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,-c(1)]
#only look at chromosome 16
HEM.plotdata <- HEM.plotdata[which(HEM.plotdata$Chrom=='16'),]


#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in c(4:15)){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}

#combine pos and neg by group
for (i in c(4:15)){
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
}

#also keep non-sig SNPs as grey dots
for (i in c(4:15)){
  assign(paste("ns.HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
  assign(paste("ns.HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}
for (i in c(4:15)){
  assign(paste("ns.HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("ns.HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("ns.HEMneg.", names(HEM.plotdata[i]), sep=""))))
}
for (i in c(4:15)){
  mydf <- paste("ns.HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[5] <- "Effect"
  assign(mydf, renamedf)
  myblob <- rep("NS", nrow(get(mydf)))
  assign(mydf, cbind(get(mydf), Trait = myblob))
}
HEM.NS <- rbind(ns.HEM.LA410, ns.HEM.LA480, ns.HEM.LA1547, ns.HEM.LA1589, ns.HEM.LA1684, ns.HEM.LA2093, ns.HEM.LA2176, ns.HEM.LA2706, ns.HEM.LA3008, ns.HEM.LA3475, ns.HEM.LA4345, ns.HEM.LA4355)

#then combine
#4:15
for (i in c(4:15)){
  mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[5] <- "Effect"
  assign(mydf, renamedf)
  myblob <- rep(names(HEM.plotdata[i]), nrow(get(mydf)))
  assign(mydf, cbind(get(mydf), Trait = myblob))
}

HEM.topSNPs <- rbind(HEM.LA410, HEM.LA480, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)

HEM.topSNPs <- rbind(HEM.NS, HEM.topSNPs)
#create a custom color scale
#myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
myColors <- c("gray97", "#999999", "#292929","#684800" ,"#CBA22A", "#63B2D3", "#1FA69D", "#57B761", "#DAD94C","#2B869D", "#EE82EE", "#D2652D", "#CC79A7")
names(myColors) <- levels(HEM.topSNPs$Trait)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#get the start position of chromosome 16
min(HEM.topSNPs$Index)
max(HEM.topSNPs$Index)
max(HEM.topSNPs$Index) - min(HEM.topSNPs$Index)

#narrow window: +- 1 kb
HEM.topSNPs$Chr16Index <- HEM.topSNPs$Index - min(HEM.topSNPs$Index) + 1
min(HEM.topSNPs$Chr16Index)
#now get target region within chromosome 16
#my feature: about 1kb
#345785 to 346542
#and I'll add 2kb on each side
HEM.topSNPsSM <- HEM.topSNPs[which(HEM.topSNPs$Chr16Index < 347542),]
HEM.topSNPsSM <- HEM.topSNPsSM[which(HEM.topSNPsSM$Chr16Index > 344785),]

HEM.topSNPs.P <- HEM.topSNPsSM
jpeg("paper/plots/FigR8/Sl_LesionSize_trueMAF20_NA10_lowTR.gene01Chr16.ManhattanPlot.jpg", width=8, height=5, units='in', res=600)
  ggplot(HEM.topSNPs.P, aes(x=Chr16Index, y=Effect))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Trait)))+
         labs(list(y="SNP Effect Estimate"))+
         guides(col = guide_legend(nrow = 8, title="SNP position"))+
         theme(legend.position="none")+
         scale_x_continuous(name="SNP position")+
    #exon 1 345785	345965
         geom_rect(mapping=aes(ymin=-1.6e-05, ymax=-1.4e-05, xmin=345785, xmax=345965))+
    #exon 2 346028	346542
         geom_rect(mapping=aes(ymin=-1.6e-05, ymax=-1.4e-05, xmin=346028, xmax=346542))+
    #breaks = c(1), labels = c("1")
         expand_limits(y=0)
dev.off()
