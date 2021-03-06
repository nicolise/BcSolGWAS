#Nicole E Soltis
#081117
#plot of SNPs along gene of interest

#12_singleGeneManhattan.R
#---------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF20_20NA/SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF20_20NA/SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")

#take the SNPs over the threshold for each phenotype

TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}


names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,-c(1)]
#only look at chromosome 16
#this includes all contigs (chrom and seg are separated)
HEM.plotdata <- HEM.plotdata[which(HEM.plotdata$Chrom=='16'),]

#get the start position of chromosome 16
min(HEM.plotdata$Index)
max(HEM.plotdata$Index)
max(HEM.plotdata$Index) - min(HEM.plotdata$Index)

#narrow window: +- 1 kb
HEM.plotdata$Chr16Index <- HEM.plotdata$Index - min(HEM.plotdata$Index) + 1
min(HEM.plotdata$Chr16Index)
#now get target region within chromosome 16
#my feature: about 1kb
#345785 to 346542
#and I'll add 2kb on each side
HEM.plotdataSM <- HEM.plotdata[which(HEM.plotdata$Pos < 347542),]
HEM.plotdataSM <- HEM.plotdataSM[which(HEM.plotdataSM$Pos > 344785),]

#trying a bigger window (8kb) to find missing phenos
#HEM.plotdataSM <- HEM.plotdata[which(HEM.plotdata$Pos < 355042),]
#HEM.plotdataSM <- HEM.plotdataSM[which(HEM.plotdataSM$Pos > 337285),]

#and remove any duplicated POS here
SNPlist_8a <- HEM.plotdataSM[!duplicated(HEM.plotdataSM$Pos), ]

#SNPlist <- as.data.frame(SNPlist_8a$Pos)
#write.csv(SNPlist,"data/genome/chr16_analysis/SNPlistFig8a.csv")

HEM.plotdata <- HEM.plotdataSM

#add a new variable so we can plot in order of seg.pos
HEM.plotdata$Seg.Pos <- paste(HEM.plotdata$Segment, HEM.plotdata$Pos, sep='.')
HEM.plotdata <- HEM.plotdata[order(HEM.plotdata$Segment, HEM.plotdata$Pos),]

#Make newindex for chromosome 16
HEM.plotdata$newindex = NA
lastbase = 0
for (i in unique(HEM.plotdata$Segment)) {
  print(i)
  if (i<1) {
    HEM.plotdata[HEM.plotdata$Segment==i, ]$newindex=HEM.plotdata[HEM.plotdata$Segment==i, ]$Pos
  }	else {
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Segment==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Segment==i, ]$newindex=HEM.plotdata[HEM.plotdata$Segment==i, ]$Pos+lastbase-min(subset(HEM.plotdata,HEM.plotdata$Segment==i)$Pos)
  }
}

getmyindex <- HEM.plotdata[,c(3,22)]

#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in c(4:15)){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,newindex,i)))
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,newindex,i)))
}

#combine pos and neg by group
for (i in c(4:15)){
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
}

#also keep non-sig SNPs as grey dots
for (i in c(4:15)){
  assign(paste("ns.HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,newindex,i)))
  assign(paste("ns.HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,newindex,i)))
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
HEM.topSNPsB <- rbind(HEM.LA410, HEM.LA480, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)


HEM.topSNPs <- rbind(HEM.NS, HEM.topSNPsB)

#create a custom color scale
#myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
myColors <- c("gray85", "#999999", "#292929","#684800" ,"#CBA22A", "#63B2D3", "#1FA69D", "#57B761", "#DAD94C","#2B869D", "#EE82EE", "#D2652D", "#CC79A7")
names(myColors) <- levels(HEM.topSNPs$Trait)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#modify colors so that wild are oranges and domesticated are blues
levels(HEM.topSNPs$Trait)
#"NS"     "LA410"  "LA480"  "LA1547" "LA1589" "LA1684" "LA2093"
#"LA2176" "LA2706" "LA3008" "LA3475" "LA4345" "LA4355"
#NS, D, W, W, W, W, W, W, D, D, D, D, D
#oranges: coral1, deeppink3, yellow1, darkorange, red4, pink1
#blues: darkblue, dodgerblue1, blueviolet, lawngreen, seagreen4, mediumorchid1
myColors <- c("gray85", "darkblue", "coral1", "deeppink3", "goldenrod2", "darkorange", "red4", "pink1", "dodgerblue1", "blueviolet", "lawngreen", "seagreen4", "mediumorchid1")

names(myColors) <- levels(HEM.topSNPs$Trait)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

HEM.topSNPs.P <- HEM.topSNPs

#add lines for plot
myblocks <- read.csv("data/genome/chr16_analysis/BlockBoundaries.csv")
block2snp <- read.csv("data/genome/chr16_analysis/plink/fig8aMatch/MatchDrawLines.csv")
names(myblocks)[1] <- "SNPnum"
mylines <- merge(myblocks, block2snp, by="SNPnum")
#add indexing here!
mylines <- mylines[,c(1,3)]
mylines <- merge(mylines, getmyindex, by="Pos")

names(HEM.topSNPs.P)
HEM.topSNPs.P <- HEM.topSNPs.P[HEM.topSNPs.P$Segment=="0",]

#jpeg("paper/plots/FigR8/Sl_LesionSize_trueMAF20_NA10_lowTR.gene01Chr16.ManhattanPlot.jpg", width=7, height=5, units='in', res=600)
jpeg("paper/plots/FigR8/Sl_LesionSize_trueMAF20_NA10_lowTR.gene01Chr16.ManhattanPlot.jpg", width=7, height=5, units='in', res=600)
  ggplot(HEM.topSNPs.P, aes(x=newindex, y=100*Effect))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Trait)))+
         labs(list(y=expression(paste("Estimated SNP Effect Size (",mm^{2},")"))))+
         guides(col = guide_legend(nrow = 8, title="SNP position"))+
         theme(legend.position="none")+
         scale_y_continuous(breaks=c(3e-03, 2e-03, 1e-03, 0, -1e-03, -2e-03, -3e-03))+
         scale_x_continuous(name="SNP position on Chr 16 (kb)", limits=c(344500, 350500), 
                            #this is by Pos
                            #breaks=c(344500, 345000, 345500, 346000, 346500, 347000, 347500), 
                            #this is by newindex
                            breaks=c(345000, 346000, 347000, 348000, 349000, 350000),
                            labels=c("345", "346", "347", "348", "349", "350"))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
    #exon 1 345785	345965 -- confirmed for newindex 345785 348638
    geom_rect(mapping=aes(ymin=-0.2e-03, ymax=0.2e-03, xmin=345785, xmax=348638), alpha=0.01, fill="darkturquoise")+
    #exon 2 346028	346542 -- confirmed for newindex 348701 349215
    geom_rect(mapping=aes(ymin=-0.2e-03, ymax=0.2e-03, xmin=348701, xmax=349215), alpha=0.01, fill="darkturquoise")+
    geom_vline(xintercept=c(347518,344825,347519,347575,347613,347730,347517,347820,347918,348343,348571,348621,348944,349406,349407,349409,349441
),linetype="dotted")+
    expand_limits(y=0)
dev.off()
