#Nicole E Soltis 
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input File: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv and .Thresh.csv
#Output File: results/Domestication_TopSNPs_SegLong.csv, results/Domestication_TopSNPs_SegWide.csv
#this goes into gene annotation and then venn diagrams
#Plots: R6_topSNPs
############################################################################

#Load packages
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.Thresh.final.csv")

#take the top 50 over the threshold for each phenotype
TH95pos <- HEM.thresh[5,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}
TH95neg <- HEM.thresh[1,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}
TH99pos <- HEM.thresh[7,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[3,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}
TH999pos <- HEM.thresh[8,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}
TH999neg <- HEM.thresh[4,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

names(HEM.plotdata)
HEM.plotdata <- HEM.plotdata[,c(2:14,16,17,15)]
#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in c(4:15)){
assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}

#skip to below for top 50
#for top 1000 only
for (i in c(4:15)){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
#negative
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), tail(arrange(get(paste("HEMneg.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
  #combine pos and neg by group
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
#then combine all
mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
renamedf <- get(mydf)
colnames(renamedf)[5] <- "Effect"
assign(mydf, renamedf)
myblob <- rep(names(HEM.plotdata[i]), 1000)
assign(mydf, cbind(get(mydf), Trait = myblob))
}

HEM.topSNPs <- rbind(HEM.LA410, HEM.LA480, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)

library(ggplot2)
plot1 <- ggplot(HEM.topSNPs, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Trait)))+ theme_bw()

#make it wide format
#currently long format : Chrom, Segment, Pos, Index, Effect, Trait
#write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")
write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1000SNPs_SegLong_trueMAF20_20NA.csv")

TopSNP.wide.DM <- reshape(HEM.topSNPs, 
                         timevar = "Trait",
                         idvar = c("Chrom","Segment","Pos","Index"),
                         direction = "wide")

write.csv(TopSNP.wide.DM, "results/12Plants_Top1000SNPs_SegWide_trueMAF20_20NA.csv")

#-------------------------------------------------------------
#plot for R6
#top SNPs overlap (figure R6)
for (i in c(4:15)){
  assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[i]), sep=""))[,5])), n=100))
  #negative
  assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), tail(arrange(get(paste("HEMneg.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))[,5])), n=100))
  #combine pos and neg by group
  assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
  #then combine all
  mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
  renamedf <- get(mydf)
  colnames(renamedf)[5] <- "Effect"
  assign(mydf, renamedf)
  myblob <- rep(names(HEM.plotdata[i]), nrow(renamedf))
  assign(mydf, cbind(get(mydf), Trait = myblob))
}
#now scale effect size for each plant
#as a loop
for (i in c(4:15)){
  dfname <- paste("HEM.", names(HEM.plotdata[i]), sep="")
  mydf <- get(dfname)
  maxfx <- max(abs(mydf$Effect))
  mydf$Effect.s <- mydf$Effect / maxfx
  assign(dfname, mydf)
}

HEM.topSNPs <- rbind(HEM.LA410, HEM.LA0480, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)
Top50SNP <- HEM.topSNPs
#two greys look the same. add purple.
myColors <- c("#999999", "#292929","#684800" ,"#CBA22A", "#63B2D3", "#1FA69D", "#57B761", "#DAD94C","#2B869D", "#EE82EE", "#D2652D", "#CC79A7")
names(myColors) <- levels(Top50SNP$Trait)
colScale <- scale_colour_manual(name = "Plant",values = myColors)

plot1 <- ggplot(Top50SNP, aes(x=Index, y=Effect.s))
jpeg("paper/plots/ActualPaper/FigR6/FigR6_largeFxPlantSNPs_NA10_top200.jpg", width=8, height=4, units='in', res=600)
plot1 + geom_point(aes(color=factor(Trait)), size=3, alpha=3/4)+ colScale+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_y_continuous(name="Relative SNP Effect Estimate")+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
dev.off()
