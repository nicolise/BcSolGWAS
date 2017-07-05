#Nicole E Soltis 
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input File: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv and .Thresh.csv
#Output File: results/Domestication_TopSNPs_SegLong.csv, results/Domestication_TopSNPs_SegWide.csv
#this goes into gene annotation and then venn diagrams
#Plots: NONE
############################################################################

#Load packages
library(plyr); library(ggplot2); library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")

HEM.plotdata <- HEM.plotdata[,c(2:14,16,17,15)]

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
#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in c(4:15)){
assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), subset(HEM.plotdata, HEM.plotdata[i] < get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i)))
}

# #for top 1000 only
# for (i in c(4:15)){
#   assign(paste("HEMpos.", names(HEM.plotdata[i]), sep=""), head(arrange(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMpos.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
#   assign(paste("HEMneg.", names(HEM.plotdata[i]), sep=""), tail(arrange(get(paste("HEMneg.", names(HEM.plotdata[i]), sep="")), desc(get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))[,5])), n=500))
# }

#combine pos and neg by group
for (i in c(4:15)){
assign(paste("HEM.", names(HEM.plotdata[i]), sep=""), rbind(get(paste("HEMpos.", names(HEM.plotdata[i]), sep="")),get(paste("HEMneg.", names(HEM.plotdata[i]), sep=""))))
}

#then combine
for (i in c(4:15)){
mydf <- paste("HEM.", names(HEM.plotdata[i]), sep="")
renamedf <- get(mydf)
colnames(renamedf)[5] <- "Effect"
assign(mydf, renamedf)
myblob <- rep(names(HEM.plotdata[i]), nrow(get(mydf)))
assign(mydf, cbind(get(mydf), Trait = myblob))
}

HEM.topSNPs <- rbind(HEM.LA410, HEM.LA0480, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)

#total sig SNPs per trait here
table(HEM.topSNPs$Trait)



#omit LA0480 because it's so high compared to the rest
HEM.topplots <- rbind(HEM.LA410, HEM.LA1547, HEM.LA1589, HEM.LA1684, HEM.LA2093, HEM.LA2176, HEM.LA2706, HEM.LA3008, HEM.LA3475, HEM.LA4345, HEM.LA4355)
library(ggplot2)
plot1 <- ggplot(HEM.topplots, aes(x=Index, y=Effect))
plot1 + geom_point(aes(color=factor(Trait)))+ theme_bw()

names(HEM.topplots)

#have to do this with top 1000 per trait only
TopSNP.wide.DM <- reshape(HEM.topSNPs, 
                          timevar = "Trait",
                          idvar = c("Chrom","Segment","Pos","Index"),
                          direction = "wide")

#count total number of traits per SNP
TopSNP.wide.DM[is.na(TopSNP.wide.DM)] <- 0
TopSNP.wide.DM$TotPos <- rowSums(TopSNP.wide.DM[,c(5:16)] > 0)
TopSNP.wide.DM$TotNeg <- rowSums(TopSNP.wide.DM[,c(5:16)] < 0)
TopSNP.wide.DM$TotTraits <- TopSNP.wide.DM$TotNeg + TopSNP.wide.DM$TotPos

hist(TopSNP.wide.DM$TotTraits)
table(TopSNP.wide.DM$TotTraits)

jpeg("paper/plots/ActualPaper/FigR7/R7a_topSNPssOverlap_12Plants.jpg", width=8, height=5, units='in', res=600)
ggplot(TopSNP.wide.DM, aes(TopSNP.wide.DM$TotTraits)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of SNPs")+
  scale_x_continuous(name= "Plant Genotypes per Candidate SNP", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12), limits = c(0, 12))
dev.off()

#make it wide format
#currently long format : Chrom, Segment, Pos, Index, Effect, Trait
#write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")
write.csv(HEM.plotdata, "data/GWAS_files/05_annotation/TrueMAF_NAs/12Plants_Top1000SNPs_SegLong_trueMAF20_10NA.csv")

write.csv(TopSNP.wide.DM, "data/GWAS_files/05_annotation/12Plants_Top1000SNPs_SegWide_trueMAF20_10NA.csv")
