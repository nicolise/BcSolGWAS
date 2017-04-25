#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
############################################################################
###Plotting the HEM results

#NEED TO CHECK CONTIGS FOR THIS

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_12plants_trueMAF20_20NA.HEM.PlotFormat.csv")

HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_12plants_trueMAF20_20NA.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1)]

TH95pos <- HEM.thresh[1,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}
TH95neg <- HEM.thresh[5,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}
TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}
TH999pos <- HEM.thresh[4,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}
TH999neg <- HEM.thresh[8,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

#create a custom color scale
myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#get midpoint positions per chromosome
((max(HEM.plotdata[ which(HEM.plotdata$Chrom=='14'),]$Index) - min(HEM.plotdata[ which(HEM.plotdata$Chrom=='14'),]$Index))/2+min(HEM.plotdata[ which(HEM.plotdata$Chrom=='14'),]$Index))

#get length per chromosome segment
max(HEM.plotdata[which(HEM.plotdata$Chrom.Seg.F=='16.7'),]$Index) - min(HEM.plotdata[which(HEM.plotdata$Chrom.Seg.F=='16.7'),]$Index)

#this is the correct list for NA10
#c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579)

#this is the correct list for NA20
#c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  9006772)

#greyscale version
#4 to 15
for (i in 7){
  #jpeg(paste("paper/plots/ActualPaper/bw_Sl_LesionSize_trueMAF20_NA10_lowTR_", names(HEM.plotdata[7]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
jpeg(paste("plots/paper/bw_Sl_LesionSize_trueMAF20_NA20_lowTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
  plot(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="SNP Effect Estimate", title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
    geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
    geom_text(aes(0,get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", vjust = 1.2, hjust=.05), col = "black")+
   geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
    geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
   geom_text(aes(0,get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", vjust = 1.5, hjust=.05), col = "black")+
    theme(legend.position="none")+
    #NA20 chromosomes
    scale_x_continuous(name="Chromosome", breaks = c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  9006772), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    #NA10 chromosomes
   # scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    expand_limits(y=0))
  dev.off()
}

#highTR BW
  #plot 2706 as example, [9]
  #with labeling removed for paper
#jpeg(paste("paper/plots/ActualPaper/FigR5/Routs/pm_BW_Sl_MAF20_highTR_", names(HEM.plotdata[8]), ".ManhattanPlot.jpg", sep=""), width=7.5, height=4, units='in', res=600)
#4 to 15
for (i in 4:15){
  jpeg(paste("plots/paper/bw_Sl_LesionSize_trueMAF20_NA20_hiTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
  plot(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="SNP Effect Estimate", title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    geom_hline(yintercept=get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
    geom_hline(yintercept=get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
    geom_text(aes(0,get(paste("TH999neg_", names(HEM.plotdata[i]), sep="")), label = "99.9% Threshold", vjust = 1.2, hjust = .05), col = "black")+
    theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  9006772), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
dev.off()
}
#stop [i]

#color highTR
for (i in 4:15){
jpeg(paste("plots/MultiPlot/NewModel0711b/color_Sl_LesionSize_MAF20_highTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
plot(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[i]))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH999_", names(HEM.plotdata[i]), sep=""))) +
  geom_hline(yintercept=-0.0004520958)+
  geom_text(aes(0,get(paste("TH999_", names(HEM.plotdata[i]), sep="")), label =
  "99.9% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
dev.off()
}

#color lowTR
for (y in 4:15){
  jpeg(paste("plots/MultiPlot/NewModel0711b/Sl_LesionSize_MAF20_lowTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
  plot(ggplot(HEM.plotdata, aes(x=Index, y=abs(HEM.plotdata[i])))+
    theme_bw()+
    colScale+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
    geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
    geom_text(aes(0,get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", vjust = 1.2, hjust = .05), col = "black")+
    geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[i]), sep="", lty=2))) +
    geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[i]), sep="", lty=2))) +
    geom_text(aes(0,get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", vjust = 1.2, hjust = .05), col = "black")+
    expand_limits(y=0))
  dev.off()}