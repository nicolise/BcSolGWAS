#Nicole E Soltis
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
#setwd("~/Projects/BcSolGWAS/")

#Input file: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv AND Sl_DomesticationLS_MAF20.HEM.Thresh.csv
#Output file: NONE
#Plots: Basic and greyscale domestication Manhattan plots
#go to script 09_bigRR_meta for domestication color plot (fig R8)
###########################################################################
#Plotting the HEM results

#Load plotting package
library(ggplot2); library(grid); library(plyr)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.PlotFormat.csv")
HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1)]

TH99pos <- HEM.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}

TH999pos <- HEM.thresh[4,]
for (i in 2:ncol(TH999pos)){
  assign(paste("TH999pos_", names(TH999pos[i]), sep=""),as.numeric(TH999pos[i]))
}

TH95pos <- HEM.thresh[1,]
for (i in 2:ncol(TH95pos)){
  assign(paste("TH95pos_", names(TH95pos[i]), sep=""),as.numeric(TH95pos[i]))
}

TH99neg <- HEM.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}

TH999neg <- HEM.thresh[8,]
for (i in 2:ncol(TH999neg)){
  assign(paste("TH999neg_", names(TH999neg[i]), sep=""),as.numeric(TH999neg[i]))
}

TH95neg <- HEM.thresh[5,]
for (i in 2:ncol(TH95neg)){
  assign(paste("TH95neg_", names(TH95neg[i]), sep=""),as.numeric(TH95neg[i]))
}

#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#this is the correct list for NA10
#c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579)

#this is the correct list for NA20
#c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  9006772)

#lowTR
for (i in 4:6){
jpeg(paste("plots/paper/SlBc_trueMAF20_10NA_lowTR_manhattan",names(HEM.plotdata[i]),".jpg", sep=""), width=8, height=5, units='in', res=600)
print(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
  theme_bw()+
  colScale+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=2) +
  geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=2) +
  geom_text(aes(0,get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", vjust = 1.2, hjust = .05), col = "black")+
  geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
  geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
  geom_text(aes(0,get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
dev.off()
}

#high TR
for (i in 4:6){
jpeg(paste("plots/paper/SlBc_trueMAF20_10NA_highTR_manhattan",names(HEM.plotdata[i]),".jpg", sep=""), width=8, height=4, units='in', res=600)
plot(ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[,i]))+
  theme_bw()+
  colScale+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title="Phenotype"))+
  geom_hline(yintercept=get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
  geom_hline(yintercept=(get(paste("TH999neg_", names(HEM.plotdata[i]), sep=""))), lty=2) +
  geom_text(aes(-0.01,get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), label = "99.9% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
theme(legend.position="none"))
dev.off()
}

#meta analysis plot
jpeg("plots/MultiPlot/trueMAF/SlBc_trueMAF20_highTR_metaplot.jpg", width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index))+
  theme_bw()+
  #colScale+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title="Phenotype"))+
  geom_point(aes(x=Index, y=(Domesticated), color = "Domesticated"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Wild), color = "Wild"), alpha=1/2)+
  geom_point(aes(x=Index, y=(DmWoD), color = "(D-W)/D"), alpha=1/2)+
  geom_hline(yintercept=get(paste("TH999pos_", names(HEM.plotdata[4]), sep=""))) +
  geom_hline(yintercept=(get(paste("TH999neg_", names(HEM.plotdata[4]), sep="")))) +
  geom_text(aes(-0.01,get(paste("TH999pos_", names(HEM.plotdata[4]), sep="")), label = "99.9% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
#theme(legend.position="none")
dev.off()
