#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
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

HEM.plotdata <- HEM.plotdata.og

HEM.plotdata <- HEM.plotdata[,-c(1)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.Thresh.final.csv")
HEM.thresh <- HEM.thresh[,-c(1)]

#average 99% threshold
mean(6.674198e-06, 1.493645e-05, 1.290865e-05, 8.490481e-06, 2.055714e-05, 1.492797e-05, 1.714183e-05, 1.756520e-05, 1.212998e-05, 1.376836e-05, 1.534867e-05, 1.443944e-05, 6.973896e-06, 1.520926e-05, 1.309197e-05, 8.455100e-06, 2.134363e-05, 1.549896e-05, 1.698335e-05, 1.761662e-05, 1.283858e-05, 1.418813e-05, 1.574277e-05, 1.534149e-05)

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
#((max(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index) - min(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index))/2+min(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index))

#get length per chromosome segment
#max(HEM.plotdata[which(HEM.plotdata$Chrom.Seg.F=='16.7'),]$Index) - min(HEM.plotdata[which(HEM.plotdata$Chrom.Seg.F=='16.7'),]$Index)

#this is the correct list for NA10
#c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579)

#this is the correct list for NA20
#c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682,  38915229)

df.hlt <- data.frame(HEM.plotdata[HEM.plotdata$Index %in% c(1099438,1836245,31154483,33853054,36555407,41350409,5703231,6589294,7955289,11188054,22332692,24530790),]) 
df.hlt <- df.hlt[,c("Index", "LA2093")]
df.hlt$myY <- 100*df.hlt$LA2093
df.hlt <- df.hlt[,c(1,3)]
names(df.hlt)[2] <- "myY"

HEM.plotdata$myY <- 100*HEM.plotdata[,7]

#greyscale version
#4 to 15
#for (i in c(7)){
i <- 7
  jpeg(paste("paper/Submissions/PlantCell/TPC revision/Figures/Fig4a_bw_Sl_LesionSize_trueMAF20_NA10_lowTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=7.5, height=5, units='in', res=600)
  plot(
    ggplot(HEM.plotdata, aes(x=Index, y=myY))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Chrom)))+
         labs(list(y=expression(paste("Estimated Effect Size (",mm^{2},")")), title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
         theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
         guides(col = guide_legend(nrow = 8, title="Chromosome"))+
         geom_hline(yintercept=100*get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
         geom_hline(yintercept=100*get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
         geom_text(aes(0,100*get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", size=14, vjust = 2, hjust=.05), col = "black")+
         geom_hline(yintercept=100*get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_hline(yintercept=100*get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
         geom_text(aes(0,100*get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", size=14, vjust = 4.5, hjust=.05), col = "black")+
         theme(legend.position="none")+
         theme(panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         #NA10 chromosomes
         scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
         geom_vline(xintercept=1099438, lty=2)+
         geom_vline(xintercept=1836245, lty=2)+
         geom_vline(xintercept=31154483, lty=2)+
         geom_vline(xintercept=33853054, lty=2)+
         geom_vline(xintercept=36555407, lty=2)+
         geom_vline(xintercept=41350409, lty=2)+
         geom_vline(xintercept=5703231, lty=2)+
         geom_vline(xintercept=6589294, lty=2)+
         geom_vline(xintercept=7955289, lty=2)+
         geom_vline(xintercept=11188054, lty=2)+
         geom_vline(xintercept=22332692, lty=2)+
         geom_vline(xintercept=24530790, lty=2)+
         geom_point(data=df.hlt,  pch=21, colour="black")
    #
  )
  dev.off()
#}

#NA20 chromosomes
#scale_x_continuous(name="Chromosome", breaks = c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682, 38915229), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+

#highTR BW
  #plot 2706 as example, [9]
  #with labeling removed for paper
#jpeg(paste("paper/plots/ActualPaper/FigR5/Routs/pm_BW_Sl_MAF20_highTR_", names(HEM.plotdata[8]), ".ManhattanPlot.jpg", sep=""), width=7.5, height=4, units='in', res=600)
#4 to 15
for (i in c(15)){
  jpeg(paste("plots/paper/bw_Sl_LesionSize_trueMAF20_NA10_hiTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
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
  scale_x_continuous(name="Chromosome", breaks = c(1677869, 5250031, 9006772, 11066464, 13584641, 17192569, 20021437, 22387756, 24411342, 26784999, 28587689, 30133032, 31892533, 34009041, 35807682, 38915229), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0))
#dev.off()
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