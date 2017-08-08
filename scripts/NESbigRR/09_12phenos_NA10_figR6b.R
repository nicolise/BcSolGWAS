#Nicole E Soltis
#Summary figure for BcSlGWAS 
#081916

#input: Sl_LesionSize_MAF20.HEM.csv
#plot: FigR7_Summary_99Thresh_ManhattanPlot.jpg
#note that plot: FigR6_topSNPsoverlap is in 08_peaks script
library(tidyr)
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
setwd("~/Documents/GitRepos/BcSolGWAS")
SNPlist  <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")
Thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_12plants_trueMAF20_10NA.HEM.Thresh.final.csv")

names(SNPlist)
SNPlist <- SNPlist[,c(2:14,16,17,15)]
names(Thresh)
Thresh <- Thresh[,c(2:14)]

#positives
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
write.csv(HOSNP, "data/GWAS_files/05_annotation/TrueMAF_NAs/12Plants_HiOverlapSNPs_trueMAF20_10NA.csv")

SUMM.plot <- SNPbin
#draw the plots!!!
#write.csv(SUMM.plot, "data/genome/SummaryManhattanPlot_data.csv")

#create a custom color scale
myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
#names(myColors) <- levels(HEM.plotdata$Chrom)
library(ggplot2)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#make plots 
#for poster figures, width=8, height=4
 jpeg("paper/plots/ActualPaper/FigR7/FigR7a_Summary_99Thresh_ManhattanPlot_NA10.jpg", width=7.5, height=5, units='in', res=600)
#SUMMtemp <- subset(SUMM.plot[SUMM.plot$Chrom==c(5,6),])
  ggplot(SUMM.plot, aes(x=Index, y=SUMMpos))+
    colScale+ #remove for rainbow plot
    theme_bw()+
#    scale_x_continuous(breaks = ticks)+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="Plant Genotypes per Significant SNP", x="Chromosome position"))+
    #nrow=8
    theme(legend.position="none")+
    guides(col = guide_legend(nrow = 8, title=element_blank()))+
    scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    scale_y_continuous(breaks= c(0,2,4,6,8,10,12), labels=c("0","2","4","6","8","10","12"))
  dev.off()

SUMMhi <- subset(SUMM.plot, SUMM.plot$SUMM == 1)

subset(HEM.plotdata, HEM.plotdata[i] > get(paste("TH999pos_", names(HEM.plotdata[i]), sep="")), select=c(Chrom,Segment,Pos,Index,i))
  
names(SUMM.plot)
TopSNPs <- SUMM.plot[which(SUMM.plot$SUMM > 1),]
#write.csv(TopSNPs, "data/TopSNPs_ALLtomato.csv"

#now plot with only chr 1 // chr 2

jpeg("paper/plots/ActualPaper/FigR7/FigR7a_Summary_99Thresh_ManhattanPlot_NA10_Chr9.jpg", width=7.5, height=5, units='in', res=600)
#SUMMtemp <- subset(SUMM.plot[SUMM.plot$Chrom==c(5,6),])
ggplot(SUMM.plot, aes(x=Index, y=SUMM))+
  colScale+ #remove for rainbow plot
  theme_bw()+
  #    scale_x_continuous(breaks = ticks)+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="Number of Significant SNPs Across Plants", x="Chromosome position"))+
  #nrow=8
  theme(legend.position="none")+
  guides(col = guide_legend(nrow = 8, title=element_blank()))+
  #name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"),
  scale_x_continuous(limits=c(23192915,25584347))+
  scale_y_continuous(breaks= c(0,2,4,6,8,10,12), labels=c("0","2","4","6","8","10","12"))
dev.off()

max(SUMM.plot[ which(SUMM.plot$Chrom=='9'),]$Index)
min(SUMM.plot[ which(SUMM.plot$Chrom=='9'),]$Index)
