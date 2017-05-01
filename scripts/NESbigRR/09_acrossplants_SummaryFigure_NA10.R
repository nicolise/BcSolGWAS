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

#just postitives
Threshval <- Thresh[7,]
names(Threshval)
SNPbin <- SNPlist
#I'll do this with a loop over the columns later
#fixing this to include negative values
#positives
Thresh_1547 <- Threshval[,"LA1547"]
SNPbin$LA1547pos <- ifelse(SNPbin$LA1547 > Thresh_1547, 1, 0)
Thresh_1589 <- Threshval[,"LA1589"]
SNPbin$LA1589pos <- ifelse(SNPbin$LA1589 > Thresh_1589, 1, 0)
Thresh_1684 <- Threshval[,"LA1684"]
SNPbin$LA1684pos <- ifelse(SNPbin$LA1684 > Thresh_1684, 1, 0)
Thresh_2093 <- Threshval[,"LA2093"]
SNPbin$LA2093pos <- ifelse(SNPbin$LA2093 > Thresh_2093, 1, 0)
Thresh_2176 <- Threshval[,"LA2176"]
SNPbin$LA2176pos <- ifelse(SNPbin$LA2176 > Thresh_2176, 1, 0)
Thresh_2706 <- Threshval[,"LA2706"]
SNPbin$LA2706pos <- ifelse(SNPbin$LA2706 > Thresh_2706, 1, 0)
Thresh_3008 <- Threshval[,"LA3008"]
SNPbin$LA3008pos <- ifelse(SNPbin$LA3008 > Thresh_3008, 1, 0)
Thresh_3475 <- Threshval[,"LA3475"]
SNPbin$LA3475pos <- ifelse(SNPbin$LA3475 > Thresh_3475, 1, 0)
Thresh_410 <- Threshval[,"LA410"]
SNPbin$LA410pos <- ifelse(SNPbin$LA410 > Thresh_410, 1, 0)
Thresh_4345 <- Threshval[,"LA4345"]
SNPbin$LA4345pos <- ifelse(SNPbin$LA4345 > Thresh_4345, 1, 0)
Thresh_4355 <- Threshval[,"LA4355"]
SNPbin$LA4355pos <- ifelse(SNPbin$LA4355 > Thresh_4355, 1, 0)
Thresh_0480 <- Threshval[,"LA0480"]
SNPbin$LA0480pos <- ifelse(SNPbin$LA0480 > Thresh_0480, 1, 0)
#negatives
#but if the positive > 0, force negative to zero.
Threshval <- Thresh[3,]
names(Threshval)
Thresh_1547 <- Threshval[,"LA1547"]
SNPbin$LA1547neg <- ifelse(SNPbin$LA1547pos > 0, 0, ifelse (SNPbin$LA1547 < Thresh_1547, 1, 0))
Thresh_1589 <- Threshval[,"LA1589"]
SNPbin$LA1589neg <- ifelse(SNPbin$LA1589pos > 0, 0, ifelse (SNPbin$LA1589 < Thresh_1589, 1, 0))
Thresh_1684 <- Threshval[,"LA1684"]
SNPbin$LA1684neg <- ifelse(SNPbin$LA1684pos > 0, 0, ifelse (SNPbin$LA1684 < Thresh_1684, 1, 0))
Thresh_2093 <- Threshval[,"LA2093"]
SNPbin$LA2093neg <- ifelse(SNPbin$LA2093pos > 0, 0, ifelse (SNPbin$LA2093 < Thresh_2093, 1, 0))
Thresh_2176 <- Threshval[,"LA2176"]
SNPbin$LA2176neg <- ifelse(SNPbin$LA2176pos > 0, 0, ifelse (SNPbin$LA2176 < Thresh_2176, 1, 0))
Thresh_2706 <- Threshval[,"LA2706"]
SNPbin$LA2706neg <- ifelse(SNPbin$LA2706pos > 0, 0, ifelse (SNPbin$LA2706 < Thresh_2706, 1, 0))
Thresh_3008 <- Threshval[,"LA3008"]
SNPbin$LA3008neg <- ifelse(SNPbin$LA3008pos > 0, 0, ifelse (SNPbin$LA3008 < Thresh_3008, 1, 0))
Thresh_3475 <- Threshval[,"LA3475"]
SNPbin$LA3475neg <- ifelse(SNPbin$LA3475pos > 0, 0, ifelse (SNPbin$LA3475 < Thresh_3475, 1, 0))
Thresh_410 <- Threshval[,"LA410"]
SNPbin$LA410neg <- ifelse(SNPbin$LA410pos > 0, 0, ifelse (SNPbin$LA410 < Thresh_410, 1, 0))
Thresh_4345 <- Threshval[,"LA4345"]
SNPbin$LA4345neg <- ifelse(SNPbin$LA4345pos > 0, 0, ifelse (SNPbin$LA4345 < Thresh_4345, 1, 0))
Thresh_4355 <- Threshval[,"LA4355"]
SNPbin$LA4355neg <- ifelse(SNPbin$LA4355pos > 0, 0, ifelse (SNPbin$LA4355 < Thresh_4355, 1, 0))
Thresh_0480 <- Threshval[,"LA0480"]
SNPbin$LA0480neg <- ifelse(SNPbin$LA0480pos > 0, 0, ifelse (SNPbin$LA0480 < Thresh_0480, 1, 0))

#now create a summation column
names(SNPbin)
SNPbin$SUMMpos <- rowSums(SNPbin[,c(17:28)])
SNPbin$SUMMneg <- rowSums(SNPbin[,c(29:40)])
SNPbin$SUMMneg
SNPbin$SUMM <- SNPbin$SUMMneg + SNPbin$SUMMpos

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
  ggplot(SUMM.plot, aes(x=Index, y=SUMM))+
    colScale+ #remove for rainbow plot
    theme_bw()+
#    scale_x_continuous(breaks = ticks)+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="Number of Significant SNPs Across Plants", x="Chromosome position"))+
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
