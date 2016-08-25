#Summary figure for BcSlGWAS 
#081916
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
SNPlist <- read.csv("data/GWAS_files/04_bigRRoutput/NewModel0711/Sl_LesionSize_MAF20.HEM.PlotFormat.csv")
Thresh <- read.csv("data/GWAS_files/04_bigRRoutput/NewModel0711/Sl_LesionSize_MAF20.HEM.Thresh.csv")

names(SNPlist)
SNPlist <- SNPlist[,c(3:16)]
names(Thresh)
Thresh <- Thresh[,c(3:15)]

#0.999 Threshold: LA0410, LA0480, LA1589, LA1684, LA2093, LA2176, LA2706
#SNPlist2 <- SNPlist[,c("Chrom", "Pos", "LA410", "LA480", "LA1589", "LA1684", "LA2093", "LA2176", "LA2706")]
#but can do with all. Will just get all zeros for the others.


Threshval <- Thresh[4,]
names(Threshval)

SNPbin <- SNPlist
#I'll do this with a loop over the columns later
Thresh_1547 <- Threshval[,"LA1547"]
SNPbin$LA1547 <- ifelse(SNPlist$LA1547 > Thresh_1547, 1, 0)
Thresh_1589 <- Threshval[,"LA1589"]
SNPbin$LA1589 <- ifelse(SNPlist$LA1589 > Thresh_1589, 1, 0)
Thresh_1684 <- Threshval[,"LA1684"]
SNPbin$LA1684 <- ifelse(SNPlist$LA1684 > Thresh_1684, 1, 0)
Thresh_2093 <- Threshval[,"LA2093"]
SNPbin$LA2093 <- ifelse(SNPlist$LA2093 > Thresh_2093, 1, 0)
Thresh_2176 <- Threshval[,"LA2176"]
SNPbin$LA2176 <- ifelse(SNPlist$LA2176 > Thresh_2176, 1, 0)
Thresh_2706 <- Threshval[,"LA2706"]
SNPbin$LA2706 <- ifelse(SNPlist$LA2706 > Thresh_2706, 1, 0)
Thresh_3008 <- Threshval[,"LA3008"]
SNPbin$LA3008 <- ifelse(SNPlist$LA3008 > Thresh_3008, 1, 0)
Thresh_3475 <- Threshval[,"LA3475"]
SNPbin$LA3475 <- ifelse(SNPlist$LA3475 > Thresh_3475, 1, 0)
Thresh_410 <- Threshval[,"LA410"]
SNPbin$LA410 <- ifelse(SNPlist$LA410 > Thresh_410, 1, 0)
Thresh_4345 <- Threshval[,"LA4345"]
SNPbin$LA4345 <- ifelse(SNPlist$LA4345 > Thresh_4345, 1, 0)
Thresh_4355 <- Threshval[,"LA4355"]
SNPbin$LA4355 <- ifelse(SNPlist$LA4355 > Thresh_4355, 1, 0)
Thresh_480 <- Threshval[,"LA480"]
SNPbin$LA480 <- ifelse(SNPlist$LA480 > Thresh_480, 1, 0)

#now create a summation column
names(SNPbin)
SNPbin$SUMM <- rowSums(SNPbin[,c(3:14)])
SNPbin$SUMM

SUMM.plot <- SNPbin
#draw the plots!!!
# #Reformat Chromosomes and Positions
SUMM.plot$Chrom <- gsub("Chromosome", "", SUMM.plot$Chrom)
SUMM.plot$Chrom <- as.numeric(as.character(SUMM.plot$Chrom))
SUMM.plot$Pos <- as.numeric(as.character(SUMM.plot$Pos))

#sort dataframe rows in order of Chrom, then Pos
SUMM.plot <- SUMM.plot[with(SUMM.plot, order(Chrom, Pos)), ]

#Make plotting variables
SUMM.plot$Index = NA
ticks = NULL
lastbase = 0

#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(SUMM.plot$Chrom)) {
  print(i)
  if (i==1) {
    SUMM.plot[SUMM.plot$Chrom==i, ]$Index=SUMM.plot[SUMM.plot$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=+lastbase+max(subset(SUMM.plot,SUMM.plot$Chrom==i-1)$Pos, 1)
    SUMM.plot[SUMM.plot$Chrom==i, ]$Index=SUMM.plot[SUMM.plot$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, SUMM.plot[SUMM.plot$Chrom==i, ]$Index[floor(length(SUMM.plot[SUMM.plot$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(SUMM.plot$Index),max(SUMM.plot$Index))

#make plots 
library(ggplot2)
  jpeg("plots/MultiPlot/Summary_999Thresh_ManhattanPlot.jpg", width=8, height=4, units='in', res=600)
  ggplot(SUMM.plot, aes(x=Index, y=SUMM))+
    theme_bw()+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="Frequency of SNPs over Threshold (99.9%)", x="Genome position"))+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))
  dev.off()
  
names(SUMM.plot)
TopSNPs <- SUMM.plot[which(SUMM.plot$SUMM > 1),]
write.csv(TopSNPs, "data/TopSNPs_ALLtomato.csv")
