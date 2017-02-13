#Nicole E Soltis
#Summary figure for BcSlGWAS 
#081916

#input: Sl_LesionSize_MAF20.HEM.csv
#plot: FigR6_Summary_999Thresh_ManhattanPlot.jpg
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
setwd("~/Documents/GitRepos/BcSolGWAS")
SNPlist <- read.csv("data/GWAS_files/04_bigRRoutput/NewModel0711/Sl_LesionSize_MAF20.HEM.csv")

Thresh <- SNPlist[c(1:4),]
SNPlist <- SNPlist[-c(1:4),]
SNPlist2 <- SNPlist

library(tidyr)

#first, need to change Chromosome1.number to Chromosome1.0.number etc.

#I'll change everything to 1.0.(1.)number

SNPlist2$X.1 <- gsub(pattern = "Chromosome1\\.", replacement = "Chromosome1.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome2\\.", replacement = "Chromosome2.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome3\\.", replacement = "Chromosome3.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome4\\.", replacement = "Chromosome4.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome5\\.", replacement = "Chromosome5.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome6\\.", replacement = "Chromosome6.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome7\\.", replacement = "Chromosome7.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome8\\.", replacement = "Chromosome8.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome9\\.", replacement = "Chromosome9.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome10\\.", replacement = "Chromosome10.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome11\\.", replacement = "Chromosome11.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome12\\.", replacement = "Chromosome12.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome13\\.", replacement = "Chromosome13.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome14\\.", replacement = "Chromosome14.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome15\\.", replacement = "Chromosome15.0.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.", replacement = "Chromosome16.0.", SNPlist2$X.1)

#Then everything with 1.0.(1.)number I'll change to 1.(1.)number

SNPlist2$X.1 <- gsub(pattern = "Chromosome1\\.0\\.1\\.", replacement = "Chromosome1.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome2\\.0\\.1\\.", replacement = "Chromosome2.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome2\\.0\\.2\\.", replacement = "Chromosome2.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome2\\.0\\.3\\.", replacement = "Chromosome2.3.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome2\\.0\\.4\\.", replacement = "Chromosome2.4.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome2\\.0\\.5\\.", replacement = "Chromosome2.5.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome3\\.0\\.1\\.", replacement = "Chromosome3.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome3\\.0\\.2\\.", replacement = "Chromosome3.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome4\\.0\\.1\\.", replacement = "Chromosome4.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome5\\.0\\.1\\.", replacement = "Chromosome5.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome6\\.0\\.1\\.", replacement = "Chromosome6.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome6\\.0\\.2\\.", replacement = "Chromosome6.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome6\\.0\\.3\\.", replacement = "Chromosome6.3.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome7\\.0\\.1\\.", replacement = "Chromosome7.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome7\\.0\\.2\\.", replacement = "Chromosome7.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome8\\.0\\.1\\.", replacement = "Chromosome8.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome8\\.0\\.2\\.", replacement = "Chromosome8.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome9\\.0\\.1\\.", replacement = "Chromosome9.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome10\\.0\\.1\\.", replacement = "Chromosome10.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome11\\.0\\.1\\.", replacement = "Chromosome11.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome12\\.0\\.1\\.", replacement = "Chromosome12.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome13\\.0\\.1\\.", replacement = "Chromosome13.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome13\\.0\\.2\\.", replacement = "Chromosome13.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome14\\.0\\.1\\.", replacement = "Chromosome14.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome14\\.0\\.2\\.", replacement = "Chromosome14.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome15\\.0\\.1\\.", replacement = "Chromosome15.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome15\\.0\\.2\\.", replacement = "Chromosome15.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome15\\.0\\.3\\.", replacement = "Chromosome15.3.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome15\\.0\\.4\\.", replacement = "Chromosome15.4.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.1\\.", replacement = "Chromosome16.1.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.2\\.", replacement = "Chromosome16.2.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.3\\.", replacement = "Chromosome16.3.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.4\\.", replacement = "Chromosome16.4.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.5\\.", replacement = "Chromosome16.5.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.6\\.", replacement = "Chromosome16.6.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.7\\.", replacement = "Chromosome16.7.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.8\\.", replacement = "Chromosome16.8.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.9\\.", replacement = "Chromosome16.9.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.10\\.", replacement = "Chromosome16.10.", SNPlist2$X.1)
SNPlist2$X.1 <- gsub(pattern = "Chromosome16\\.0\\.11\\.", replacement = "Chromosome16.11.", SNPlist2$X.1)

SNPlist3 <- separate (SNPlist2, X.1, into = c("Chrom", "Segment", "Pos") )

SNPlist <- SNPlist3
names(SNPlist)
SNPlist <- SNPlist[,c(2:16)]
names(Thresh)
Thresh <- Thresh[,c(2:14)]

#0.999 Threshold: LA0410, LA0480, LA1589, LA1684, LA2093, LA2176, LA2706
#SNPlist2 <- SNPlist[,c("Chrom", "Pos", "LA410", "LA480", "LA1589", "LA1684", "LA2093", "LA2176", "LA2706")]
#but can do with all. Will just get all zeros for the others.

Threshval <- Thresh[4,]
names(Threshval)

SNPbin <- SNPlist
#I'll do this with a loop over the columns later
#fixing this to include negative values
Thresh_1547 <- Threshval[,"LA1547"]
SNPbin$LA1547 <- ifelse(abs(SNPlist$LA1547) > Thresh_1547, 1, 0)
Thresh_1589 <- Threshval[,"LA1589"]
SNPbin$LA1589 <- ifelse(abs(SNPlist$LA1589) > Thresh_1589, 1, 0)
Thresh_1684 <- Threshval[,"LA1684"]
SNPbin$LA1684 <- ifelse(abs(SNPlist$LA1684) > Thresh_1684, 1, 0)
Thresh_2093 <- Threshval[,"LA2093"]
SNPbin$LA2093 <- ifelse(abs(SNPlist$LA2093) > Thresh_2093, 1, 0)
Thresh_2176 <- Threshval[,"LA2176"]
SNPbin$LA2176 <- ifelse(abs(SNPlist$LA2176) > Thresh_2176, 1, 0)
Thresh_2706 <- Threshval[,"LA2706"]
SNPbin$LA2706 <- ifelse(abs(SNPlist$LA2706) > Thresh_2706, 1, 0)
Thresh_3008 <- Threshval[,"LA3008"]
SNPbin$LA3008 <- ifelse(abs(SNPlist$LA3008) > Thresh_3008, 1, 0)
Thresh_3475 <- Threshval[,"LA3475"]
SNPbin$LA3475 <- ifelse(abs(SNPlist$LA3475) > Thresh_3475, 1, 0)
Thresh_410 <- Threshval[,"LA410"]
SNPbin$LA410 <- ifelse(abs(SNPlist$LA410) > Thresh_410, 1, 0)
Thresh_4345 <- Threshval[,"LA4345"]
SNPbin$LA4345 <- ifelse(abs(SNPlist$LA4345) > Thresh_4345, 1, 0)
Thresh_4355 <- Threshval[,"LA4355"]
SNPbin$LA4355 <- ifelse(abs(SNPlist$LA4355) > Thresh_4355, 1, 0)
Thresh_480 <- Threshval[,"LA480"]
SNPbin$LA480 <- ifelse(abs(SNPlist$LA480) > Thresh_480, 1, 0)

#now create a summation column
names(SNPbin)
SNPbin$SUMM <- rowSums(SNPbin[,c(4:15)])
SNPbin$SUMM

SUMM.plot <- SNPbin
#draw the plots!!!
# #Reformat Chromosomes and Positions
SUMM.plot$Chrom <- gsub("Chromosome", "", SUMM.plot$Chrom)
SUMM.plot$Chrom <- as.numeric(as.character(SUMM.plot$Chrom))
SUMM.plot$Segment <- as.numeric(as.character(SUMM.plot$Segment))
SUMM.plot$Pos <- as.numeric(as.character(SUMM.plot$Pos))

#sort dataframe rows in order of Chrom, then Pos
SUMM.plot <- SUMM.plot[with(SUMM.plot, order(Chrom, Segment, Pos)), ]
#now make segments line up consecutively
SUMM.plot$Chrom.Seg <- paste(SUMM.plot$Chrom, SUMM.plot$Segment, sep=".")
SUMM.plot$Chrom.Seg <- as.numeric(SUMM.plot$Chrom.Seg)

write.csv(SUMM.plot, "data/genome/SummaryManhattanPlot_data.csv")

#let's try making the chrom.seg integers so that R isn't confused
unique(SUMM.plot$Chrom.Seg)
data$scode <- revalue(data$sex, c("M"="1", "F"="2"))
library(plyr)
SUMM.plot$Chrom.Seg.F <- as.factor(SUMM.plot$Chrom.Seg)
unique(SUMM.plot$Chrom.Seg.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6,  6.1, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 13.2, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

SUMM.plot$Chrom.Seg.Int <- recode.vars$newvals[match(SUMM.plot$Chrom.Seg.F, recode.vars$OGvals)]
unique(SUMM.plot$Chrom.Seg.Int)

#Make plotting variables
SUMM.plot$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing

for (i in unique(SUMM.plot$Chrom.Seg.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of SUMM.plot rows with Chromosome 1, set Index variable for each row to equal Pos.
    SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Index=SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
      #current lastbase counter plus the maxiumum position of chromosome i-1
      #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(SUMM.plot,SUMM.plot$Chrom.Seg.Int==i-1)$Pos, 1)
    #and then for the subset of SUMM.plot rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Index=SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Pos+lastbase
  }
  #set ticks to be a list of existing ticks, plus the current Index
  #floor rounds it down to the nearest whole number
 # ticks=c(ticks, SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Index[floor(length(SUMM.plot[SUMM.plot$Chrom==i, ]$Index)/2)+1])
  
  ticks=c(ticks, SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Index[floor(length(SUMM.plot[SUMM.plot$Chrom.Seg.Int==i, ]$Index)/2)+1])
}
ticklim=c(min(SUMM.plot$Index),max(SUMM.plot$Index))

#troubleshoot
hist(SUMM.plot$Index, breaks=100)
#gaps are between: 3 to 3.1, 6 to 6.3
#the order of chr segments for 1, 2, 4, are all wrong
#chr segments 3 are in the right order
#I'm pretty sure that segments are moving based on chr 1, 2, 3... only

##get midpoint positions per chromosome
((max(SUMM.plot[ which(SUMM.plot$Chrom=='16'),]$Index) - min(SUMM.plot[ which(SUMM.plot$Chrom=='16'),]$Index))/2+min(SUMM.plot[ which(SUMM.plot$Chrom=='16'),]$Index))
#c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687)

#create a custom color scale
myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
#names(myColors) <- levels(HEM.plotdata$Chrom)
library(ggplot2)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#make plots 
#for poster figures, width=8, height=4
 jpeg("plots/MultiPlot/meta/FigR6_Summary_999Thresh_ManhattanPlot.jpg", width=7.5, height=4.4, units='in', res=600)
#SUMMtemp <- subset(SUMM.plot[SUMM.plot$Chrom==c(5,6),])
  ggplot(SUMM.plot, aes(x=Index, y=SUMM))+
    colScale+ #remove for rainbow plot
    theme_bw()+
#    scale_x_continuous(breaks = ticks)+
    geom_point(aes(color = factor(Chrom)))+
    labs(list(y="Number of Significant SNPs Across Plants", x="Chromosome position"))+
    #nrow=8
    theme(legend.position="none")+
    guides(col = guide_legend(nrow = 8, title="Chromosome"))+
    scale_x_continuous(name="Chromosome Position", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
  dev.off()
  
names(SUMM.plot)
TopSNPs <- SUMM.plot[which(SUMM.plot$SUMM > 1),]
write.csv(TopSNPs, "data/TopSNPs_ALLtomato.csv")
