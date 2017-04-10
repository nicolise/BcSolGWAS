#Nicole E Soltis
#040317
#plot of SNPs/ kb along genome

#---------------------------------------------
rm(list=ls())
library(plyr); library(tidyr); library(ggplot2)
setwd("~/Documents/GitRepos/BcSolGWAS/data/GWAS_files/")
setwd("~/Projects/BcSolGWAS/")
#read in MAF20 SNPs file
SolSNPs <- read.csv("data/GWAS_files/03_bigRRinput/binSNP_bigRR_MAF20hp_Sol.csv")
TestSNPs <- read.csv("02_csvPrep/hp_binaryMAF20.csv")
#make dataframe with 1 observation (sum of SNPs) per kb
SolSNPs <- TestSNPs
SolSNPs$Pos <- as.character(SolSNPs$POS)#ensure that position data is not in scientific notation
#SolSNPs <- SolSNPs[-c(1,3)]
names(SolSNPs)

#unique(SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome1$", replacement = "Chromosome1.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome2$", replacement = "Chromosome2.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome3$", replacement = "Chromosome3.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome4$", replacement = "Chromosome4.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome5$", replacement = "Chromosome5.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome6$", replacement = "Chromosome6.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome7$", replacement = "Chromosome7.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome8$", replacement = "Chromosome8.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome9$", replacement = "Chromosome9.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome10$", replacement = "Chromosome10.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome11$", replacement = "Chromosome11.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome12$", replacement = "Chromosome12.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome13$", replacement = "Chromosome13.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome14$", replacement = "Chromosome14.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome15$", replacement = "Chromosome15.0", SolSNPs$X.CHROM)
SolSNPs$X.CHROM <- gsub(pattern = "Chromosome16$", replacement = "Chromosome16.0", SolSNPs$X.CHROM)
unique(SolSNPs$X.CHROM)

#SolSNPs2 <- separate(SolSNPs, X.CHROM, into = c("Chromosome", "Segment") )
#double check
#unique(SolSNPs2$Chromosome)
#unique(SolSNPs2$Segment)

#new to edit
#Reformat Chromosomes and Positions
#SolSNPs2$Chromosome <- gsub("Chromosome", "", SolSNPs2$Chromosome)
SolSNPs2 <- SolSNPs
SolSNPs2$X.CHROM <- gsub("Chromosome", "", SolSNPs2$X.CHROM)
SolSNPs2$CHROM <- as.numeric(as.character(SolSNPs2$X.CHROM))
#SolSNPs2$Chromosome <- as.numeric(as.character(SolSNPs2$Chromosome))
#SolSNPs2$Segment <- as.numeric(as.character(SolSNPs2$Segment))
SolSNPs2$Pos <- as.numeric(as.character(SolSNPs2$Pos))

#sort dataframe rows in order of Chrom, then Pos
SolSNPs <- SolSNPs2
SolSNPs<-SolSNPs[with(SolSNPs, order(X.CHROM, Pos)), ]
#SolSNPs <- SolSNPs[with(SolSNPs, order(Chromosome, Segment, Pos)), ]

#now make segments line up consecutively
#SolSNPs$Chrom.Seg <- paste(SolSNPs$Chromosome, SolSNPs$Segment, sep=".")
#SolSNPs$Chrom.Seg <- as.numeric(SolSNPs$Chrom.Seg)

#let's try making the chrom.seg integers so that R isn't confused
#unique(SolSNPs$Chrom.Seg)
unique(SolSNPs$X.CHROM)
#SolSNPs$Chrom.Seg.F <- as.factor(SolSNPs$Chrom.Seg)
SolSNPs$Chrom.Seg.F <- as.factor(SolSNPs$X.CHROM)
unique(SolSNPs$Chrom.Seg.F)
recode.vars <- data.frame(OGvals = c("1.0",1.1,"2.0",2.1,2.2,2.3,2.4,2.5, "3.0", 3.1, 3.2, "4.0",  4.1, "5.0",  5.1, "6.0",  6.1, 6.2, 6.3, "7.0",  7.1, 7.2, "8.0",  8.1, 8.2, "9.0",  9.1, "10.0", 10.1, "11.0", 11.1, "12.0", 12.1, "13.0", 13.1, 13.2, "14.0", 14.1, 14.2, "15.0", 15.1, 15.2, 15.3, 15.4, "16.0", 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

SolSNPs$Chrom.Seg.Int <- recode.vars$newvals[match(SolSNPs$Chrom.Seg.F, recode.vars$OGvals)]
unique(SolSNPs$Chrom.Seg.Int)

SolSNPs <- SolSNPs[complete.cases(SolSNPs),]

#Make plotting variables
SolSNPs$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing

for (i in unique(SolSNPs$Chrom.Seg.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of SolSNPs rows with Chromosome 1, set Index variable for each row to equal Pos.
    SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Index=SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(SolSNPs,SolSNPs$Chrom.Seg.Int==i-1)$Pos, 1)
    #and then for the subset of SolSNPs rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Index=SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Pos+lastbase
  }
  #set ticks to be a list of existing ticks, plus the current Index
  #floor rounds it down to the nearest whole number
  # ticks=c(ticks, SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Index[floor(length(SolSNPs[SolSNPs$Chrom==i, ]$Index)/2)+1])
  
  ticks=c(ticks, SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Index[floor(length(SolSNPs[SolSNPs$Chrom.Seg.Int==i, ]$Index)/2)+1])
}
ticklim=c(min(SolSNPs$Index),max(SolSNPs$Index))

names(SolSNPs)
#All groups (4:6)
#keep only: SNPs over 99.9% Threshold
#assign(paste("HEM.", names(SolSNPs[6]), sep=""), subset(SolSNPs, abs(SolSNPs[6]) > get(paste("TH999_", names(SolSNPs[6]), sep="")), select=c(Chrom,Segment,Pos,Index,6)))
names(SolSNPs)
SolSNPs$Count <- rowSums(SolSNPs[5:101])
jpeg("plot/SNPcount.jpg", width=8, height=4, units='in', res=600)
hist(SolSNPs$Count)
dev.off()
SolSNPs$Freq <- (SolSNPs$Count)/95
jpeg("plot/SNPfreq.jpg", width=8, height=4, units='in', res=600)
hist(SolSNPs$Freq)
dev.off()
jpeg("plot/MAF20hist.jpg", width=8, height=4, units='in', res=600)
SolSNPs$MAF <- ifelse(SolSNPs$Freq > 0.5, 1-SolSNPs$Freq, SolSNPs$Freq) 
hist(SolSNPs$MAF)
dev.off()
  #dataframe$periodframe <- ifelse(dataframe$year > 1991,"post-1991", "pre-1991")

jpeg("plot/MAF20.jpg", width=8, height=4, units='in', res=600)
ggplot(SolSNPs, aes(x=Index, y=Freq))+
 #scale_y_continuous(limits = c(0,1.1))+
  theme_bw()+
  geom_point(aes(color = factor(X.CHROM)), alpha=0.2)
dev.off()
#   labs(list(y="SNP Frequency", x="Chromosome position", title=paste("Lesion Size on ", names(HEM.plotdata[9]))))+
#   guides(col = guide_legend(nrow = 8, title="Chromosome"))+
#   geom_hline(yintercept=get(paste("TH999_", names(HEM.plotdata[9]), sep=""))) +
#   geom_hline(yintercept=-0.0004520958)+
#   geom_text(aes(0,get(paste("TH999_", names(HEM.plotdata[9]), sep="")), label =
#                   ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")+
#   theme(legend.position="none")+
#   scale_x_continuous(name=element_blank(), breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
#   expand_limits(y=-0.001)
# dev.off()

SummDat <- ddply(ModDat, c("PlGenoNm", "Igeno", "Species", "ExpBlock"), summarise,
                 mLS   = mean(Scale.LS),
                 sdLS = sd(Scale.LS))

#use Index to plot x axis