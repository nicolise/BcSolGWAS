#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS")
############################################################################
###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/Sl_LesionSize.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/Sl_LesionSize.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]
TH99 <- HEM.thresh[3,]
TH99_LA0410 <- as.numeric(TH99[2])
TH99_LA0480 <- as.numeric(TH99[3])
TH99_LA1547 <- as.numeric(TH99[4])
TH99_LA1589 <- as.numeric(TH99[5])
TH99_LA1684 <- as.numeric(TH99[6])
TH99_LA2093 <- as.numeric(TH99[7])
TH99_LA2176 <- as.numeric(TH99[8])
TH99_LA2706 <- as.numeric(TH99[9])
TH99_LA3008 <- as.numeric(TH99[10])
TH99_LA3475 <- as.numeric(TH99[11])
TH99_LA4345 <- as.numeric(TH99[12])
TH99_LA4355 <- as.numeric(TH99[13])

TH999 <- HEM.thresh[4,]
TH999_LA0410 <- as.numeric(TH999[2])
TH999_LA0480 <- as.numeric(TH999[3])
TH999_LA1547 <- as.numeric(TH999[4])
TH999_LA1589 <- as.numeric(TH999[5])
TH999_LA1684 <- as.numeric(TH999[6])
TH999_LA2093 <- as.numeric(TH999[7])
TH999_LA2176 <- as.numeric(TH999[8])
TH999_LA2706 <- as.numeric(TH999[9])
TH999_LA3008 <- as.numeric(TH999[10])
TH999_LA3475 <- as.numeric(TH999[11])
TH999_LA4345 <- as.numeric(TH999[12])
TH999_LA4355 <- as.numeric(TH999[13])

TH95 <- HEM.thresh[1,]
TH95_LA0410 <- as.numeric(TH95[2])
TH95_LA0480 <- as.numeric(TH95[3])
TH95_LA1547 <- as.numeric(TH95[4])
TH95_LA1589 <- as.numeric(TH95[5])
TH95_LA1684 <- as.numeric(TH95[6])
TH95_LA2093 <- as.numeric(TH95[7])
TH95_LA2176 <- as.numeric(TH95[8])
TH95_LA2706 <- as.numeric(TH95[9])
TH95_LA3008 <- as.numeric(TH95[10])
TH95_LA3475 <- as.numeric(TH95[11])
TH95_LA4345 <- as.numeric(TH95[12])
TH95_LA4355 <- as.numeric(TH95[13])


# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata2 <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Pos)), ]
range(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Pos)
min(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Pos)
max(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Pos)
range(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Pos)
min(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Pos)
max(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Pos)
#range = min to max
#head and tail not equivalent

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequencial		-- accurate position indexing
##RFF code isn't working... replaced "Pos" with "pos"
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#make plots for each phenotype
#try a log-scaled plot
HEM.plotdata$logLA0410 <- log(HEM.plotdata$LA0410 + 1)
l.TH99_LA0410 <- log(TH99_LA0410 + 1)
l.TH95_LA0410 <- log(TH95_LA0410 + 1)
#plot it!

jpeg("plots/Sl_LesionSize_LA0410.low.ManhattanPlot.jpg")
qplot(Index,LA0410, data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA0410", colour=factor(Chrom)) +
  scale_y_continuous(limits=c(-0.002,0.002))+
 geom_hline(yintercept=TH999_LA0410) +
  geom_text(aes(0,TH999_LA0410, label = ".999 Threshold", vjust = 1.5, hjust = .05), col = "black") +
 geom_hline(yintercept=TH99_LA0410) +
  geom_text(aes(0,TH99_LA0410, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black")
dev.off()

jpeg("plots/Sl_LesionSize_LA1547.low.ManhattanPlot.jpg")
qplot(Index,LA1547, data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA1547", colour=factor(Chrom)) +
  scale_y_continuous(limits=c(-0.000025,0.000025))+
  geom_hline(yintercept=TH95_LA1547) +
  geom_text(aes(0,TH95_LA1547, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

#how many SNPs are above a certain threshhold?
sum(HEM.plotdata$LA0410 >= TH999_LA0410) #spits out number
#491

#see top 20 SNPs for each
HEM.plotdata$LA0410
require(data.table)
d <- data.table(HEM.plotdata, key="LA0410")
dprint <- d[, head(.SD, 20), by=LA0410]

library(reshape2)
attach(topSNPs)
wideSNPs <- (dcast(topSNPs, Num~Geno, value.var=c('Chr.Pos'), fun=mean))
pairs(wideSNPs)

#view list of Indices
m <- ggplot(HEM.plotdata, aes(x=Index, fill=as.factor(Chrom)))
m + geom_bar() 

library(plyr)
maxLA0410 <- HEM.plotdata[order(HEM.plotdata$LA0410,decreasing=T)[1:20],c("Index","LA0410")]
maxLA0410 <- rename(maxLA0410, c("LA0410" = "Effects"))
maxLA0410$Geno <- "LA0410"
maxLA0480 <- HEM.plotdata[order(HEM.plotdata$LA0480,decreasing=T)[1:20],c("Index","LA0480")]
maxLA0480 <- rename(maxLA0480, c("LA0480" = "Effects"))
maxLA0480$Geno <- "LA0480"
maxLA1547 <- HEM.plotdata[order(HEM.plotdata$LA1547,decreasing=T)[1:20],c("Index","LA1547")]
maxLA1547 <- rename(maxLA1547, c("LA1547" = "Effects"))
maxLA1547$Geno <- "LA1547"
maxLA1589 <- HEM.plotdata[order(HEM.plotdata$LA1589,decreasing=T)[1:20],c("Index","LA1589")]
maxLA1589 <- rename(maxLA1589, c("LA1589" = "Effects"))
maxLA1589$Geno <- "LA1589"
maxLA1684 <- HEM.plotdata[order(HEM.plotdata$LA1684,decreasing=T)[1:20],c("Index","LA1684")]
maxLA1684 <- rename(maxLA1684, c("LA1684" = "Effects"))
maxLA1684$Geno <- "LA1684"
maxLA2093 <- HEM.plotdata[order(HEM.plotdata$LA2093,decreasing=T)[1:20],c("Index","LA2093")]
maxLA2093 <- rename(maxLA2093, c("LA2093" = "Effects"))
maxLA2093$Geno <- "LA2093"
maxLA2176 <- HEM.plotdata[order(HEM.plotdata$LA2176,decreasing=T)[1:20],c("Index","LA2176")]
maxLA2176 <- rename(maxLA2176, c("LA2176" = "Effects"))
maxLA2176$Geno <- "LA2176"
maxLA2706 <- HEM.plotdata[order(HEM.plotdata$LA2706,decreasing=T)[1:20],c("Index","LA2706")]
maxLA2706 <- rename(maxLA2706, c("LA2706" = "Effects"))
maxLA2706$Geno <- "LA2706"
maxLA3008 <- HEM.plotdata[order(HEM.plotdata$LA3008,decreasing=T)[1:20],c("Index","LA3008")]
maxLA3008 <- rename(maxLA3008, c("LA3008" = "Effects"))
maxLA3008$Geno <- "LA3008"
maxLA3475 <- HEM.plotdata[order(HEM.plotdata$LA3475,decreasing=T)[1:20],c("Index","LA3475")]
maxLA3475<- rename(maxLA3475, c("LA3475" = "Effects"))
maxLA3475$Geno <- "LA3475"
maxLA4345 <- HEM.plotdata[order(HEM.plotdata$LA4345,decreasing=T)[1:20],c("Index","LA4345")]
maxLA4345 <- rename(maxLA4345, c("LA4345" = "Effects"))
maxLA4345$Geno <- "LA4345"
maxLA4355 <- HEM.plotdata[order(HEM.plotdata$LA4355,decreasing=T)[1:20],c("Index","LA4355")]
maxLA4355 <- rename(maxLA4355, c("LA4355" = "Effects"))
maxLA4355$Geno <- "LA4355"
topSNPs <- rbind(maxLA0410,maxLA0480,maxLA1547,maxLA1589,maxLA1684,maxLA2093,maxLA2176,maxLA2706,maxLA3008,maxLA3475,maxLA4345,maxLA4355)
topSNPs$Num <- rep(c(1:20),12)

library(reshape2)
attach(topSNPs)
wideSNPs <- (dcast(topSNPs, Num~Geno, value.var=c('Index'), fun=mean))
pairs(wideSNPs)
#still not the same SNPs across genomes

table(unique(topSNPs$Index))
hist(topSNPs$Index, breaks=200)

#try running to make the plot with Effect Size scaled -- see lower peaks

library(plyr)
#check: are top 20 SNPs/ plant geno driven by the organic subset of isolates?
#1. list top 20 SNP positions per geno by Chrom, Pos - from  HEM.plotdata
#2. go back to binSNP_bigRR_MAF20hp.csv and compare organics vs. random subset of 5 isolates vs. B05.10 at these SNPs

maxLA0410 <- HEM.plotdata[order(HEM.plotdata$LA0410,decreasing=T)[1:20],c("Chrom","Pos","Index","LA0410")]
maxLA0410 <- rename(maxLA0410, c("LA0410" = "Effects"))
maxLA0410$Geno <- "LA0410"
maxLA0480 <- HEM.plotdata[order(HEM.plotdata$LA0480,decreasing=T)[1:20],c("Chrom","Pos","Index","LA0480")]
maxLA0480 <- rename(maxLA0480, c("LA0480" = "Effects"))
maxLA0480$Geno <- "LA0480"
maxLA1547 <- HEM.plotdata[order(HEM.plotdata$LA1547,decreasing=T)[1:20],c("Chrom","Pos","Index","LA1547")]
maxLA1547 <- rename(maxLA1547, c("LA1547" = "Effects"))
maxLA1547$Geno <- "LA1547"
maxLA1589 <- HEM.plotdata[order(HEM.plotdata$LA1589,decreasing=T)[1:20],c("Chrom","Pos","Index","LA1589")]
maxLA1589 <- rename(maxLA1589, c("LA1589" = "Effects"))
maxLA1589$Geno <- "LA1589"
maxLA1684 <- HEM.plotdata[order(HEM.plotdata$LA1684,decreasing=T)[1:20],c("Chrom","Pos","Index","LA1684")]
maxLA1684 <- rename(maxLA1684, c("LA1684" = "Effects"))
maxLA1684$Geno <- "LA1684"
maxLA2093 <- HEM.plotdata[order(HEM.plotdata$LA2093,decreasing=T)[1:20],c("Chrom","Pos", "Index", "LA2093")]
maxLA2093 <- rename(maxLA2093, c("LA2093" = "Effects"))
maxLA2093$Geno <- "LA2093"
maxLA2176 <- HEM.plotdata[order(HEM.plotdata$LA2176,decreasing=T)[1:20],c("Chrom","Pos","Index","LA2176")]
maxLA2176 <- rename(maxLA2176, c("LA2176" = "Effects"))
maxLA2176$Geno <- "LA2176"
maxLA2706 <- HEM.plotdata[order(HEM.plotdata$LA2706,decreasing=T)[1:20],c("Chrom","Pos","Index","LA2706")]
maxLA2706 <- rename(maxLA2706, c("LA2706" = "Effects"))
maxLA2706$Geno <- "LA2706"
maxLA3008 <- HEM.plotdata[order(HEM.plotdata$LA3008,decreasing=T)[1:20],c("Chrom","Pos","Index","LA3008")]
maxLA3008 <- rename(maxLA3008, c("LA3008" = "Effects"))
maxLA3008$Geno <- "LA3008"
maxLA3475 <- HEM.plotdata[order(HEM.plotdata$LA3475,decreasing=T)[1:20],c("Chrom","Pos","Index","LA3475")]
maxLA3475<- rename(maxLA3475, c("LA3475" = "Effects"))
maxLA3475$Geno <- "LA3475"
maxLA4345 <- HEM.plotdata[order(HEM.plotdata$LA4345,decreasing=T)[1:20],c("Chrom","Pos","Index","LA4345")]
maxLA4345 <- rename(maxLA4345, c("LA4345" = "Effects"))
maxLA4345$Geno <- "LA4345"
maxLA4355 <- HEM.plotdata[order(HEM.plotdata$LA4355,decreasing=T)[1:20],c("Chrom","Pos","Index","LA4355")]
maxLA4355 <- rename(maxLA4355, c("LA4355" = "Effects"))
maxLA4355$Geno <- "LA4355"
topSNPs <- rbind(maxLA0410,maxLA0480,maxLA1547,maxLA1589,maxLA1684,maxLA2093,maxLA2176,maxLA2706,maxLA3008,maxLA3475,maxLA4345,maxLA4355)
topSNPs$Num <- rep(c(1:20),12)

#First need to modify Chrom.Pos without the Chrom.SubChrom.Pos
SNPs <- read.csv("data/GWAS_files/03_bigRRinput/binSNP_bigRR_MAF20hp.csv", row.names = 1)
SNPs$X.CHROM <- gsub(pattern = "Chromosome1\\.[0-9]", replacement = "Chromosome1", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome2\\.[0-9]", replacement = "Chromosome2", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome3\\.[0-9]", replacement = "Chromosome3", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome4\\.[0-9]", replacement = "Chromosome4", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome5\\.[0-9]", replacement = "Chromosome5", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome6\\.[0-9]", replacement = "Chromosome6", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome7\\.[0-9]", replacement = "Chromosome7", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome8\\.[0-9]", replacement = "Chromosome8", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome9\\.[0-9]", replacement = "Chromosome9", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome10\\.[0-9]", replacement = "Chromosome10", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome11\\.[0-9]", replacement = "Chromosome11", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome12\\.[0-9]", replacement = "Chromosome12", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome13\\.[0-9]", replacement = "Chromosome13", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome14\\.[0-9]", replacement = "Chromosome14", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome15\\.[0-9]", replacement = "Chromosome15", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome16\\.[0-9]", replacement = "Chromosome16", SNPs$X.CHROM)
SNPs$X.CHROM <- gsub(pattern = "Chromosome16[0-9]", replacement = "Chromosome16", SNPs$X.CHROM)

#Reformat Chromosomes and Positions
SNPs$Chrom <- gsub("Chromosome", "", SNPs$X.CHROM)
SNPs$Chrom <- as.numeric(as.character(SNPs$Chrom))
SNPs$Pos <- as.numeric(as.character(SNPs$POS))

#summary table for Chi-square
SNPchi <- SNPs[,c(1:3,95,96)]
myG1 <- names(SNPs) %in% c("X2.04.04","X2.04.21","X2.04.18","X2.04.14","X2.04.20","X2.04.17","X2.04.03","X2.04.09","X2.04.11","X2.04.08","X2.04.12")
myomit <- names(SNPs) %in% c("X.CHROM","POS","REF","Chrom","Pos", "X2.04.04","X2.04.21","X2.04.18","X2.04.14","X2.04.20","X2.04.17","X2.04.03","X2.04.09","X2.04.11","X2.04.08","X2.04.12")

#make a chi-square table
#categorize by isolate
SNPchi$sumG1 <- rowSums(SNPs[,myG1])
SNPchi$sumAll <- rowSums(SNPs[,!myomit])
#categorize by locus
SNPchi$Chr.Pos <- paste(SNPchi$Chrom, SNPchi$Pos, sep=".")
SNPchi2 <- SNPchi[SNPchi$Chr.Pos %in% topSNPs$Chr.Pos,]
`%ni%` = Negate(`%in%`) 
SNPchi3 <- SNPchi[SNPchi$Chr.Pos %ni% topSNPs$Chr.Pos,]
sum(SNPchi2$sumG1)
sum(SNPchi2$sumAll)
sum(SNPchi3$sumG1)
sum(SNPchi3$sumAll)

myChiSq <- read.table(text = '
otherIsos OGIsos
7261 1452
12770285 2244380
', header = TRUE, stringsAsFactors = FALSE)
row.names(myChiSq) <- c("highouts","others")
#asks: is SNP outlier independent of isolate group?
chisq.test(myChiSq)


#Subset the SNPs of interest
#keeping only certain genotypes
SNPsub <- SNPs[, sample(ncol(SNPs[,4:94]), 12)]
SNPsub <- cbind(SNPs[,c(95,96,3)], SNPsub)
SNPsubC1 <- cbind(SNPsub, SNPs[,c("X2.04.04","X2.04.21","X2.04.18","X2.04.14","X2.04.20","X2.04.17","X2.04.03","X2.04.09","X2.04.11","X2.04.08","X2.04.12")])
SNPsubC2 <- cbind(SNPsub, SNPs[,c("X1.02.01","X1.02.02","X1.02.06","X1.02.20","X1.02.03")])
SNPsubC3 <- cbind(SNPsub, SNPs[,c("X1.01.06","X1.01.01","X1.01.02","X1.01.03","X1.01.04")])
SNPsubC4 <- cbind(SNPsub, SNPs[,c("X1.04.03","X1.02.15","X1.04.04","X1.04.05")])
SNPsubC5 <- cbind(SNPsub, SNPs[,c("Peachy","X1.04.20","X1.04.25","X1.04.21","X1.02.16")])

#only keep SNPsubs that match Chr.Pos for topSNPs
SNPsubC1$Chr.Pos <- paste(SNPsubC1$Chrom, SNPsubC5$Pos, sep=".")
topSNPs$Chr.Pos <- paste(topSNPs$Chrom, topSNPs$Pos, sep = ".")
SNPsubMatch <- SNPsubC1
#why am I losing so many? Should have 240 matches, only get 198
SNPsubMatch <- SNPsubMatch[SNPsubMatch$Chr.Pos %in% topSNPs$Chr.Pos,]
topSNPMatch <- topSNPs[topSNPs$Chr.Pos %in% SNPsub$Chr.Pos,] #this has all 240 so that's good...
SNPsubMatch$Chr.Pos <- as.numeric(SNPsubMatch$Chr.Pos)
topSNPMatch$Chr.Pos <- as.numeric(topSNPMatch$Chr.Pos)

#add Geno column to SNPsubMatch based on matching value in Index with topSNPs
topSNPmini <- topSNPs[,c("Chr.Pos","Geno")]
SNPsubMatch <- merge(SNPsubMatch, topSNPmini, by="Chr.Pos")


write.csv(SNPsubMatch, "data/GWAS_files/04_bigRRoutput/High20subsetC1.csv")
write.csv(topSNPMatch, "data/GWAS_files/04_bigRRoutput/High20SNPs.csv")

#now compare organics vs. random subset of 5 isolates vs. B05.10 at these SNPs
#add a column of sums per group per locus
SNPsubMatch$sumG1 <- rowSums(SNPsubMatch[,c("X2.04.04","X2.04.21","X2.04.18","X2.04.14","X2.04.20","X2.04.17","X2.04.03","X2.04.09","X2.04.11","X2.04.08","X2.04.12")]) 
SNPsubMatch$sumG2 <- rowSums(SNPsubMatch[,c("Ausubel","Gallo1", "KatieTomato","X1.04.20","X1.03.18")])
