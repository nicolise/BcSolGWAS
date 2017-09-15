#Nicole E Soltis
#091317
#plot of SNPs along gene of interest
#and now: haplotype plots

#14_C_HAP.FILE.R
#------------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

# first bring Chr16_A.qassoc.hap into excel or libreoffice calc. Do fixed-width import -- tab separated didn't work in R directly. Save as .csv
myHAP.1 <-  read.csv("data/genome/chr16_analysis/haps/chr16_A.qassoc.hap.csv", sep=",")
myHAP.2 <-  read.csv("data/genome/chr16_analysis/haps/chr16_A.qassoc_test.hap.csv", sep=",")

#columns for HAP.FILE:
#1. ASSOC (+, -) <- use STAT > 0 or STAT < 0
#ifelse(mySNP.FILE$LA1547 > 0, '+', '-')
#2. G.PVAL <- take average of I.PVAL for the whole LOCUS
#3. I.PVAL <- use P
#4. many columns, named by haplotype LOCUS
#Haplotypes are presented in a step-wise fashion with the major allele given as 1 and the minor allele as 2; haplotype variants for a set of SNPs should be grouped. SNP labels in HAP.FILE must be the same as in SNP.FILE, and only SNPs with corresponding haplotypes need to be included.

myHAP <- myHAP.2

#let's try to do this just for Win1
head(myHAP, 7)

#subset by LOCUS
myHAP.win1 <- myHAP[which(myHAP$LOCUS=="WIN1"),]
myHAP.win1
#find the major allele for HAPLOTYPE.1
tail(names(sort(table(myHAP.win1$HAPLOTYPE.1))), 1)
#now code a new variable with 1 if HAPLOTYPE.2 matches this, 2 if it does not
myHAP.win1$RECODE.1 <- ifelse(myHAP.win1$HAPLOTYPE.1 == tail(names(sort(table(myHAP.win1$HAPLOTYPE.1))), 1), 1, 2)
#then H.2
myHAP.win1$RECODE.2 <- ifelse(myHAP.win1$HAPLOTYPE.2 == tail(names(sort(table(myHAP.win1$HAPLOTYPE.2))), 1), 1, 2)
#then H.3
myHAP.win1$RECODE.3 <- ifelse(myHAP.win1$HAPLOTYPE.3 == tail(names(sort(table(myHAP.win1$HAPLOTYPE.3))), 1), 1, 2)

myHAP.append <- myHAP 
#make a nonsense vector with the right elements
myHAP.append$RECODE.1 <- c(rep(1:2, each=1468),1)
myHAP.append$RECODE.2 <- c(rep(1:2, each=1468),1)
myHAP.append$RECODE.3 <- c(rep(1:2, each=1468),1)
#now I'll just make it short
myHAP.append <- myHAP.append[1,]
#okay, now let's try to make this a loop
myWins<-unique(myHAP$LOCUS)
for (i in 1:length(myWins)){ 
  temp <- myHAP[myHAP$LOCUS==myWins[i],]
  #more things to do with temp
  temp$RECODE.1 <- ifelse(temp$HAPLOTYPE.1 == tail(names(sort(table(temp$HAPLOTYPE.1))), 1), 1, 2)
  temp$RECODE.2 <- ifelse(temp$HAPLOTYPE.2 == tail(names(sort(table(temp$HAPLOTYPE.2))), 1), 1, 2)
  temp$RECODE.3 <- ifelse(temp$HAPLOTYPE.3 == tail(names(sort(table(temp$HAPLOTYPE.3))), 1), 1, 2)
  myHAP.append <- rbind(myHAP.append, temp)
}
#and remove nonsense row 1
myHAP.append <- myHAP.append[-c(1),]

#kinda sorta make the dataframe we need
#columns for HAP.FILE:
#1. ASSOC (+, -) <- use STAT > 0 or STAT < 0
#2. G.PVAL <- take average of I.PVAL for the whole LOCUS
#3. I.PVAL <- use P
#4. many columns, named by haplotype LOCUS
#Haplotypes are presented in a step-wise fashion with the major allele given as 1 and the minor allele as 2; haplotype variants for a set of SNPs should be grouped. SNP labels in HAP.FILE must be the same as in SNP.FILE, and only SNPs with corresponding haplotypes need to be included.

myHAP.FILE <- myHAP.append[,c("STAT", "LOCUS","P", "HAPLOTYPE.1", "HAPLOTYPE.2", "HAPLOTYPE.3", "RECODE.1", "RECODE.2", "RECODE.3")]
myHAP.FILE$ASSOC <- ifelse(myHAP.FILE$STAT > 0, '+', '-')
myHAP.FILE$I.PVAL <- as.numeric(myHAP.FILE$P)

myHAP.FILE$G.PVAL <- mean(myHAP.FILE$I.PVAL[myHAP.FILE$LOCUS == "WIN1"])
blah <- by(myHAP.FILE$I.PVAL, myHAP.FILE$LOCUS, mean)

#new problem: how to fill in the SNPs in this staggered pattern
