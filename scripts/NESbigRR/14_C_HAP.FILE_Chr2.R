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
myHAP.1 <-  read.csv("data/genome/chr2_analysis/haps/chr2_A.qassoc.hap.csv", sep=",")
myHAP.2 <-  read.csv("data/genome/chr2_analysis/haps/chr2_A.qassoc_test.hap.csv", sep=",")

#split HAPLOTYPE into 3 columns
splitHaps <- as.data.frame(do.call(rbind, strsplit(as.character(myHAP.2$HAPLOTYPE), '\\.')))
names(splitHaps)[1] <- "HAPLOTYPE.1"
names(splitHaps)[2] <- "HAPLOTYPE.2"
names(splitHaps)[3] <- "HAPLOTYPE.3"
myHAP.2 <- cbind(myHAP.2, splitHaps)
myHAP.2 <- myHAP.2[,c(1,9,10,11,3:8)]

#columns for HAP.FILE:
#1. ASSOC (+, -) <- use STAT > 0 or STAT < 0
#ifelse(mySNP.FILE$LA1547 > 0, '+', '-')
#2. G.PVAL <- take average of I.PVAL for the whole LOCUS
#3. I.PVAL <- use P
#4. many columns, named by haplotype LOCUS
#Haplotypes are presented in a step-wise fashion with the major allele given as 1 and the minor allele as 2; haplotype variants for a set of SNPs should be grouped. SNP labels in HAP.FILE must be the same as in SNP.FILE, and only SNPs with corresponding haplotypes need to be included.

myHAP <- myHAP.2
#make a nonsense vector with the right elements (just needs to include 1 and 2 in any order)
myHAP$RECODE.1 <- c(rep(1:2, each=59))
myHAP$RECODE.2 <- c(rep(1:2, each=59))
myHAP$RECODE.3 <- c(rep(1:2, each=59))
#now I'll just make it short
myHAP.append <- myHAP 
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

#split NPS into 3 SNPs
splitSNPs <- as.data.frame(do.call(rbind, strsplit(as.character(myHAP.append$NPS), '\\|')))
names(splitSNPs)[1] <- "SNPS.1"
names(splitSNPs)[2] <- "SNPS.2"
names(splitSNPs)[3] <- "SNPS.3"
myHAP.append <- cbind(myHAP.append, splitSNPs)

#kinda sorta make the dataframe we need
#columns for HAP.FILE:
#1. ASSOC (+, -) <- use STAT > 0 or STAT < 0
#2. G.PVAL <- take average of I.PVAL for the whole LOCUS
#3. I.PVAL <- use P
#4. many columns, named by haplotype LOCUS
#Haplotypes are presented in a step-wise fashion with the major allele given as 1 and the minor allele as 2; haplotype variants for a set of SNPs should be grouped. SNP labels in HAP.FILE must be the same as in SNP.FILE, and only SNPs with corresponding haplotypes need to be included.
names(myHAP.append)[9] <- "P"
myHAP.FILE <- myHAP.append[,c("STAT", "LOCUS","P", "HAPLOTYPE.1", "HAPLOTYPE.2", "HAPLOTYPE.3", "RECODE.1", "RECODE.2", "RECODE.3", "SNPS.1","SNPS.2","SNPS.3")]
myHAP.FILE$ASSOC <- ifelse(myHAP.FILE$STAT > 0, '+', '-')
myHAP.FILE$I.PVAL <- as.numeric(as.character(myHAP.FILE$P))
#this doesn't work: only lists each mean once

mylist <- as.data.frame(tapply(myHAP.FILE$I.PVAL, myHAP.FILE$LOCUS, mean, na.rm=TRUE))
mylist$LOCUS <- rownames(mylist)
names(mylist)[1] <- "G.PVAL"
myHAP.FILE <- merge(myHAP.FILE, mylist, by="LOCUS")

#new problem: how to fill in the SNPs in this staggered pattern
unique(myHAP.FILE$SNPS.1)
head(myHAP.FILE)
HAP.FILE <- myHAP.FILE[,c("LOCUS","ASSOC", "G.PVAL", "I.PVAL", "SNPS.1", "SNPS.2", "SNPS.3", "RECODE.1","RECODE.2","RECODE.3")]
#sort by SNPs
HAP.FILE <- HAP.FILE[with(HAP.FILE, order(SNPS.1, SNPS.2, SNPS.3)), ]
#now abbreviate
HAP.FILE$SNPs <- paste(HAP.FILE$SNPS.1, HAP.FILE$SNPS.2, HAP.FILE$SNPS.3, sep=" ")
HAP.FILE$HAPs <- paste(HAP.FILE$RECODE.1, HAP.FILE$RECODE.2, HAP.FILE$RECODE.3, sep=" ")
HAP.FILE$myrows <- c(1:118)
#HAP.FILE <- HAP.FILE[,-c(4:9)]

HAP.FILE.s <- HAP.FILE[,c("SNPs", "HAPs", "myrows")]

#copy values of SNPs to be column names for HAPs
bleh <- reshape(HAP.FILE.s, idvar = "myrows", timevar = "SNPs", direction = "wide")

#and then maybe if it works, merge other variables back in
HAP.FILE.r <- HAP.FILE[,c("ASSOC", "G.PVAL", "I.PVAL", "myrows")]
HAP.FILE.r$G.PVAL <- as.numeric(HAP.FILE.r$G.PVAL)
HAP.FILE.fin <- merge(HAP.FILE.r, bleh, by = "myrows")
HAP.FILE.fin <- HAP.FILE.fin[,-c(1)]
write.csv(HAP.FILE.fin, "data/genome/chr2_analysis/haps/HAP.FILE.csv")
