#Nicole E Soltis
#091317
#plot of SNPs along gene of interest
#and now: haplotype plots

#14_C_FAM_importphenos.R
#-------------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

myFAM <-  read.delim("data/genome/chr16_analysis/plink/dummyPheno/myCHR16_A.fam",sep=" ",header=FALSE)

#check same phenos as 04_runbigRR_indplants.R
myPhenos <- read.csv("data/GWAS_files/03_bigRRinput/NewModel0711/Sl_Pheno_bigRR_trueMAF20_10NA.csv")

SNPnames <- read.csv("data/GWAS_files/02_csvPrep/Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,4)]

names(SNPnames)[1] <- "Igeno"
names(SNPnames)[2] <- "V2"

myPhenos <- merge(myPhenos,SNPnames, by="Igeno", all=TRUE)
#remove empty isolates
myPhenos <- myPhenos[c(1:94),]
addPhenos <- myPhenos[,c("V2","LA1547")]

#match files on myFAM$V2 (original )
myFAM.2 <- merge(myFAM,addPhenos, by="V2")
#fill in missing phenos with average
mean(myFAM.2$LA1547)
myFAM.2 <- rbind(myFAM.2, c("B236", "FAM1", 0, 0, 1, 1, 0.5044))
myFAM.2 <- rbind(myFAM.2, c("B4", "FAM1", 0, 0, 1, 1, 0.5044))
myFAM.2 <- rbind(myFAM.2, c("B78", "FAM1", 0, 0, 1, 1, 0.5044))

match(myFAM$V2,myFAM.2$V2)


#write FAM with real phenotype
write.table(myFAM, file="data/genome/chr16_analysis/plink/myCHR16_A.fam", quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
