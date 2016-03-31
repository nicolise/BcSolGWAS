#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")

SNPs <- read.csv("hp_charMAF20.csv", row.names = 1)
#SNPsDF <- SNPs
#SNPsDF <- SNPsDF[c(1:2),]
#write.csv(SNPsDF, "SNPgenos.csv")
SNPs_rename <- SNPs

SNPnames <- read.csv("Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]
Phenos <- read.csv("LSMforbigRR_est.csv")

#change names from genotype file to match phenotype file
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

#now only keep genotypes and phenotypes that match
miniSNPs <- SNPs_rename[c(1:3),]
#write.csv(miniSNPs, "miniSNP_practice.csv")
miniSNPs <- as.data.frame(t(miniSNPs))

miniPhenos <- subset(Phenos, Igeno %in% SNPs_rename[0,])
