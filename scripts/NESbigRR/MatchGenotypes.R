#Nicole E Soltis
#match Bc isolates from genotype data to Bc isolates in phenotype data

#-----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
#for laptop setwd("~/Projects/BcSolGWAS/data/genome")
SNPs <- read.csv("hp_charMAF20.csv", row.names = 1)
#SNPs <- read.csv("miniSNP_practice.csv") 
#SNPsDF <- SNPs
#SNPsDF <- SNPsDF[c(1:2),]
#write.csv(SNPsDF, "SNPgenos.csv")
SNPs_rename <- SNPs

SNPnames <- read.csv("Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]
Phenos <- read.csv("LSMforbigRR_est_Xnames.csv")

#change names from genotype file to match phenotype file
#File SNPs_rename has columns of isolate genotypes that I want to rename
#File SNPnames is a list indexing original names to match (SNPname), and actual name to replace it (GenoRename)
names( SNPs_rename ) <- SNPnames[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

## now only keep genotypes and phenotypes that match

#save practice file
#miniSNPs <- SNPs_rename[c(1:3),]
#write.csv(miniSNPs, "miniSNP_practice.csv")

#only keep phenotype rows that match SNP names
SNPMatch <- as.data.frame(names(SNPs))
PhenoMatch <- Phenos
PhenoMatch <- PhenoMatch[PhenoMatch$Igeno %in% SNPMatch$"names(SNPs)", ]

#only keep SNP rows that match phenotype names
PhenoMt <- t(as.data.frame(Phenos[,1]))
PhenoMt <- as.data.frame(t(Phenos[,1]))

names( SNPs_rename ) <- PhenoMatch[ match( names( SNPs_rename ) , SNPnames[ , 'SNPname' ] ) ,'GenoRename' ] 

#extra things
#miniSNPs <- as.data.frame(t(miniSNPs))
#miniPhenos <- subset(Phenos, Igeno %in% SNPs_rename[0,])
