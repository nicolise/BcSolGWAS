#03_FAM_addPhenotypes
#Nicole E Soltis
#072717
#add actual phenotype column to .fam file
#--------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")

myFam <- read.delim("GEMMA_files/02_csvPrep/fulldata/PLINK_55chrom/dpcharMAF20NA10_out.fam", sep=" ")
myPhenos <- read.csv("GEMMA_files/02_csvPrep/Sl_Pheno_bigRR_trueMAF20_10NA.csv")

#first: add column with isolate nicknames from Suzi to phenos
#then: merge myPhenos file to myFam

#add all 12 phenotypes. Each time I run GEMMA, specify which column (col 6 = 1, col 7 = 2, etc.)

myFam[,6] <- 

write.delim(myPED, "GEMMA_files/02_csvPrep/fulldata/PED_dpcharMAF20NA10_55chrom.tsv", quote = FALSE, col.names = F, row.names = FALSE, sep = "\t")
