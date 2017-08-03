#04_troubleshoot

#troubleshoot missingness
#--------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")

#this fails, not sure why
#myPed <- read.delim("GEMMA_files/02_csvPrep/fulldata/01_PLINK_55chrom/dpcharMAF20NA10.ped", sep=" ")
myPed <- read.csv("GEMMA_files/02_csvPrep/fulldata/PED_dpcharMAF20NA10.csv")
