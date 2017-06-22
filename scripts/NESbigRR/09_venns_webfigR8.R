#Nicole E Soltis 
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input File: Sl_DomesticationLS_MAF20.HEM.PlotFormat.csv and .Thresh.csv
#Output File: results/Domestication_TopSNPs_SegLong.csv, results/Domestication_TopSNPs_SegWide.csv
#this takes gene annotation and goes into venn diagrams

DomestSNP <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1kSNPs_SegWide_trueMAF20_10NA.csv")

AnnotSNP <- read.csv(data/)
