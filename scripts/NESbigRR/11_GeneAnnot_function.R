#11_GeneAnnot_function
#Nicole E Soltis
#052217

#Input: SNP data/GWAS_files/05_annotation/Domest_topSNPs_genestoAnnot.csv 
#and Gene data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv
#from script 10_GeneAnnot_10NA_venns_figR8.R
#Output:
#Figures: None
#---------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
library(dplyr)
SNPsForAnnot <- read.csv("data/GWAS_files/05_annotation/SlBcDomest_10NA_topSNPs_genestoAnnot.csv")
GenesForAnnot <- read.csv("data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv")
FuncAnnot <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_use.csv")

head(FuncAnnot)

Funcs <- FuncAnnot[,c(2,3,4)]
Funcs <- Funcs[!duplicated(Funcs[,1]),]

head(SNPsForAnnot)
SNPsForAnnot$GENE <- SNPsForAnnot$geneID

#this is the full list including all SNPs
AnnotSNPs <- merge(SNPsForAnnot, Funcs, by="GENE")

#now within each gene, merge back onto AnnotGenes
colnames(Funcs)[1] <- "geneID"
GenesForAnnot <- GenesForAnnot[,-c(1)]
DoGenAnt2 <- merge(GenesForAnnot, Funcs, by="geneID")
#keep only the first mention of each gene
DoGenAnt2 <- DoGenAnt2[!duplicated(DoGenAnt2[,1]),]

#summarize annotations
DoGenAnt2 %>%
  group_by(PFAM_DESCRIPTION) %>%
  summarize(count_genes = count(geneID, na.rm = TRUE))

mytable <- as.data.frame(table(droplevels(DoGenAnt2$PFAM_DESCRIPTION)))
write.csv(mytable, "data/GWAS_files/05_annotation/Domest_NA10_FunctionsSummary.csv")
