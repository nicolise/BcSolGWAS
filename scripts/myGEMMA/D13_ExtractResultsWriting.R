#Nicole E Soltis

#042018

#grabbing results for writing
#---------------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")

#---------------------------------------------------------------------------------------
#How many SNPs analyzed per phenotype in GEMMA?
#check log files in D_04_randphenos

#---------------------------------------------------------------------------------------

#How many SNPs significant per plant pheno in GEMMA?
#for 99% Thr
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut99Thr_kmat1.csv")
#for 99.9% Thr
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut999Thr_kmat1.csv")

table(pheno.bin$X1_LA0410_betabin)
table(pheno.bin$X2_LA0480_betabin)
table(pheno.bin$X3_LA1547_betabin)
table(pheno.bin$X4_LA1589_betabin)
table(pheno.bin$X5_LA1684_betabin)
table(pheno.bin$X6_LA2093_betabin)
table(pheno.bin$X7_LA2176_betabin)
table(pheno.bin$X8_LA2706_betabin)
table(pheno.bin$X9_LA3008_betabin)
table(pheno.bin$X10_LA3475_betabin)
table(pheno.bin$X11_LA4345_betabin)
table(pheno.bin$X12_LA4355_betabin)

#---------------------------------------------------------------------------------------
#How many SNPs significant across multiple traits?

#99% Thr
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut99Thr_kmat1.csv")

#99.9% Thr
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut999Thr_kmat1.csv")

table(pheno.bin$SUMM)
#--------------------------------------------------------------------------------------
#How many genes significant across multiple traits?

#99% Thr

DWS.G <- read.csv("D_08_results/domestplants_Genes4Func_99thr.csv")
#this subsets the "high overlap" SNPs: SNPs that are sig in > 6 phenotypes
pheno12.G <- read.csv("D_08_results/hi12plants_Genes4Func_99thr.csv")
#this is just top 1000 SNPs per phenotype
pheno12HO.G <- read.csv("D_08_results/hi12HOplants_Genes4Func_99thr.csv")

#99.9% Thr

DWS.G <- read.csv("D_08_results/domestplants_Genes4Func_999thr.csv")
pheno12.G <- read.csv("D_08_results/hi12plants_Genes4Func_999thr.csv")
pheno12HO.G <- read.csv("D_08_results/hi12HOplants_Genes4Func_999thr_bad.csv")

#only keep one record per gene
plant12.genes<- pheno12.G[!duplicated(pheno12.G$geneID),]
plantHO.genes <- pheno12HO.G[!duplicated(pheno12HO.G$geneID),]
domest.genes <- DWS.G[!duplicated(DWS.G$geneID),]

table(plant12.genes$TotPhenos)
table(plantHO.genes$TotPhenos)
#-----------------------------------------------------------------------------------
#gene overlap from GEMMA to bigRR?
pheno12.overlap.Ant <- read.csv("D_08_results/Annotated_pheno12_overlap_999thr.csv")

pheno12.overlap.Ant <- read.csv("D_08_results/Annotated_pheno12_overlap_99thr.csv")

#-----------------------------------------------------------------------------------
#function overlap from GEMMA to bigRR?
#summarize by functional annotation

annotsSUMM <- as.data.frame(table(pheno12.overlap.Ant$PFAM_DESCRIPTION))