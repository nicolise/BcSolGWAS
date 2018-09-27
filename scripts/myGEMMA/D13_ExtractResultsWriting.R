#Nicole E Soltis

#042018

#grabbing results for writing
#---------------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")

#---------------------------------------------------------------------------------------
#mean lesion size and SE for each tomato genotype?
setwd("~/Projects/BcSolGWAS")
ModDat <- read.csv("data/preGWAS/BcSolGWAS_PhenotypePlotData_ModDat.csv")
SummDat <- aggregate(ModDat[, "Scale.LS"], list(ModDat$PlGenoNm), mean)
SummDat$sd <- aggregate(ModDat[, "Scale.LS"], list(ModDat$PlGenoNm), sd)
SummDat$n <- aggregate(ModDat[, "Scale.LS"], list(ModDat$PlGenoNm), length)
SummDat$se <- SummDat$sd.x / (SummDat$n.x^0.5)
#---------------------------------------------------------------------------------------
#genome info for B05.10 GEMMA?
#237,878 SNPs at MAF 0.20 or greater and  less than 10% missing SNP calls

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
#How many SNPs significant across multiple traits in GEMMA?

#99% Thr
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut99Thr_kmat1.csv")

#99.9% Thr
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut999Thr_kmat1.csv")

table(pheno.bin$SUMM)

#and for domestication?
#99.9% Thr
myGEMMA <- read.csv("data/GEMMA_files/D_08_results/GEMMA_peaksDWS_kmat1_999thr.csv")

#99% Thr
myGEMMA <- read.csv("data/GEMMA_files/D_08_results/GEMMA_peaksDWS_kmat1_99thr.csv")

#domestication SNP-level Venn
table(myGEMMA$TotTraits)
#--------------------------------------------------------------------------------------
#How many genes significant across multiple traits in GEMMA?

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

#domestication gene-level Venn for GEMMA here
table(domest.genes$TotTraits)
#-----------------------------------------------------------------------------------
how many genes significant across multiple traits in bigRR?
setwd("~/Projects/BcSolGWAS")
HOplant2 <- read.csv("data/GWAS_files/05_annotation/window2kb/12plants_HO_genesTOANNOT.csv")
FullGenes <- read.csv("data/GWAS_files/05_annotation/OverlapCount_plant12topgenes_99thr.csv")
table(FullGenes$TotTraits)
#--------------------------------------------------------------------------------------
#gene overlap from GEMMA to bigRR?
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")
#99.9% Thr

DoGenOverlap <- read.csv("D_08_results/AllDOfuncs_byGene_999thr.csv")
DoGenSumm <- DoGenOverlap[,c("T4vanKan.BROAD", "BcinB0510gene","TotTraits.b","TotTraits.G","PFAM_DESCRIPTION")]
DoGenSumm <- DoGenSumm[!duplicated(DoGenSumm$T4vanKan.BROAD),]
write.csv(DoGenSumm, "D_09_tables/Overlap_Domest_tableS3.csv")

HoGenOverlap <- read.csv("D_08_results/AllHOannots_byGene_999thr.csv")
Gen12Overlap <- read.csv("D_08_results/All12annots_byGene_999thr.csv") #this is for only >1 pheno genes by each method
Gen12Summ <- Gen12Overlap[,c("T4vanKan.BROAD", "BcinB0510gene","PFAM_DESCRIPTION","TotPhenos.G","TotPhenos.b","TotPhenos.O")]
Gen12Summ <- Gen12Summ[!duplicated(Gen12Summ$T4vanKan.BROAD),]
write.csv(Gen12Summ, "D_09_tables/Overlap_12pGenes_tableS3.csv")
#-----------------------------------------------------------------------------------
#function overlap from GEMMA to bigRR?
#summarize by functional annotation

setwd("~/Projects/BcSolGWAS/data/GEMMA_files")
#99.9% Thr
pheno12.overlap.Ant <- read.csv("D_08_results/HiOverlap12p_FuncOverrep_999thr_tops.csv")

Domest.overlap.Ant <- read.csv("D_08_results/Domestication_funcOverrep_999thr_tops.csv")