#Nicole E Soltis
#B08_GEMMA_bigRR_overlaps

#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")

DWS.G <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/domestplants_Genes4Func.csv")
pheno12.G <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/hi12plants_Genes4Func.csv")
pheno12HO.G <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/hi12HOplants_Genes4Func.csv")

DWS.b <- read.csv("data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv")
pheno12.b <- read.csv("data/GWAS_files/05_annotation/window2kb/hi12plants_genesTOANNOT.csv")
pheno12HO.b <- read.csv("data/GWAS_files/05_annotation/window2kb/12plants_HO_genesTOANNOT.csv")

names(DWS.b)
unique(DWS.b$geneID)
names(DWS.G)
unique(DWS.G$geneID)
names(pheno12.b)
names(pheno12HO.b)
names(pheno12.G)
names(pheno12HO.G)
#extract just gene ids to get BCINs from BotPortal
gene.D.G <- unique(DWS.G$geneID)
gene.D.G <- gsub("gene:", "", gene.D.G)
gene.D.G <- as.data.frame(gene.D.G)
names(gene.D.G)[1] <- "geneID"

gene.D.b <- as.data.frame(unique(DWS.b$geneID))
names(gene.D.b)[1] <- "geneID"

gene.12.G <- unique(pheno12.G$geneID)
gene.12.G <- gsub("gene:", "", gene.12.G)
gene.12.G <- as.data.frame(gene.12.G)
names(gene.12.G)[1] <- "geneID"

gene.12.b <- as.data.frame(unique(pheno12.b$geneID))
names(gene.12.b)[1] <- "geneID"

gene.HO.G <- unique(pheno12HO.G$geneID)
gene.HO.G <- gsub("gene:", "", gene.HO.G)
gene.HO.G <- as.data.frame(gene.HO.G)
names(gene.HO.G)[1] <- "geneID"

gene.HO.b <- as.data.frame(unique(pheno12HO.b$geneID))
names(gene.HO.b)[1] <- "geneID"

gene.G <- rbind(gene.D.G, gene.12.G, gene.HO.G)
gene.G <- unique(gene.G)
gene.b <- rbind(gene.D.b, gene.12.b, gene.HO.b)
gene.b <- unique(gene.b)

write.csv(gene.G, "data/GEMMA_files/05_compMethods/AllGEMMA_list.csv")
write.csv(gene.b, "data/GEMMA_files/05_compMethods/AllbigRR_list.csv")
#this file indexes T4 gene names to PFAM annotation terms
FuncAnnot <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")
