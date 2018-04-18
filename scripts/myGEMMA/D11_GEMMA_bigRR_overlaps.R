#Nicole E Soltis
#B08_GEMMA_bigRR_overlaps

#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")

## doing overlap with 99 thr for now. If there's lots (seems unlikely...) I can try this with the 999 thr
DWS.G <- read.csv("data/GEMMA_files/D_08_results/domestplants_Genes4Func_99thr.csv")
pheno12.G <- read.csv("data/GEMMA_files/D_08_results/hi12plants_Genes4Func_99thr.csv")
pheno12HO.G <- read.csv("data/GEMMA_files/D_08_results/hi12HOplants_Genes4Func_99thr.csv")

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
gene.G <- unique(gene.G) #total of ~6k at 99% thr
gene.b <- rbind(gene.D.b, gene.12.b, gene.HO.b)
gene.b <- unique(gene.b)

write.csv(gene.G, "data/GEMMA_files/D_08_results/AllGEMMA_geneslist_99thr.csv")
write.csv(gene.b, "data/GEMMA_files/D_08_results/AllbigRR_geneslist.csv")
#this file indexes T4 gene names to PFAM annotation terms
FuncAnnot <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")

#this file is the gene names from bigRR T4 genes in botportal -- use to match to GEMMA B05.10 genes
matchgenes_botport <- read.csv("data/GEMMA_files/05_compMethods/allBotPort_bigRR_geneLookup.csv")
#and for all functions: messy lookup list here. would need to reformat for R
#data/GEMMA_files/05_compMethods/allBotPort_bigRR_functionLookup.xlsx

#add botport names to gene.b
names(gene.b)[1] <- "T4vanKan.BROAD"
names(gene.G)[1] <- "BcinB0510gene"

#just sanity check: does merge work for the botport list and the gene names it was generated from (gene.b)
named.genes <- merge(gene.b, matchgenes_botport, by = "T4vanKan.BROAD") #yes

library(dplyr)
#now check for overlaps within phenotypes
#domesticated first because most interesting
DWS.b <- DWS.b[,-c(1)]
names(DWS.G)
DWS.b <- DWS.b %>% rename(T4vanKan.BROAD = geneID, mean_Domest.b = mean_Domest, mean_Wild.b = mean_Wild, mean_Sens.b = mean_DmWoD, TotTraits.b = TotTraits)
DWS.G <- DWS.G[,-c(1)]
DWS.G <- DWS.G %>% rename(BcinB0510gene = geneID, min_Domest.G = min_Domest, min_Wild.G = min_Wild, min_Sens.G = min_DmWoD, TotTraits.G = TotTraits)
DWS.G$BcinB0510gene <- gsub("gene:", "", DWS.G$BcinB0510gene)

DWS.overlap <- merge(matchgenes_botport, DWS.b, by = "T4vanKan.BROAD")
DWS.overlap <- merge(DWS.overlap, DWS.G, by = "BcinB0510gene")

venn.dwsover<- DWS.overlap[unique(DWS.overlap$BcinB0510gene),]
table(venn.dwsover$TotTraits.G)
table(venn.dwsover$TotTraits.b)

#recode traits for gene-level overlap across methods
venn.dwsover$TotTraits.o <- paste(venn.dwsover$TotTraits.b, venn.dwsover$TotTraits.G, sep="")
table(venn.dwsover$TotTraits.o)

#annotate this list!
names(FuncAnnot)[2] <- "T4vanKan.BROAD"
DWS.overlap.Ant <- merge(DWS.overlap, FuncAnnot, by = "T4vanKan.BROAD")
write.csv(DWS.overlap.Ant, "data/GEMMA_files/D_08_results/Annotated_DWS_overlap.csv")

#now gene.12
