#Nicole E Soltis
#B08_GEMMA_bigRR_overlaps

#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")

## doing overlap with 99 thr for now. If there's lots (seems unlikely...) I can try this with the 999 thr
DWS.G <- read.csv("data/GEMMA_files/D_08_results/domestplants_Genes4Func_999thr.csv")
pheno12.G <- read.csv("data/GEMMA_files/D_08_results/hi12plants_Genes4Func_999thr.csv")
pheno12HO.G <- read.csv("data/GEMMA_files/D_08_results/hi12HOplants_Genes4Func_999thr.csv")

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

##check output file name
#write.csv(gene.G, "data/GEMMA_files/D_08_results/AllGEMMA_geneslist_99thr.csv")
#write.csv(gene.b, "data/GEMMA_files/D_08_results/AllbigRR_geneslist.csv")

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
#--------------------------------------------------------------------------------------
#domesticated first because most interesting
DWS.b <- DWS.b[,-c(1)]
names(DWS.G)
DWS.b <- DWS.b %>% rename(T4vanKan.BROAD = geneID, mean_Domest.b = mean_Domest, mean_Wild.b = mean_Wild, mean_Sens.b = mean_DmWoD, TotTraits.b = TotTraits)
DWS.G <- DWS.G[,-c(1)]
DWS.G <- DWS.G %>% rename(BcinB0510gene = geneID, min_Domest.G = min_Domest, min_Wild.G = min_Wild, min_Sens.G = min_DmWoD, TotTraits.G = TotTraits)
DWS.G$BcinB0510gene <- gsub("gene:", "", DWS.G$BcinB0510gene)

#check for GEMMA genes that are named non-B05.10
DWS.nomatch <- DWS.G[(is.na(match(DWS.G$BcinB0510gene, DWS.overlap$BcinB0510gene))),]
#mismatches of format ENSRNA049514122 or EFBcG00000013384 -- no matches in botport
#done here

DWS.overlap <- merge(matchgenes_botport, DWS.b, by = "T4vanKan.BROAD")
DWS.overlap <- merge(DWS.overlap, DWS.G, by = "BcinB0510gene")

venn.dwsover<- DWS.overlap[!duplicated(DWS.overlap$BcinB0510gene),]
table(venn.dwsover$TotTraits.G)
table(venn.dwsover$TotTraits.b)

#recode traits for gene-level overlap across methods
venn.dwsover$TotTraits.o <- paste(venn.dwsover$TotTraits.b, venn.dwsover$TotTraits.G, sep="")
table(venn.dwsover$TotTraits.o)

#annotate this list!
##check this output file name
write.csv(DWS.overlap, "data/GEMMA_files/D_08_results/ToAnnot_DWS_overlap_999thr.csv")

#---------------------------------------------------------------------------
#now pheno12
pheno12.b <- pheno12.b[,-c(1)]
pheno12.b <- pheno12.b %>% rename(T4vanKan.BROAD = geneID, LA0410.b = tot_LA0410, LA0480.b = tot_LA0480, LA1547.b = tot_LA1547, LA1589.b = tot_LA1589, LA1684.b = tot_LA1684, LA2093.b = tot_LA2093, LA2176.b = tot_LA2176, LA2706.b = tot_LA2706, LA3008.b = tot_LA3008, LA3475.b = tot_LA3475, LA4345.b = tot_LA4345, LA4355.b = tot_LA4355, TotPhenos.b = TotPhenos)
pheno12.G <- pheno12.G[,-c(1)]
pheno12.G <- pheno12.G %>% rename(BcinB0510gene = geneID, LA0410.G = tot_LA0410, LA0480.G = tot_LA0480, LA1547.G = tot_LA1547, LA1589.G = tot_LA1589, LA1684.G = tot_LA1684, LA2093.G = tot_LA2093, LA2176.G = tot_LA2176, LA2706.G = tot_LA2706, LA3008.G = tot_LA3008, LA3475.G = tot_LA3475, LA4345.G = tot_LA4345, LA4355.G = tot_LA4355, TotPhenos.G = TotPhenos)
pheno12.G$BcinB0510gene <- gsub("gene:", "", pheno12.G$BcinB0510gene)

pheno12.overlap <- merge(matchgenes_botport, pheno12.b, by = "T4vanKan.BROAD")
pheno12.overlap <- merge(pheno12.overlap, pheno12.G, by = "BcinB0510gene")

venn.pheno12over <- pheno12.overlap[!duplicated(pheno12.overlap$BcinB0510gene),]
table(venn.pheno12over$TotPhenos.G)
table(venn.pheno12over$TotPhenos.b)

#recode traits for gene-level overlap across methods
venn.pheno12over$TotPhenos.o <- venn.pheno12over$TotPhenos.b + venn.pheno12over$TotPhenos.G
table(venn.pheno12over$TotPhenos.o)

#annotate this list!
##check this output file name
write.csv(pheno12.overlap, "data/GEMMA_files/D_08_results/ToAnnot_pheno12_overlap_999thr.csv")

#---------------------------------------------------------------------------
#now pheno12HO
pheno12HO.b <- pheno12HO.b[,-c(1)]
pheno12HO.b <- pheno12HO.b %>% rename(T4vanKan.BROAD = geneID, LA0410.b = tot_LA0410, LA0480.b = tot_LA0480, LA1547.b = tot_LA1547, LA1589.b = tot_LA1589, LA1684.b = tot_LA1684, LA2093.b = tot_LA2093, LA2176.b = tot_LA2176, LA2706.b = tot_LA2706, LA3008.b = tot_LA3008, LA3475.b = tot_LA3475, LA4345.b = tot_LA4345, LA4355.b = tot_LA4355, TotPhenos.b = TotPhenos)
pheno12HO.G <- pheno12HO.G[,-c(1)]
pheno12HO.G <- pheno12HO.G %>% rename(BcinB0510gene = geneID, LA0410.G = tot_LA0410, LA0480.G = tot_LA0480, LA1547.G = tot_LA1547, LA1589.G = tot_LA1589, LA1684.G = tot_LA1684, LA2093.G = tot_LA2093, LA2176.G = tot_LA2176, LA2706.G = tot_LA2706, LA3008.G = tot_LA3008, LA3475.G = tot_LA3475, LA4345.G = tot_LA4345, LA4355.G = tot_LA4355, TotPhenos.G = TotPhenos)
pheno12HO.G$BcinB0510gene <- gsub("gene:", "", pheno12HO.G$BcinB0510gene)

pheno12HO.overlap <- merge(matchgenes_botport, pheno12HO.b, by = "T4vanKan.BROAD")
pheno12HO.overlap <- merge(pheno12HO.overlap, pheno12HO.G, by = "BcinB0510gene")
venn.pheno12HOover <- pheno12HO.overlap[!duplicated(pheno12HO.overlap$BcinB0510gene),]
table(venn.pheno12HOover$TotPhenos.G)
table(venn.pheno12HOover$TotPhenos.b)
venn.pheno12HOover$TotPhenos.o <- venn.pheno12HOover$TotPhenos.b + venn.pheno12HOover$TotPhenos.G
table(venn.pheno12HOover$TotPhenos.o)

#annotate this list!
##check this output file name
write.csv(pheno12HO.overlap, "data/GEMMA_files/D_08_results/ToAnnot_HOpheno12_overlap_999thr.csv")