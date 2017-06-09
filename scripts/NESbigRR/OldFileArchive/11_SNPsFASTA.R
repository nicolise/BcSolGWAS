#Nicole E Soltis
#031417
#11_SNPsFASTA.R 
#script to take target SNP lists and T4 whole-genome FASTA to output FASTA files with only candidate SNPs included. These lists can then be used with BLAST2GO for hierarchical functional annotation
#Inputs: from 10_GeneAnnotations.R: TopSNPs_domest_justGenes.csv AND TopSNPs_IndPlants_justGenes.csv AND from T4 WGS: botrytis_cinerea__t4__1_genes.fasta
#
#-----------------------------------------------------------------------------
#first read in all relevant files
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")

#read in gene-level SNP files
DomestGenes <- read.csv("GWAS_files/05_annotation/Domesticated/TopSNPs_domest_justGenes.csv")
PlantGenes <- read.csv("GWAS_files/05_annotation/IndPlants/TopSNPs_IndPlants_justGenes.csv")
  
#read in fasta file
library("Biostrings")
FullFASTA <- readDNAStringSet("genome/WGS/botrytis_cinerea__t4__1_genes.fasta")
#turn it into a dataframe to make it easy to work with
seq_name <- names(FullFASTA)
sequence <- paste(FullFASTA)
FASTAdf <- data.frame(seq_name, sequence)
#now only take characters before | pipe into a new variable to match with SNPs "Gene"
names(DomestGenes)
names(PlantGenes)

x <- "BcT4_1 | Botrytis cinerea (T4) BcT4_1 (1968 nt)"
x <- gsub("\\|.*", "", x)
x <- gsub(" ", "", x, fixed = TRUE)
names(FASTAdf)
FASTAdf$Gene <- FASTAdf$seq_name
FASTAdf$Gene <- gsub("\\|.*", "", FASTAdf$Gene)
FASTAdf$Gene <- gsub(" ", "", FASTAdf$Gene, fixed = TRUE)
names(FASTAdf)
FASTAgetseqs <- FASTAdf[,c("Gene", "sequence")]

#now add sequences onto SNP lists
names(DomestGenes)
DomestAnnot <- as.data.frame(DomestGenes[,c("Gene")])
DomestAnnot <- rename(DomestAnnot, "DomestGenes[, c(\"Gene\")]" = "Gene")
DomestAnnot <- merge(DomestAnnot, FASTAgetseqs, by="Gene")

PlantAnnot <- as.data.frame(PlantGenes[,c("Gene")])
PlantAnnot <- rename(PlantAnnot, "PlantGenes[, c(\"Gene\")]" = "Gene")
PlantAnnot <- merge(PlantAnnot, FASTAgetseqs, by="Gene")

#now convert back to FASTA
names(DomestAnnot)
#install.packages("seqinr")
library(seqinr)
write.fasta(as.list(PlantAnnot$sequence),PlantAnnot$Gene,"GWAS_files/05_annotation/IndPlants/IndPlantGenes.fasta")
write.fasta(as.list(DomestAnnot$sequence),DomestAnnot$Gene,"GWAS_files/05_annotation/Domesticated/DomestGenes.fasta")
