#Nicole E Soltis
#031417
#11_SNPsFASTA.R 
#script to take target SNP lists and T4 whole-genome FASTA to output FASTA files with only candidate SNPs included. These lists can then be used with BLAST2GO for hierarchical functional annotation
#Inputs: from 10_GeneAnnotations.R: TopSNPs_domest_justGenes.csv AND TopSNPs_IndPlants_justGenes.csv AND from T4 WGS: botrytis_cinerea__t4__1_genes.fasta
#
#-----------------------------------------------------------------------------
#first read in all relevant files
rm(list=ls())
setwd("~/Projects/BcSolGWAS")

#read in gene-level SNP files

DomestGenes <- read.csv("data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv")
DomestGenes <- DomestGenes[,-c(1)]
#PlantGenes <- read.csv("GWAS_files/05_annotation/IndPlants/TopSNPs_IndPlants_justGenes.csv")
  
#read in fasta file
library("Biostrings")
FullFASTA <- readDNAStringSet("data/genome/WGS/botrytis_cinerea__t4__1_genes.fasta")
#turn it into a dataframe to make it easy to work with
seq_name <- names(FullFASTA)
sequence <- paste(FullFASTA)
FASTAdf <- data.frame(seq_name, sequence)
#now only take characters before | pipe into a new variable to match with SNPs "Gene"
names(DomestGenes)

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
DomestAnnot <- DomestGenes

#individual domestication phenotypes
stop("feed in selected phenotype here")
DomestAnnot.dom <- subset(DomestAnnot, DomestAnnot[ ,2] != 0) 
DomestAnnot.wild <- subset(DomestAnnot, DomestAnnot[,3] != 0)
DomestAnnot.sens <- subset(DomestAnnot, DomestAnnot[,4] != 0)

DomestAnnot <- DomestAnnot.sens

#all domestication phenotypes:
DomestAnnot <- as.data.frame(DomestAnnot[,c("geneID")])
DomestAnnot <- rename(DomestAnnot, "DomestAnnot[, c(\"geneID\")]" = "Gene")
DomestAnnot <- merge(DomestAnnot, FASTAgetseqs, by="Gene")

#now convert back to FASTA
names(DomestAnnot)
#install.packages("seqinr")
library(seqinr)
write.fasta(as.list(DomestAnnot$sequence),DomestAnnot$Gene,"data/GWAS_files/05_annotation/DomestGenes_sens_NA10.fasta")
