#Nicole E Soltis 
#041017
#local annotation of trueMAF genes using Suzi's lists

#----------------------------------------------------
#prepare work environment
remove(list=ls())
setwd("~/Projects/BcSolGWAS")
#read in annotation file
DomestAnt <- read.csv("results/Domestication_TopSNPs_SegWide_trueMAF.csv")
DomestAnt_oldMAF <- read.csv("results/Domestication_TopSNPs_SegWide_annot.csv")
GeneAnnots <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_use.csv")
#IndPlAnt <- read.csv("data/GWAS_files/05_annotation/IndPlants/Plants_TopSNPs_SegLong_0224_annotated.csv")


#check overlap between SNP lists
#how many SNPs are shared between both?
DomestAnt$Chrom.Seg.Pos <- paste(DomestAnt$Chrom,".",DomestAnt$Segment,".",DomestAnt$Pos)
DomestAnt_oldMAF$Chrom.Seg.Pos <- paste(DomestAnt_oldMAF$Chromosome,".",DomestAnt_oldMAF$Segment,".",DomestAnt_oldMAF$Pos)
length(intersect(DomestAnt$Chrom.Seg.Pos, DomestAnt_oldMAF$Chrom.Seg.Pos))
#1015 retained

#how many SNPs are removed for trueMAF?

#are there any novel SNPs for trueMAF?
#if no, I can match up annotations myself
#if yes, I need Suzi to help with annotations
