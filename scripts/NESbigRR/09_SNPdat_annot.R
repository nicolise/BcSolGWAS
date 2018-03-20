#Nicole E Soltis
#format SNP peaks for SNPdat annotation - functions and genes
#100217

#from bigRR 09_SNPdat_annot.R
#-----------------------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#load files
plant12snp <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/12Plants_Top1000SNPs_SegLong_trueMAF20_10NA.csv")
#do need wide top 1000 for gene annotation and figure S3b
plant12snp.w <- read.csv("data/GWAS_files/05_annotation/12Plants_Top1000SNPs_SegWide_trueMAF20_10NA.csv")
plant12snp.HO <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/12Plants_HiOverlapSNPs_trueMAF20_10NA.csv")

domestsnp <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")

#format files for SNPdat 
names(plant12snp.HO)
plant12snp.HO$chromosome.id <- paste(plant12snp.HO$Chrom, plant12snp.HO$Segment, sep=".")
plant12snp.HO$chromosome.id <- paste("Chromosome",plant12snp.HO$chromosome.id, sep="")
plant12snp.HO$position <- plant12snp.HO$Pos
plant12snp.HO$mutation <- "A"
plant12snp.HO.snpdat <- plant12snp.HO[,c("chromosome.id", "position", "mutation")]
#get rid of ".0"s
plant12snp.HO.snpdat$chromosome.id <- gsub("\\.0$", "", plant12snp.HO.snpdat$chromosome.id)
unique(plant12snp.HO.snpdat$chromosome.id)
write.table(plant12snp.HO.snpdat, file="data/SNPdat_Annotate/MyAnnots/plant12snp.HO.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

names(domestsnp)
domestsnp$chromosome.id <- paste(domestsnp$Chrom, domestsnp$Segment, sep=".")
domestsnp$chromosome.id <- paste("Chromosome",domestsnp$chromosome.id, sep="")
domestsnp$position <- domestsnp$Pos
domestsnp$mutation <- "A"
domestsnp.snpdat <- domestsnp[,c("chromosome.id", "position", "mutation")]
#get rid of ".0"s
domestsnp.snpdat$chromosome.id <- gsub("\\.0$", "", domestsnp.snpdat$chromosome.id)
write.table(domestsnp.snpdat, file="data/SNPdat_Annotate/MyAnnots/domestsnp.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

names(plant12snp.w)
plant12snp.w$chromosome.id <- paste(plant12snp.w$Chrom, plant12snp.w$Segment, sep=".")
plant12snp.w$chromosome.id <- paste("Chromosome",plant12snp.w$chromosome.id, sep="")
plant12snp.w$position <- plant12snp.w$Pos
plant12snp.w$mutation <- "A"
plant12snp.w.snpdat <- plant12snp.w[,c("chromosome.id", "position", "mutation")]
#get rid of ".0"s
plant12snp.w.snpdat$chromosome.id <- gsub("\\.0$", "", plant12snp.w.snpdat$chromosome.id)
write.table(plant12snp.w.snpdat, file="data/SNPdat_Annotate/MyAnnots/plant12snp.top1k.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
