#Nicole E Soltis
#format SNP peaks for SNPdat annotation - functions and genes
#100217

#from bigRR 09_SNPdat_annot.R
#-----------------------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#load files
#this includes all SNP with p < 0.01 in the top 1000 (small p values) for at least 1/11 genotypes (1/12 failed)
plant12snp <- read.csv("paper/plots/addGEMMA/12Plants_top1kSNPs_MAF20_10NA_GEMMA_kmat1.csv")
#do need wide top 1000 for gene annotation and figure S3b: make that in this script

#this includes all SNP with p < 0.01 across > 6/11 genotypes (1/12 failed)
plant12snp.HO <- read.csv("paper/plots/addGEMMA/12Plants_HiOverlapSNPs_trueMAF20_10NA_GEMMA_kmat1.csv")

#this includes all SNP with p < 0.01 for D, W, or S
domestsnp <- read.csv("data/GEMMA_files/04_analysis/GEMMA_peaksDWS_kmat1.csv")

#format files for SNPdat 
names(plant12snp.HO)
plant12snp.HO$chromosome.id <- paste("Chromosome",plant12snp.HO$chr, sep="")
plant12snp.HO$position <- plant12snp.HO$ps
plant12snp.HO$mutation <- "A"
plant12snp.HO.snpdat <- plant12snp.HO[,c("chromosome.id", "position", "mutation")]
write.table(plant12snp.HO.snpdat, file="paper/plots/addGEMMA/SNPdat_toAnnot/plant12snp.HO.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#need to generate this file next
names(domestsnp)
domestsnp$chromosome.id <- paste("Chromosome",domestsnp$chromosome.id, sep="")
domestsnp$position <- domestsnp$ps
domestsnp$mutation <- "A"
domestsnp.snpdat <- domestsnp[,c("chromosome.id", "position", "mutation")]
write.table(domestsnp.snpdat, file="paper/plots/addGEMMA/SNPdat_toAnnot/domestsnp.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#need to convert plant12snp from long to wide?
#basically already in correct format
names(plant12snp)
plant12snp$chromosome.id <- paste("Chromosome",plant12snp$chr, sep="")
plant12snp$position <- plant12snp$ps
plant12snp$mutation <- "A"
plant12snp.snpdat <- plant12snp[,c("chromosome.id", "position", "mutation")]
write.table(plant12snp.snpdat, file="paper/plots/addGEMMA/SNPdat_toAnnot/plant12snp.top1k.FORPERL.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
