#10_troubleshootGenes

#SNP list seems to be catching too many genes (10-25% of total)

#check 1: how long are genes in annotation file?

rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
AllGenes <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")

hist(AllGenes$LENGTH) #this is PFAM_STOP - PFAM_START
LongGenes <- AllGenes[AllGenes$LENGTH > 2000,]
hist(LongGenes$LENGTH)


#check 2: how large are windows around SNPs?
SNP12windows <- read.csv("data/SNPdat_Annotate/12Plants_Top1000SNPs_SegWide_trueMAF20_10NA_FORPERL.output.csv")

SNP12windows$Distance.to.nearest.feature <- as.numeric(SNP12windows$Distance.to.nearest.feature)
hist(SNP12windows$Distance.to.nearest.feature)

#could remove all SNPs with Distance > 1000 to nearest future
SNPdist <- SNP12windows[SNP12windows$Distance.to.nearest.feature > 0,]
#max currently is 4168 bp

#check 3: are some SNPs catching multiple genes?
