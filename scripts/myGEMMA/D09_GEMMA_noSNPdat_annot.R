#Nicole E Soltis
#format SNP peaks for SNPdat annotation - functions and genes
#100217

#from bigRR 09_SNPdat_annot.R
#-----------------------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#load files
#this includes all SNP with p < 99th percentile threshold in the top 1000 (small p values) for at least 1/12 genotypes
## can replace for 99.9% thr - check file
plant12.snp <- read.csv("data/GEMMA_files/D_08_results/12Plants_top1kSNPs_MAF20_10NA_GEMMA_kmat1_99thr.csv")
#do need wide top 1000 for gene annotation and figure S3b: make that in this script

#this includes all SNP with p < 99th %ile threshold across > 6/12 genotypes
##check file for threshold
plantHO.snp <- read.csv("data/GEMMA_files/D_08_results/12Plants_HiOverlapSNPs_trueMAF20_10NA_GEMMA_1kpermut99Thr_kmat1.csv")

#this includes all SNP with p < 0.01 for D, W, or S
## check file for threshold
domest.snp <- read.csv("data/GEMMA_files/D_08_results/GEMMA_peaksDWS_kmat1_99thr.csv")

#SNPdat is being stupid so I'll try to do it myself
setwd("~/Projects")
my.gtf <- read.table("BcGenome/data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE) #replace missing gene names with NAs

my.gtf <- my.gtf[,1:14]

#calculate gene center
#calculate distance gene center to SNP
#add gene with min distance
#range +-1 kb around each snp: lowrange toprange
#match snp chromosome.id to gene V1
#1:18 but have no sig SNPs on chr 17, 18

#do need to keep p-scores: 
#for snp-level overlap, have "SUMM" to say how many phenotypes are significantly associated with a particular SNP
#but for gene-level overlap, need to check again in next script (B07) how many phenotypes have a significant association with any SNP in that GENE
plant12.snp <- plant12.snp[,c(2,3,45,58,11,14,17,20,23,26,29,32,35,38,41,44)]

#associate each plant SNP with nearest gene from my.gtf (this is B05.10 gene annotation)
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant12.sub <- plant12.snp[plant12.snp$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2

  for (i in c(1:nrow(plant12.sub))){
    this.snp <- as.numeric(plant12.sub[i,2])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant12.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, plant12.genes <- all.genes, plant12.genes <- rbind(plant12.genes, all.genes))
}
#plant12.genes now has all SNPs with gene annotations 

plantHO.snp <- plantHO.snp[,c(2,3,58,45,11,14,17,20,23,26,29,32,35,38,41,44)]
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant.sub <- plantHO.snp[plantHO.snp$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (i in c(1:nrow(plant.sub))){
    this.snp <- as.numeric(plant.sub[i,2])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, plantHO.genes <- all.genes, plantHO.genes <- rbind(plantHO.genes, all.genes))
}
#plantHO.genes now has all SNPs with gene annotations

domest.snp <- domest.snp[,c(3,4,11,12,6,8,10)]
#domesticated!
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant.sub <- domest.snp[domest.snp$chr==j,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (i in c(1:nrow(plant.sub))){
    this.snp <- as.numeric(plant.sub[i,2])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant.sub[i,], this.gene)
    ifelse(i == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(j == 1, domest.genes <- all.genes, domest.genes <- rbind(domest.genes, all.genes))
}
#domest.genes now has all SNPs with gene annotations

#now only keep genes if nearest end is within 1kb of SNP
domest.genes$closest.end <- pmin(abs(domest.genes$ps - domest.genes$V4),abs(domest.genes$ps - domest.genes$V5)) 
domest.genes.sub <- domest.genes[domest.genes$closest.end < 1000,]

plant12.genes$closest.end <- pmin(abs(plant12.genes$ps - plant12.genes$V4),abs(plant12.genes$ps - plant12.genes$V5)) 
plant12.genes.sub <- plant12.genes[plant12.genes$closest.end < 1000,]

plantHO.genes$closest.end <- pmin(abs(plantHO.genes$ps - plantHO.genes$V4),abs(plantHO.genes$ps - plantHO.genes$V5)) 
plantHO.genes.sub <- plantHO.genes[plantHO.genes$closest.end < 1000,]


## check file names for threshold
write.csv(domest.genes.sub, "BcSolGWAS/data/GEMMA_files/D_08_results/toannot_domestgenes_99thr.csv")
write.csv(plantHO.genes.sub, "BcSolGWAS/data/GEMMA_files/D_08_results/toannot_plant12HOgenes_99thr.csv")
write.csv(plant12.genes.sub, "BcSolGWAS/data/GEMMA_files/D_08_results/toannot_plant12topgenes_99thr.csv")
