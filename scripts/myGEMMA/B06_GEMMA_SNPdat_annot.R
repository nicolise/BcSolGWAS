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
domestsnp$chromosome.id <- paste("Chromosome",domestsnp$chr, sep="")
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


#SNPdat is being stupid so I'll try to do it myself
setwd("~/Projects")
my.gtf <- read.table("BcGenome/data/ensembl/B05.10/extractedgff/Botrytis_cinerea.ASM83294v1.38.chrom.gtf", fill=TRUE) #replace missing gene names with NAs
plant12.snp <- read.table("BcSolGWAS/paper/plots/addGEMMA/SNPdat_toAnnot/plant12snp.top1k.FORPERL.txt", header=TRUE)
plantHO.snp <- read.table("BcSolGWAS/paper/plots/addGEMMA/SNPdat_toAnnot/plant12snp.HO.FORPERL.txt", header=TRUE)
domest.snp <- read.table("BcSolGWAS/paper/plots/addGEMMA/SNPdat_toAnnot/domestsnp.FORPERL.txt", header=TRUE)

my.gtf <- my.gtf[,1:14]

#calculate gene center
#calculate distance gene center to SNNP
#add gene with min distance
#range +-1 kb around each snp: lowrange toprange
#match snp chromosome.id to gene V1
#1:18 but have no sig SNPs on chr 17, 18
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant12.sub <- plant12.snp[plant12.snp$chromosome.id==paste("Chromosome",j,sep=""),]
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


for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant.sub <- plantHO.snp[plantHO.snp$chromosome.id==paste("Chromosome",j,sep=""),]
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

#domesticated!
for (j in 1:16){
  gtf.sub <- my.gtf[my.gtf$V1==paste("Chromosome",j,sep=""),]
  plant.sub <- domest.snp[domest.snp$chromosome.id==paste("Chromosome",j,sep=""),]
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

#now only keep genes if nearest end is within 1kb of SNP
domest.genes$closest.end <- pmin(abs(domest.genes$position - domest.genes$V4),abs(domest.genes$position - domest.genes$V5)) 
domest.genes.sub <- domest.genes[domest.genes$closest.end < 1000,]

plant12.genes$closest.end <- pmin(abs(plant12.genes$position - plant12.genes$V4),abs(plant12.genes$position - plant12.genes$V5)) 
plant12.genes.sub <- plant12.genes[plant12.genes$closest.end < 1000,]

plantHO.genes$closest.end <- pmin(abs(plantHO.genes$position - plantHO.genes$V4),abs(plantHO.genes$position - plantHO.genes$V5)) 
plantHO.genes.sub <- plantHO.genes[plantHO.genes$closest.end < 1000,]

write.csv(domest.genes.sub, "BcSolGWAS/paper/plots/addGEMMA/SNPdat_toAnnot/domestgenes.csv")
write.csv(plantHO.genes.sub, "BcSolGWAS/paper/plots/addGEMMA/SNPdat_toAnnot/plant12HOgenes.csv")
write.csv(plant12.genes.sub, "BcSolGWAS/paper/plots/addGEMMA/SNPdat_toAnnot/plant12topgenes.csv")
