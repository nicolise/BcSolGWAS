#Nicole E Soltis
#051118
#hand-annotating ANY OVERLAP gene list for Figure R5b (previously only annotated top 1000 SNPs/ phenotype or High Overlap (>6 phenotypes) SNPs)
#-----------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
#this is from the 99% positive threshold
fullSNP <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/12Plants_AnyOverlapSNPs_trueMAF20_10NA.csv")

#SNPdat is being stupid so I'll try to do it myself
setwd("~/Projects")
#here's the T4 gtf
my.gtf <- read.table("BcSolGWAS/data/SNPdat_Annotate/genes_Chromosome.gtf", fill=TRUE) #replace missing gene names with NAs

#keep only uniquely named genes
num.genes <- my.gtf[unique(my.gtf$V12),]

#calculate gene center
#calculate distance gene center to SNP
#add gene with min distance
#range +-1 kb around each snp: lowrange toprange
#match snp chromosome.id to gene V1
#1:18 but have no sig SNPs on chr 17, 18

#do need to keep p-scores: 
#for snp-level overlap, have "SUMM" to say how many phenotypes are significantly associated with a particular SNP
#but for gene-level overlap, need to check again in next script (B07) how many phenotypes have a significant association with any SNP in that GENE
plant12.snp <- fullSNP[,c(2:4,17,18:41)]
#do these each because I'm lazy
plant12.snp$LA0480 <- plant12.snp$LA0480neg + plant12.snp$LA0480pos
plant12.snp$LA0410 <- plant12.snp$LA410neg + plant12.snp$LA410pos
plant12.snp$LA1547 <- plant12.snp$LA1547neg + plant12.snp$LA1547pos
plant12.snp$LA1589 <- plant12.snp$LA1589neg + plant12.snp$LA1589pos
plant12.snp$LA1684 <- plant12.snp$LA1684neg + plant12.snp$LA1684pos
plant12.snp$LA2093 <- plant12.snp$LA2093neg + plant12.snp$LA2093pos
plant12.snp$LA2176 <- plant12.snp$LA2176neg + plant12.snp$LA2176pos
plant12.snp$LA2706 <- plant12.snp$LA2706neg + plant12.snp$LA2706pos
plant12.snp$LA3008 <- plant12.snp$LA3008neg + plant12.snp$LA3008pos
plant12.snp$LA3475 <- plant12.snp$LA3475neg + plant12.snp$LA3475pos
plant12.snp$LA4345 <- plant12.snp$LA4345neg + plant12.snp$LA4345pos
plant12.snp$LA4355 <- plant12.snp$LA4355neg + plant12.snp$LA4355pos
plant12.snp <- plant12.snp[,c(1:4,29:40)]

#need Chr.Segment
#Make plotting variables
plant12.snp <- plant12.snp[order(plant12.snp$Chrom, plant12.snp$Segment, plant12.snp$Pos),]
plant12.snp$Chrom.Seg <- paste("Chromosome",plant12.snp$Chrom,".",plant12.snp$Segment, sep="")
plant12.snp$Chrom.Seg <- gsub("\\.0", "", plant12.snp$Chrom.Seg)

#associate each plant SNP with nearest gene from my.gtf (this is B05.10 gene annotation)
blah <- Sys.time()
Sys.time()
for (k in 1:56){
  currChr <- unique(my.gtf$V1)[k]
  gtf.sub <- my.gtf[my.gtf$V1==currChr,]
  plant12.sub <- plant12.snp[plant12.snp$Chrom.Seg==currChr,]
  #within these matched sets...
  gtf.sub$midgene <- (gtf.sub$V4 + gtf.sub$V5)/2
  
  for (m in c(1:nrow(plant12.sub))){
    this.snp <- as.numeric(plant12.sub[m,"Pos"])
    gtf.sub$genedist <- abs(gtf.sub$midgene - this.snp)
    my.closest.gene <- which(gtf.sub$genedist == min(gtf.sub$genedist))
    this.gene <- gtf.sub[my.closest.gene,]
    this.line <- cbind(plant12.sub[m,], this.gene)
    ifelse(m == 1, all.genes <- this.line, all.genes <- rbind(all.genes, this.line))
  }
  ifelse(k == 1, plant12.genes <- all.genes, plant12.genes <- rbind(plant12.genes, all.genes))
  print(Sys.time())
}
blah
Sys.time() #took 20 mins total, not bad

#plant12.genes now has all SNPs with gene annotations 
plant12.genes$closest.end <- pmin(abs(plant12.genes$Pos - plant12.genes$V4),abs(plant12.genes$Pos - plant12.genes$V5)) 
plant12.genes.sub <- plant12.genes[plant12.genes$closest.end < 1000,]

## check file names for threshold
write.csv(plant12.genes.sub, "BcSolGWAS/data/GWAS_files/05_annotation/toannot_plant12topgenes_99thr.csv")

#determine phenotype overlap for these genes 
#first, summarize by gene
#just number of sig SNPs and gene name
plant12g.plot <- plant12.genes.sub[,c(5:16,29)]
plant12g.plot <- plant12g.plot[!duplicated(plant12g.plot),]
#summarize by gene
library(dplyr)
plant12g.plot2 <- plant12g.plot %>%
  group_by(V12) %>%
  summarise(LA0480 = sum(LA0480), LA0410 = sum(LA0410), LA1547 = sum(LA1547), LA1589 = sum(LA1589), LA1684 = sum(LA1684), LA2093 = sum(LA2093), LA2176 = sum(LA2176), LA2706 = sum(LA2706), LA3008 = sum(LA3008), LA3475 = sum(LA3475), LA4345 = sum(LA4345), LA4355 = sum(LA4355))
plant12g.plot3 <- plant12g.plot2

#now filter to 0, 1 at phenotype level
for( i in 2:13 ){
  for (j in 1:nrow(plant12g.plot3)){
    ifelse(plant12g.plot3[j,i]>0,plant12g.plot3[j,i] <- 1, plant12g.plot3[j,i] <- 0)
  }
}

#now sum phenotypes per gene
plant12g.plot3$TotTraits <- plant12g.plot3$LA0480 + plant12g.plot3$LA0410 + plant12g.plot3$LA1547 + plant12g.plot3$LA1589 + plant12g.plot3$LA1684 + plant12g.plot3$LA2093 + plant12g.plot3$LA2176 + plant12g.plot3$LA2706 + plant12g.plot3$LA3008 + plant12g.plot3$LA3475 + plant12g.plot3$LA4345 + plant12g.plot3$LA4355

table(plant12g.plot3$TotTraits)
write.csv(plant12g.plot3, "data/GWAS_files/05_annotation/OverlapCount_plant12topgenes_99thr.csv")

plant12g.plot3 <- read.csv("data/GWAS_files/05_annotation/OverlapCount_plant12topgenes_99thr.csv")
#--------------------------------------------------------------------------------------
#plots Figure 5b

#previous: plots top 12plants genes
setwd("~/Projects/BcSolGWAS")
#here is the gene annotation overlap from only SNPdat. My R annotation finds more gene overlap: this is all I will plot.
HOplant2 <- read.csv("data/GWAS_files/05_annotation/window2kb/12plants_HO_genesTOANNOT.csv")
FullGenes <- read.csv("data/GWAS_files/05_annotation/OverlapCount_plant12topgenes_99thr.csv")
table(FullGenes$TotTraits)

#now need to add BCINs and PFAMs to these genes for Table S3

jpeg("paper/plots/addGEMMA/R5b_topGenesOverlap_IndPlants_2kbWin_Small.jpg", width=4, height=3, units='in', res=600)
ggplot(FullGenes, aes(FullGenes$TotTraits)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(7,8,9,10,11,12),labels=c(7,8,9,10,11,12), limits=c(6,13))
dev.off()

#plot for Fig R5b
IPlant2 <- FullGenes[FullGenes$TotTraits < 7,]
IPlant2 <- IPlant2[,c("V12","TotTraits")]
names(IPlant2) <- c("geneID","TotPhenos")
IPlant2 <- rbind(IPlant2, HOplant2[,c(2,15)])

jpeg("paper/plots/addGEMMA/R5b_topGenesOverlap_IndPlants_2kbWin.jpg", width=7.5, height=5, units='in', res=600)
ggplot(plant12g.plot3, aes(plant12g.plot3$TotTraits)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12))
dev.off()

jpeg("paper/plots/addGEMMA/R5b_topGenesOverlap_IndPlants_2kbWin_Small.jpg", width=4, height=3, units='in', res=600)
ggplot(HOplant2, aes(HOplant2$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(7,8,9,10,11,12),labels=c(7,8,9,10,11,12), limits=c(6,13))
dev.off()

#here's mixed plot to look... DO NOT actually include this in paper.
ggplot(IPlant2, aes(IPlant2$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12))