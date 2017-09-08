#Nicole E Soltis
#02242017
#10_GeneAnnotations.R

#Goal: summarize gene annotation files and draw venn diagrams
#Input Files: Domestication_TopSNPs_SegLong_annot.csv, Plants_TopSNPs_SegLong_0224_annotated.csv

#and draw fig R7 histogram
#------------------------------------------------------
rm(list=ls())
library(dplyr); library(ggplot2)
setwd("~/Projects/BcSolGWAS")
#read in annotation file
DomestAnt <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")
#this is just the list of genes, from SNPdat. SNP info etc. must be reattached.
DoGenes <- read.csv("data/SNPdat_Annotate/Final_annots/Domestication_TopSNPs_SegLong_trueMAF20_10NA.output.csv")
names(DomestAnt)

#includes the top 1000 SNPs per plant genotype (all >99% Thr)
IPGenes <- read.csv("data/SNPdat_Annotate/Final_annots/12Plants_Top1000SNPs_SegWide_trueMAF20_10NA_FORPERL.output.csv")
IndPlAnt <- read.csv("data/SNPdat_Annotate/Final_annots/12Plants_Top1000SNPs_SegWide_trueMAF20_10NA.csv")

#top overlap SNPs
HOGenes <- read.csv("data/SNPdat_Annotate/Final_annots/12Plants_HiOverlapSNPs_trueMAF20_10NA_output.csv")
HOSNP <- read.csv("data/SNPdat_Annotate/Final_annots/12Plants_HiOverlapSNPs_trueMAF20_10NA.csv")

#go long to wide

#for IndPlant: merge with genes
IndPlAnt$Chrom2 <- paste("CHROMOSOME",IndPlAnt$Chrom, sep='')
IndPlAnt$Chrom.Pos <- paste(IndPlAnt$Chrom2, IndPlAnt$Pos, sep='.')
IPGenes$Chrom.Pos <- paste(IPGenes$Chromosome.Number, IPGenes$SNP.Position, sep='.')
#now, need to narrow down IPGenes list to only include genes WITHIN WINDOW of SNP
#window options to try:
#1kb (500 bp each side), 2kb (1kb each side), 4kb (2kb each side)
IPGenes$Distance.to.nearest.feature <- as.numeric(IPGenes$Distance.to.nearest.feature)
IPGenes$Distance.to.nearest.feature[is.na(IPGenes$Distance.to.nearest.feature)] <-0
IPgenp5 <- IPGenes[IPGenes$Distance.to.nearest.feature<250,]
IPgen1 <- IPGenes[IPGenes$Distance.to.nearest.feature<500,]
IPgen2 <- IPGenes[IPGenes$Distance.to.nearest.feature<1000,]
IPgen4 <- IPGenes[IPGenes$Distance.to.nearest.feature<2000,]

#sub in IPgen1, IPgen2, IPgen4
IPgen <- IPgen2[,c(25,11)]
#this includes multiple SNPs per gene
IPgen <- unique(IPgen)
IndPlAnt <- IndPlAnt[,-c(1,18)]
IPant <- merge(IndPlAnt, IPgen, by="Chrom.Pos")

#then count number of phenotypes 
names(IPant)
colnames(IPant)[18] <- "geneID"

#if any SNP is > 0 for a given gene, will get sum > 0
#else, sum = 0
IPant2 <- IPant %>%
  group_by(geneID) %>%
  summarize(tot_LA0410 = sum(abs(Effect.LA410), na.rm = TRUE),
            tot_LA0480 = sum(abs(Effect.LA0480), na.rm = TRUE),
            tot_LA1547 = sum(abs(Effect.LA1547), na.rm = TRUE),
            tot_LA1589 = sum(abs(Effect.LA1589), na.rm = TRUE),
            tot_LA1684 = sum(abs(Effect.LA1684), na.rm = TRUE),
            tot_LA2093 = sum(abs(Effect.LA2093), na.rm = TRUE),
            tot_LA2176 = sum(abs(Effect.LA2176), na.rm = TRUE),
            tot_LA2706 = sum(abs(Effect.LA2706), na.rm = TRUE),
            tot_LA3008 = sum(abs(Effect.LA3008), na.rm = TRUE),
            tot_LA3475 = sum(abs(Effect.LA3475), na.rm = TRUE),
            tot_LA4345 = sum(abs(Effect.LA4345), na.rm = TRUE),
            tot_LA4355 = sum(abs(Effect.LA4355), na.rm = TRUE))

for (i in 2:13){
  fxcol = IPant2[,paste(colnames(IPant2[i]),sep='')]
  IPant2[,paste(colnames(IPant2[i]))] <- ifelse(fxcol > 0, 1, 0)
  #assign(IPant2[i], fxcol)
}

IPant2$TotPhenos <- (IPant2$tot_LA0410 + IPant2$tot_LA0480 + IPant2$tot_LA1547 + IPant2$tot_LA1589 + IPant2$tot_LA1684 + IPant2$tot_LA2093 + IPant2$tot_LA2176 + IPant2$tot_LA2706 + IPant2$tot_LA3008 + IPant2$tot_LA3475 + IPant2$tot_LA4345 + IPant2$tot_LA4355)
hist(IPant2$TotPhenos)
table(IPant2$TotPhenos)

IPant3 <- IPant2 %>%
  group_by(TotPhenos) %>%
  summarise(n = n())

write.csv(IPant2, "data/GWAS_files/05_annotation/FINAL_2kbWindow/12plants_genesTOANNOT.csv")

jpeg("paper/plots/FigR7/R7b_topGenesOverlap_IndPlants_1kbWin_Small.jpg", width=4, height=3, units='in', res=600)
ggplot(IPant2, aes(IPant2$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(6,7,8,9,10,11,12),labels=c(6,7,8,9,10,11,12), limits=c(5,13))
dev.off()
#-------------------------------------------------------------------

#just for Domestication SNPs
DomestAnt$TotTraits <- ifelse(abs(DomestAnt$Domesticated) > 0 & abs(DomestAnt$Wild) >0 & abs(DomestAnt$DmWoD) >0, "ALL", 
                           ifelse(abs(DomestAnt$Domesticated) >0 & abs(DomestAnt$Wild) >0, "DW",
                                  ifelse(abs(DomestAnt$DmWoD) >0 & abs(DomestAnt$Wild) >0, "WS",
                                         ifelse(abs(DomestAnt$Domesticated) >0 & abs(DomestAnt$DmWoD) >0, "DS",
                                                ifelse(abs(DomestAnt$Domesticated) >0, "D",
                                                       ifelse(abs(DomestAnt$Wild) >0, "W", "S"))))))
#yay! now get table for venn diagram
#remove duplicate SNPs by Index
DomestAnt <- DomestAnt[!duplicated(DomestAnt$Index),]
#SNP level venn! Consistent with the figure
table(DomestAnt$TotTraits)

#now add genes in
#first, subset by distance
DoGenes$Distance_to_nearest_feature <- as.numeric(DoGenes$Distance_to_nearest_feature)
#multiple listings per gene, with different positional information
#so: cut off list to remove distant genes (>1kb from SNP, 2kb windows) THEN remove duplicate genes
#up to 10kb, crazy
DoGenes$Distance_to_nearest_feature[is.na(DoGenes$Distance_to_nearest_feature)] <-0
hist(DoGenes$Distance_to_nearest_feature)
#DoGenp5 <- DoGenes[DoGenes$Distance_to_nearest_feature<250,]
#DoGen1 <- DoGenes[DoGenes$Distance_to_nearest_feature<500,]
DoGen2 <- DoGenes[DoGenes$Distance_to_nearest_feature<1000,]
#DoGen4 <- DoGenes[DoGenes$Distance_to_nearest_feature<2000,]

#changes this out: DoGen1, DoGen2, DoGen4
stop("replace the correct DoGen list and delete this stop")
DoGenes <- DoGen2

#keep only chrom, snp, and gene id
DoGenes <- DoGenes[,c(1,2,11)]
DoGenes$geneID <- DoGenes$gene_ID_containing_the_current_feature
DoGenes$Chrom <- DoGenes$Chromosome_Number
DoGenes$Pos <- DoGenes$SNP_Position
DoGenes$Chrom.Pos <- paste(DoGenes$Chrom, DoGenes$Pos, sep='.')
#now only Chrom, Pos, Chrom.Pos and geneID
DoGenes <- DoGenes[,c(5,6,7,4)]
#only Chrom.Pos and geneID
GeneList <- DoGenes[,c(3,4)]
#make sure those are unique: choosing only one Gene per SNP
GeneList <- GeneList[!duplicated(GeneList$Chrom.Pos),]

DomestAnt$Chrom2<- paste("CHROMOSOME",DomestAnt$Chrom, sep='')
DomestAnt$Chrom.Pos<- paste(DomestAnt$Chrom2,DomestAnt$Pos, sep='.')
DomestAnt2 <- DomestAnt[,c(2,3,4,5,6,7,11,14)]

#here is where I actually attach Gene ID
#this gives one D, W, S value for each 
DoGenAnt <- merge(DomestAnt2, GeneList, by="Chrom.Pos")

#and for annotation list:
DomestAnt <- DomestAnt[,-c(1)]
DomestAnt <- DomestAnt[,c(1,2,3,4,5,6,10,11,13)]
GenesforAnnot <- merge(DomestAnt, GeneList, by="Chrom.Pos")

#check venn values
table(GenesforAnnot$TotTraits)

write.csv(GenesforAnnot, "data/GWAS_files/05_annotation/FINAL_2kbWindow/Domest_genestoANNOT.csv")

#how to summarize per gene for Venn?
#can take average within genes. will still be 0 if no sig fx snps, nonzero if 1 or more sig fx snps

DoGenAnt <- DoGenAnt %>%
  group_by(geneID) %>%
  summarize(mean_Domest = mean(Domesticated, na.rm = TRUE),
            mean_Wild = mean(Wild, na.rm = TRUE),
            mean_DmWoD = mean(DmWoD, na.rm = TRUE))

#just for Domestication SNPs
DoGenAnt$TotTraits <- ifelse(abs(DoGenAnt$mean_Domest) > 0 & abs(DoGenAnt$mean_Wild) >0 & abs(DoGenAnt$mean_DmWoD) >0, "ALL", 
                              ifelse(abs(DoGenAnt$mean_Domest) >0 & abs(DoGenAnt$mean_Wild) >0, "DW",
                                     ifelse(abs(DoGenAnt$mean_DmWoD) >0 & abs(DoGenAnt$mean_Wild) >0, "WS",
                                            ifelse(abs(DoGenAnt$mean_Domest) >0 & abs(DoGenAnt$mean_DmWoD) >0, "DS",
                                                   ifelse(abs(DoGenAnt$mean_Domest) >0, "D",
                                                          ifelse(abs(DoGenAnt$mean_Wild) >0, "W", "S"))))))

table(DoGenAnt$TotTraits)

write.csv(DoGenAnt, "data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv")

#------------------------------------------------------------------------------------
#for HiOverlapGenes: merge with genes
#need to match on chromosome.pos for both dataframes
head(HOSNP)
head(HOGenes)
HOSNP$Chrom2 <- paste("CHROMOSOME",HOSNP$Chrom, sep='')
HOSNP$Chrom.Pos <- paste(HOSNP$Chrom2, HOSNP$Pos, sep='.')
HOGenes$Chrom.Pos <- paste(HOGenes$Chromosome.Number, HOGenes$SNP.Position, sep='.')
#now, need to narrow down IPGenes list to only include genes WITHIN WINDOW of SNP
#window options to try:
#1kb (500 bp each side), 2kb (1kb each side), 4kb (2kb each side)
#sticking with 2kb (1kb each side)
HOGenes$Distance.to.nearest.feature <- as.numeric(HOGenes$Distance.to.nearest.feature)
HOGenes$Distance.to.nearest.feature[is.na(HOGenes$Distance.to.nearest.feature)] <-0
HOgen2 <- HOGenes[HOGenes$Distance.to.nearest.feature<1000,]


#sub in IPgen1, IPgen2, IPgen4
HOgen <- HOgen2[,c(25,11)] 
#this includes multiple SNPs per gene
HOgen <- unique(HOgen) #this way each gene will only get annotated to this SNP set once. Should not lose any genes this way.
HOSNP <- HOSNP[,-c(1)]
HOplant <- merge(HOSNP, HOgen, by="Chrom.Pos") #this includes multiple SNPs per gene again
#also count number of phenotypes PER GENE
names(HOplant)
colnames(HOplant)[46] <- "geneID"
#for SNP info per gene, only keep first mention of each gene

#now summarize by gene (just counting phenotypes per gene)
HOplant2 <- HOplant %>%
  group_by(geneID) %>%
  summarize(tot_LA0410 = sum(LA410pos, LA410neg, na.rm = TRUE),
            tot_LA0480 = sum(LA0480pos, LA0480neg, na.rm = TRUE),
            tot_LA1547 = sum(LA1547pos, LA1547neg, na.rm = TRUE),
            tot_LA1589 = sum(LA1589pos, LA1589neg, na.rm = TRUE),
            tot_LA1684 = sum(LA1684pos, LA1684neg, na.rm = TRUE),
            tot_LA2093 = sum(LA2093pos, LA2093neg, na.rm = TRUE),
            tot_LA2176 = sum(LA2176pos, LA2176neg, na.rm = TRUE),
            tot_LA2706 = sum(LA2706pos, LA2706neg, na.rm = TRUE),
            tot_LA3008 = sum(LA3008pos, LA3008neg, na.rm = TRUE),
            tot_LA3475 = sum(LA3475pos, LA3475neg, na.rm = TRUE),
            tot_LA4345 = sum(LA4345pos, LA4345neg, na.rm = TRUE),
            tot_LA4355 = sum(LA4355pos, LA4345neg, na.rm = TRUE))
names(HOplant2)

for (i in 2:13){
  fxcol = HOplant2[,paste(colnames(HOplant2[i]),sep='')]
  HOplant2[,paste(colnames(HOplant2[i]))] <- ifelse(fxcol > 0, 1, 0)
  #assign(IPant2[i], fxcol)
}

HOplant2$TotPhenos <- (HOplant2$tot_LA0410 + HOplant2$tot_LA0480 + HOplant2$tot_LA1547 + HOplant2$tot_LA1589 + HOplant2$tot_LA1684 + HOplant2$tot_LA2093 + HOplant2$tot_LA2176 + HOplant2$tot_LA2706 + HOplant2$tot_LA3008 + HOplant2$tot_LA3475 + HOplant2$tot_LA4345 + HOplant2$tot_LA4355)
hist(HOplant2$TotPhenos)
table(HOplant2$TotPhenos)

#now, merge on gene annotations for SNP info:
names(HOGenes)
colnames(HOGenes)[11] <- "geneID"
finalGen <- HOGenes[!duplicated(HOGenes$geneID),]
FullGeneInfo <- merge(HOplant2, finalGen, by="geneID")

#gene level phenotype overlap for top overlapping SNPs!
HOplant3 <- HOplant2 %>%
  group_by(TotPhenos) %>%
  summarise(n = n())

write.csv(HOplant2, "data/GWAS_files/05_annotation/FINAL_2kbWindow/12plants_HO_genesTOANNOT.csv")
write.csv(FullGeneInfo, "data/GWAS_files/05_annotation/FINAL_2kbWindow/12plants_HO_geneINFO.csv")

jpeg("paper/plots/FigR7/R7b_topGenesOverlap_IndPlants_1kbWin.jpg", width=7.5, height=5, units='in', res=600)
ggplot(IPant2, aes(IPant2$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12))
dev.off()
#-----------------------------------
