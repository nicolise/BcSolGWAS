#Nicole E Soltis
#02242017
#B07_GEMMA_Genes

#from 10_GeneAnnot_10NA_venns_figR9figR7b.R
#Goal: summarize gene annotation files and draw venn diagrams
#and draw fig R7 histogram
#------------------------------------------------------
rm(list=ls())
library(dplyr); library(ggplot2)
setwd("~/Projects/BcSolGWAS")
#read in annotation file: see B06_GEMMA_SNPdat_annot.R for which
#could not get SNPdat working so modified file to annotate with nearest genes < 1kb myself
#this includes all SNP with p < 0.01 in the top 1000 (small p values) for at least 1/11 genotypes (1/12 failed)
#consistent with bigRR: distance cutoff is 2kb window around the gene body
plant12gen <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/plant12topgenes.csv")
plant12gen <- plant12gen[,-c(1,17,18,22:25,27,29:34)]

#this includes all SNP with p < 0.01 across > 6/11 genotypes (1/12 failed)
plant12gen.HO <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/plant12HOgenes.csv")
plant12gen.HO <- plant12gen.HO[,-c(1,10,18,19,23:26,28,30:35)]

#this includes all SNP with p < 0.01 for D, W, or S
domestgen <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/domestgenes.csv")
domestgen <- domestgen[,-c(1,9,10,14:17,19,21:26)]

#this includes multiple SNPs per gene


#then count number of phenotypes 
names(plant12gen)[20] <- "geneID"

#if any SNP is > 0 for a given gene, will get sum > 0
#else, sum = 0
IPlant2 <- plant12gen %>%
  group_by(geneID) %>%
  summarize(tot_LA0410 = sum(X9_LA0410_pscore < 0.01, na.rm = TRUE),
            tot_LA0480 = sum(X12_LA0480_pscore < 0.01, na.rm = TRUE),
            tot_LA1547 = sum(X1_LA1547_pscore < 0.01, na.rm = TRUE),
            #tot_LA1589 = sum(abs(Effect.LA1589), na.rm = TRUE),
            tot_LA1684 = sum(X3_LA1684_pscore < 0.01, na.rm = TRUE),
            tot_LA2093 = sum(X4_LA2093_pscore < 0.01, na.rm = TRUE),
            tot_LA2176 = sum(X5_LA2176_pscore < 0.01, na.rm = TRUE),
            tot_LA2706 = sum(X6_LA2706_pscore < 0.01, na.rm = TRUE),
            tot_LA3008 = sum(X7_LA3008_pscore < 0.01, na.rm = TRUE),
            tot_LA3475 = sum(X8_LA3475_pscore < 0.01, na.rm = TRUE),
            tot_LA4345 = sum(X10_LA4345_pscore < 0.01, na.rm = TRUE),
            tot_LA4355 = sum(X11_LA4355_pscore < 0.01, na.rm = TRUE))

IPlant3 <- IPlant2
#for tot phenos just make this 0 SNPs or 1 (any SNPs) sig per pheno
for (i in 2:12){
  fxcol = IPlant3[,paste(colnames(IPlant3[i]),sep='')]
  IPlant3[,paste(colnames(IPlant3[i]))] <- ifelse(fxcol > 0, 1, 0)
}

IPlant3$TotPhenos <- (IPlant3$tot_LA0410 + IPlant3$tot_LA0480 + IPlant3$tot_LA1547 + IPlant3$tot_LA1684 + IPlant3$tot_LA2093 + IPlant3$tot_LA2176 + IPlant3$tot_LA2706 + IPlant3$tot_LA3008 + IPlant3$tot_LA3475 + IPlant3$tot_LA4345 + IPlant3$tot_LA4355)
hist(IPlant3$TotPhenos)
table(IPlant3$TotPhenos)

IPlant4 <- IPlant3 %>%
  group_by(TotPhenos) %>%
  summarise(n = n())

write.csv(IPlant3, "paper/plots/addGEMMA/SNPdat_toAnnot/hi12plants_Genes4Func.csv")
#-------------------------------------------------------------------

#remove duplicate SNPs by Index: this just keeps (arbitrarily) the first Gene annotation per SNP
DomestAnt <- domestgen[!duplicated(domestgen$Index),]
#here's the data for SNP level Venn
table(DomestAnt$TotTraits)

#check values. some genes listed multiply (do not use for gene-level venn!)
#how to summarize per gene for Venn?
#can take average within genes. will still be 0 if no sig fx snps, nonzero if 1 or more sig fx snps
names(domestgen)[12] <- "geneID"

#only keep first gene listing per SNP
DoGenAnt <- domestgen[!duplicated(domestgen$Index),]

#this only keeps 1 value per gene
DoGenAnt <- DoGenAnt %>%
  group_by(geneID) %>%
  summarize(min_Domest = min(pscore.D, na.rm = TRUE),
            min_Wild = min(pscore.W, na.rm = TRUE),
            min_DmWoD = min(pscore.S, na.rm = TRUE))

#just for Domestication SNPs
DoGenAnt$TotTraits <- ifelse(DoGenAnt$min_Domest < 0.01 & DoGenAnt$min_Wild < 0.01 & DoGenAnt$min_DmWoD < 0.01, "ALL", 
                              ifelse(DoGenAnt$min_Domest < 0.01 & DoGenAnt$min_Wild < 0.01, "DW",
                                     ifelse(DoGenAnt$min_Wild < 0.01 & DoGenAnt$min_DmWoD < 0.01, "WS",
                                            ifelse(DoGenAnt$min_Domest < 0.01 &  DoGenAnt$min_DmWoD < 0.01, "DS",
                                                   ifelse(DoGenAnt$min_Domest < 0.01, "D",
                                                          ifelse(DoGenAnt$min_Wild < 0.01, "W", "S"))))))

table(DoGenAnt$TotTraits)


write.csv(DoGenAnt, "paper/plots/addGEMMA/SNPdat_toAnnot/domestplants_Genes4Func.csv")

#------------------------------------------------------------------------------------
#now high overlap list

names(plant12gen.HO)[20] <- "geneID"

#if any SNP is > 0 for a given gene, will get sum > 0
#else, sum = 0
library(dplyr)
IPlant2.HO <- plant12gen.HO %>%
  group_by(geneID) %>%
  summarize(tot_LA0410 = sum(X9_LA0410_pscore < 0.01, na.rm = TRUE),
            tot_LA0480 = sum(X12_LA0480_pscore < 0.01, na.rm = TRUE),
            tot_LA1547 = sum(X1_LA1547_pscore < 0.01, na.rm = TRUE),
            #tot_LA1589 = sum(abs(Effect.LA1589), na.rm = TRUE),
            tot_LA1684 = sum(X3_LA1684_pscore < 0.01, na.rm = TRUE),
            tot_LA2093 = sum(X4_LA2093_pscore < 0.01, na.rm = TRUE),
            tot_LA2176 = sum(X5_LA2176_pscore < 0.01, na.rm = TRUE),
            tot_LA2706 = sum(X6_LA2706_pscore < 0.01, na.rm = TRUE),
            tot_LA3008 = sum(X7_LA3008_pscore < 0.01, na.rm = TRUE),
            tot_LA3475 = sum(X8_LA3475_pscore < 0.01, na.rm = TRUE),
            tot_LA4345 = sum(X10_LA4345_pscore < 0.01, na.rm = TRUE),
            tot_LA4355 = sum(X11_LA4355_pscore < 0.01, na.rm = TRUE))

IPlant3.HO <- IPlant2.HO
#for tot phenos just make this 0 SNPs or 1 (any SNPs) sig per pheno
for (i in 2:12){
  fxcol = IPlant3.HO[,paste(colnames(IPlant3.HO[i]),sep='')]
  IPlant3.HO[,paste(colnames(IPlant3.HO[i]))] <- ifelse(fxcol > 0, 1, 0)
}

IPlant3.HO$TotPhenos <- (IPlant3.HO$tot_LA0410 + IPlant3.HO$tot_LA0480 + IPlant3.HO$tot_LA1547 + IPlant3.HO$tot_LA1684 + IPlant3.HO$tot_LA2093 + IPlant3.HO$tot_LA2176 + IPlant3.HO$tot_LA2706 + IPlant3.HO$tot_LA3008 + IPlant3.HO$tot_LA3475 + IPlant3.HO$tot_LA4345 + IPlant3.HO$tot_LA4355)
hist(IPlant3.HO$TotPhenos)
table(IPlant3.HO$TotPhenos)

IPlant4.HO <- IPlant3.HO %>%
  group_by(TotPhenos) %>%
  summarise(n = n())

write.csv(IPlant3.HO, "paper/plots/addGEMMA/SNPdat_toAnnot/hi12HOplants_Genes4Func.csv")
#------------------------------------------------------------------------------------
library(ggplot2)
IPlant3 <- read.csv("paper/plots/addGEMMA/SNPdat_toAnnot/hi12plants_Genes4Func.csv")
jpeg("paper/plots/addGEMMA/S3B_topgeneOverlap_12Plants_GEMMA.jpg", width=7.5, height=5, units='in', res=600)
ggplot(IPlant3, aes(IPlant3$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "Number of Genes", limits=c(0,1200))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12), limits=c(0,12))
dev.off()

jpeg("paper/plots/addGEMMA/S3B_topgeneOverlap_12Plants_GEMMA_inset.jpg", width=4, height=3, units='in', res=600)
ggplot(IPlant3, aes(IPlant3$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  scale_y_continuous(name= "", limits = c(0,170))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "", breaks=c(7,8,9,10,11,12),labels=c(7,8,9,10,11,12), limits=c(6,13))
dev.off()
