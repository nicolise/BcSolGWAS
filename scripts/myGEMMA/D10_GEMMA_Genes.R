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
## read in 99% Thr or 99.9% Thr
plant12gen <- read.csv("data/GEMMA_files/D_08_results/toannot_plant12topgenes_99thr.csv")
plant12gen <- plant12gen[,-c(1,18,19,23:26,28,30:31)]

#this includes all SNP with p < 0.01 across > 6/11 genotypes (1/12 failed)
##check file
plant12gen.HO <- read.csv("data/GEMMA_files/D_08_results/toannot_plant12HOgenes_99thr.csv")
plant12gen.HO <- plant12gen.HO[,-c(1,18,19,23:26,28,30:31)]

#this includes all SNP with p < 0.01 for D, W, or S
##check file
domestgen <- read.csv("data/GEMMA_files/D_08_results/toannot_domestgenes_99thr.csv")
domestgen <- domestgen[,-c(1,9,10,14:17,21)]

#this includes multiple SNPs per gene

#then count number of phenotypes 
names(plant12gen)[21] <- "geneID"

#if any SNP is > 0 for a given gene, will get sum > 0
#else, sum = 0

#get thresholds here 
mythrs <- read.csv("data/GEMMA_files/D_07_randOUTS/GEMMA_1krand_thresholds.csv")
#now, PosPhenos = 1 if p < avg 2500th SNP 1000x thr (99%), 0 if p > thr
# also calculate/ report in text for 99.9% Thr = 250
gethr <- mythrs[mythrs$SNPnum==2500,] ##check for correct level

IPlant2 <- plant12gen %>%
  group_by(geneID) %>%
  summarize(tot_LA3008 = sum(X9_LA3008_pscore < gethr[11,3], na.rm = TRUE),
            tot_LA4355 = sum(X12_LA4355_pscore < gethr[14,3], na.rm = TRUE),
            tot_LA0410 = sum(X1_LA0410_pscore < gethr[3,3], na.rm = TRUE),
            tot_LA0480 = sum(X2_LA0480_pscore < gethr[4,3], na.rm = TRUE),
            tot_LA1547 = sum(X3_LA1547_pscore < gethr[5,3], na.rm = TRUE),
            tot_LA1589 = sum(X4_LA1589_pscore < gethr[6,3], na.rm = TRUE),
            tot_LA1684 = sum(X5_LA1684_pscore < gethr[7,3], na.rm = TRUE),
            tot_LA2093 = sum(X6_LA2093_pscore < gethr[8,3], na.rm = TRUE),
            tot_LA2176 = sum(X7_LA2176_pscore < gethr[9,3], na.rm = TRUE),
            tot_LA2706 = sum(X8_LA2706_pscore < gethr[10,3], na.rm = TRUE),
            tot_LA3475 = sum(X10_LA3475_pscore < gethr[12,3], na.rm = TRUE),
            tot_LA4345 = sum(X11_LA4345_pscore < gethr[13,3], na.rm = TRUE))

IPlant3 <- IPlant2
#for tot phenos just make this 0 SNPs or 1 (any SNPs) sig per pheno
for (i in 2:13){
  fxcol = IPlant3[,paste(colnames(IPlant3[i]),sep='')]
  IPlant3[,paste(colnames(IPlant3[i]))] <- ifelse(fxcol > 0, 1, 0)
}

IPlant3$TotPhenos <- (IPlant3$tot_LA0410 + IPlant3$tot_LA0480 + IPlant3$tot_LA1547 + IPlant3$tot_LA1589 + IPlant3$tot_LA1684 + IPlant3$tot_LA2093 + IPlant3$tot_LA2176 + IPlant3$tot_LA2706 + IPlant3$tot_LA3008 + IPlant3$tot_LA3475 + IPlant3$tot_LA4345 + IPlant3$tot_LA4355)
hist(IPlant3$TotPhenos)
table(IPlant3$TotPhenos)

IPlant4 <- IPlant3 %>%
  group_by(TotPhenos) %>%
  summarise(n = n())

## check file name
write.csv(IPlant3, "data/GEMMA_files/D_08_results/hi12plants_Genes4Func_99thr.csv")
#-------------------------------------------------------------------

#remove duplicate SNPs by Index: this just keeps (arbitrarily) the first gene annotation per SNP
DomestAnt <- domestgen[!duplicated(domestgen$Index),]

#here's the data for SNP level Venn - one gene per SNP
table(DomestAnt$TotTraits)

#check values. some genes listed multiply (do not use for gene-level venn!)
#how to summarize per gene for Venn?
#can take average within genes. will still be 0 if no sig fx snps, nonzero if 1 or more sig fx snps
#replace V12 with Index
names(domestgen)[13] <- "geneID"

#only keep first gene listing per SNP, arbitrarily -- just for Venn
DoGenAnt <- domestgen[!duplicated(domestgen$Index),]

#this only keeps 1 value per gene
DoGenAnt <- DoGenAnt %>%
  group_by(geneID) %>%
  summarize(min_Domest = min(pscore.D, na.rm = TRUE),
            min_Wild = min(pscore.W, na.rm = TRUE),
            min_DmWoD = min(pscore.S, na.rm = TRUE))

#just for Domestication SNPs
DoGenAnt$TotTraits <- ifelse(DoGenAnt$min_Domest < gethr[2,3] & DoGenAnt$min_Wild < gethr[15,3] & DoGenAnt$min_DmWoD < gethr[1,3], "ALL", 
                              ifelse(DoGenAnt$min_Domest < gethr[2,3] & DoGenAnt$min_Wild < gethr[15,3], "DW",
                                     ifelse(DoGenAnt$min_Wild < gethr[15,3] & DoGenAnt$min_DmWoD < gethr[1,3], "WS",
                                            ifelse(DoGenAnt$min_Domest < gethr[2,3] &  DoGenAnt$min_DmWoD < gethr[1,3], "DS",
                                                   ifelse(DoGenAnt$min_Domest < gethr[2,3], "D",
                                                          ifelse(DoGenAnt$min_Wild < gethr[15,3], "W", "S"))))))

table(DoGenAnt$TotTraits)

##check output file
write.csv(DoGenAnt, "data/GEMMA_files/D_08_results/domestplants_Genes4Func_99thr.csv")

#------------------------------------------------------------------------------------
#now high overlap list
## troubleshoot this for 99.9% Thr, seems too low

#rename V12 as geneID
names(plant12gen.HO)[21] <- "geneID"

#if any SNP is > 0 for a given gene, will get sum > 0
#else, sum = 0
library(dplyr)
IPlant2.HO <- plant12gen.HO %>%
  group_by(geneID) %>%
  summarize(tot_LA3008 = sum(X9_LA3008_pscore < gethr[11,3], na.rm = TRUE),
            tot_LA4355 = sum(X12_LA4355_pscore < gethr[14,3], na.rm = TRUE),
            tot_LA0410 = sum(X1_LA0410_pscore < gethr[3,3], na.rm = TRUE),
            tot_LA0480 = sum(X2_LA0480_pscore < gethr[4,3], na.rm = TRUE),
            tot_LA1547 = sum(X3_LA1547_pscore < gethr[5,3], na.rm = TRUE),
            tot_LA1589 = sum(X4_LA1589_pscore < gethr[6,3], na.rm = TRUE),
            tot_LA1684 = sum(X5_LA1684_pscore < gethr[7,3], na.rm = TRUE),
            tot_LA2093 = sum(X6_LA2093_pscore < gethr[8,3], na.rm = TRUE),
            tot_LA2176 = sum(X7_LA2176_pscore < gethr[9,3], na.rm = TRUE),
            tot_LA2706 = sum(X8_LA2706_pscore < gethr[10,3], na.rm = TRUE),
            tot_LA3475 = sum(X10_LA3475_pscore < gethr[12,3], na.rm = TRUE),
            tot_LA4345 = sum(X11_LA4345_pscore < gethr[13,3], na.rm = TRUE))

IPlant3.HO <- IPlant2.HO
#for tot phenos just make this 0 SNPs or 1 (any SNPs) sig per pheno
for (i in 2:13){
  fxcol = IPlant3.HO[,paste(colnames(IPlant3.HO[i]),sep='')]
  IPlant3.HO[,paste(colnames(IPlant3.HO[i]))] <- ifelse(fxcol > 0, 1, 0)
}

IPlant3.HO$TotPhenos <- (IPlant3.HO$tot_LA0410 + IPlant3.HO$tot_LA0480 + IPlant3.HO$tot_LA1547 + IPlant3.HO$tot_LA1589 + IPlant3.HO$tot_LA1684 + IPlant3.HO$tot_LA2093 + IPlant3.HO$tot_LA2176 + IPlant3.HO$tot_LA2706 + IPlant3.HO$tot_LA3008 + IPlant3.HO$tot_LA3475 + IPlant3.HO$tot_LA4345 + IPlant3.HO$tot_LA4355)
hist(IPlant3.HO$TotPhenos)
table(IPlant3.HO$TotPhenos)

IPlant4.HO <- IPlant3.HO %>%
  group_by(TotPhenos) %>%
  summarise(n = n())

write.csv(IPlant3.HO, "data/GEMMA_files/D_08_results/hi12HOplants_Genes4Func_99thr.csv")
#------------------------------------------------------------------------------------
library(ggplot2)
## check input threshold
IPlant3 <- read.csv("data/GEMMA_files/D_08_results/hi12plants_Genes4Func_99thr.csv")
## check output threshold
jpeg("paper/plots/addGEMMA/S3B_topgeneOverlap_12Plants_GEMMA_99thr.jpg", width=7.5, height=5, units='in', res=600)
ggplot(IPlant3, aes(IPlant3$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  ## check axis height
  #scale_y_continuous(name= "Number of Genes", limits=c(0,2000))+
  scale_y_continuous(name= "Number of Genes", limits=c(0,1000))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "Plant Genotypes per Candidate Gene", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12), limits=c(0,13))
dev.off()

##check output value
jpeg("paper/plots/addGEMMA/S3B_topgeneOverlap_12Plants_GEMMA_inset_99thr.jpg", width=4, height=3, units='in', res=600)
ggplot(IPlant3, aes(IPlant3$TotPhenos)) + 
  geom_bar()+
  theme_bw()+
  ## check axis height
  #scale_y_continuous(name= "", limits = c(0,11))+
  scale_y_continuous(name= "", limits = c(0,200))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  scale_x_continuous(name= "", breaks=c(7,8,9,10,11,12),labels=c(7,8,9,10,11,12), limits=c(6,13))
dev.off()
