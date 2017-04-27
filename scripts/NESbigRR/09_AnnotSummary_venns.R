#Nicole E Soltis
#012417
#Summary of gene annotation data from Bcinerea x Tomato GWAS
#----------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Plots: Venn Diagrams for SNP overlaps between these lists ( Phenotypes > threshold for Domestication phenos vs. Single geno phenos vs. Both //  DmWoD vs. Domesticated vs. Wild)

#most of these data come from 06_bigRRplots_meta.R

#---------------------------------------------------------------------
TopSNP.wide.DM.NA20 <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1000SNPs_SegWide_trueMAF20_20NA.csv")
TopSNP.wide.DM.NA10 <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_Top1000SNPs_SegWide_trueMAF20_10NA.csv")
Top50SNP.wide.DM <- TopSNP.wide.DM.NA10
#this is all SNPs
names(Top50SNP.wide.DM)
a <- table(Top50SNP.wide.DM$Cat)

#domestication is first color, pale green. Wild is lilac.
myColors <- c("#2F4F4F", "#9EFA6C", "#AB82FF")

#install.packages("eulerr")
library(eulerr)
jpeg("plots/MultiPlot/meta/VennDia_domest.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(Do=149, Wi=263, Di=596, "Do&Wi" = 136, 
                "Di&Do" = 3, "Di&Wi" = 25, 
                "Do&Wi&Di" = 17))
plot(fit, fill_opacity=0.3, fill=c("#2F4F4F", "#9EFA6C", "#AB82FF"))

require(VennDiagram)
#domesticated = 17
venn.diagram(list(Domesticated=c(1:17,18:20,21:157,183:331), Wild=c(1:17,21:157,158:182,332:594), Sensitivity=c(1:17,18:20,158:182,595:1190)), fill_opacity=0.3, fill=c("#9EFA6C", "#AB82FF", "#2F4F4F"), cex = 2, cex.axis=2,filename="plots/PosNeg_VennSNPS.emf")

#now genes only
names(Top50SNP.wide.DM)

#-------------------------------------------
#more options

#can also use venn() in gplots package
#by gene
jpeg("plots/MultiPlot/meta/VennDia_domest_gene.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(Do=74, Wi=93, Di=129, "Do&Wi" = 56, 
                "Di&Do" = 40, "Di&Wi" = 45, 
                "Do&Wi&Di" = 30))
plot(fit, fill_opacity=0.3)
dev.off()

#summarize DOMESTICATED gene annotations
GeneDat <- read.csv("data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest_geneannot_forR.csv")
names(GeneDat)
GeneDat2 <- reshape(GeneDat, idvar = c("Chrom","Segment","Pos","Gene","TotalTraits","Annot1","Annot2","Class"), timevar = "Trait", direction = "wide")
write.csv(GeneDat2, "data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest_genesummarized.csv")

GeneDat3 <- GeneDat2[,c("Chrom","Segment","Pos","Gene","Annot1","Annot2","Class")]
GeneDat3 <- reshape(GeneDat3, idvar = c("Chrom","Segment","Gene","Annot1","Annot2","Class"), timevar = "Pos", direction = "wide")
write.csv(GeneDat3, "data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest_justGenes.csv")
