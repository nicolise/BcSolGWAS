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
#numbers come from "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_domestONLY_wide.csv"
install.packages("eulerr")
library(eulerr)
jpeg("plots/MultiPlot/meta/VennDia_domest.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(Do=154, Wi=223, Di=378, "Do&Wi" = 76, 
                "Di&Do" = 4, "Di&Wi" = 2, 
                "Do&Wi&Di" = 2))
plot(fit, fill_opacity=0.3)
dev.off()
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
