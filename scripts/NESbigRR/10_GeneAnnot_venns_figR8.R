#Nicole E Soltis
#02242017
#10_GeneAnnotations.R

#Goal: summarize gene annotation files and draw venn diagrams
#Input Files: Domestication_TopSNPs_SegLong_annot.csv, Plants_TopSNPs_SegLong_0224_annotated.csv
#------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
#read in annotation file
DomestAnt <- read.csv("data/GWAS_files/05_annotation/TrueMAF_NAs/Domestication_TopSNPs_SegLong_trueMAF20_10NA.csv")
#GeneAnnots <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_use.csv")
GeneAnnots <- read.csv("data/Annotate/12Plants_Top1000SNPs_SegWide_trueMAF20_20NA_geneAnnot.csv")

IndPlAnt <- read.csv("data/GWAS_files/05_annotation/IndPlants/Plants_TopSNPs_SegLong_0224_annotated.csv")
names(DomestAnt)
#go long to wide

#just for Domestication SNPs
DomestAnt$TotTraits <- ifelse(abs(DomestAnt$Domesticated) > 0 & abs(DomestAnt$Wild) >0 & abs(DomestAnt$DmWoD) >0, "ALL", 
                           ifelse(abs(DomestAnt$Effect.Domesticated) >0 & abs(DomestAnt$Effect.Wild) >0, "DW",
                                  ifelse(abs(DomestAnt$Effect.DmWoD) >0 & abs(DomestAnt$Effect.Wild) >0, "WS",
                                         ifelse(abs(DomestAnt$Effect.Domesticated) >0 & abs(DomestAnt$Effect.DmWoD) >0, "DS",
                                                ifelse(abs(DomestAnt$Effect.Domesticated) >0, "D",
                                                       ifelse(abs(DomestAnt$Effect.Wild) >0, "W", "S"))))))
#yay! now get table for venn diagram
#remove duplicates by Index
DomestAnt <- DomestAnt[!duplicated(DomestAnt$Index),]
table(DomestAnt$TotTraits)

#summarize PLANT gene annotations
#for individual SNPs include Pos and Index
PlANNOT <- IndPlAnt

names(IndPlAnt)
IndPlAnt <- IndPlAnt[,c("Chrom.Segment", "Gene", "Effect", "Plant")]
PlSNP <- reshape(IndPlAnt, idvar = c("Chrom.Segment","Gene"), timevar="Plant", direction = "wide")
PlSNP$TotTraits <- "none"
PlSNP[is.na(PlSNP)] <- 0

names(PlANNOT)
PlANNOT <- PlANNOT[,c("Gene", "Annot1", "Annot2", "Plant")]
#file for hierarchical functions
PlANNOTb <- PlANNOT
PlANNOTb$Count <- 1
PlANNOTb <- reshape(PlANNOTb, idvar = c("Annot1", "Annot2", "Gene"), timevar="Plant", direction="wide")
PlANNOTb[is.na(PlANNOTb)] <- 0
PlANNOTb$PhenoSum <- rowSums(PlANNOTb[,4:14])
write.csv(PlANNOTb, "data/GWAS_files/05_annotation/IndPlants/TopSNPs_IndPlants_justGenes.csv")

#file for gene-level functions
PlANNOTc <- reshape(PlANNOT, idvar = c("Annot1", "Annot2"), timevar="Plant", direction="wide")
PlANNOTc$TotPlants <- "none"
PlANNOTc$TotPlants <- (rowSums(!is.na(PlANNOTc)) -3)
write.csv(PlANNOTc, "paper/plots/ActualPaper/TableSx2_IndPlantFuncs.csv")

PlANNOTb <- PlANNOT[,c("Chrom.Segment", "Gene", "Annot1", "Annot2")]
#remove rows with Annot1 = #N/A
PlANNOTb <- PlANNOTb[PlANNOTb$Annot1!="#N/A" & PlANNOTb$Annot2!="N/A", ]
PlFuncs <- as.data.frame(table(PlANNOTb$Annot2))
write.csv(PlFuncs, "paper/plots/ActualPaper/TableSx_IndPlantFuncs.csv")


names(PlSNP)
#now count number of non-zero values across rows
PlSNP$TotTraits <- (rowSums(PlSNP != 0) - 3)

table(PlSNP$TotTraits)
hist(PlSNP$TotTraits)
library(ggplot2)

duplicated(PlSNP$Gene) #none, so fine.
m <- ggplot(PlSNP, aes(x=TotTraits))

jpeg("paper/plots/ActualPaper/FigR7/Routs/FigR7c_genehist_w.jpg", width=7, height=4, units='in', res=600)
#jpeg("plots/poster/FigR7c_genehist.jpg", width=6, height=4, units='in', res=600)
m + geom_bar(color="grey20")+
  #geom_bar(color="dodgerblue4", fill="dodgerblue4")+
  theme_bw()+
  scale_x_continuous(name="Plant Genotypes per Candidate Gene", breaks=1:8)+
  scale_y_continuous(name="Number of Genes", breaks=c(0,20,40,60,80,100))
dev.off()

#------------------------------------------------------
#DOMESTICATION
#add in functional annotations for domestication genes
names(GeneAnnots)
GeneAntDom <- GeneAnnots[,c(1:3)]
names(GeneAntDom)[1] <- "Gene"
names(GeneAntDom)[2] <- "Annot1"
names(GeneAntDom)[3] <- "Annot2"
#summarize annotations
table(GeneAntDom$Annot2)

GeneAntSum <- reshape(GeneAntDom, idvar = "Annot2", timevar="Gene", direction="wide")
GeneAntSum <- GeneAntDom[!duplicated(GeneAntDom$Gene),]
GenomeWideFuncs <- table(GeneAntSum$Annot2)
write.csv(GenomeWideFuncs, "paper/plots/ActualPaper/TableSx/TableSx_DomestFuncs.csv")

AllDomest <- merge(DomestAnt, GeneAntDom, by="Gene")
AllDomest <- AllDomest[,c(1,3:9)] #remove X column
AllDomestG <- AllDomest[,c("Chromosome", "Gene", "Annot1", "Annot2", "Trait", "Effect")]
GeneDat <- reshape(AllDomestG, idvar = c("Chromosome", "Gene", "Annot1", "Annot2"), timevar = "Trait", direction = "wide")
names(GeneDat)
GeneDat$TotTraits <- "none"
GeneDat[is.na(GeneDat)] <- 0
#ifelse doesn't work on NA so replace them all with 0

#Remove duplicates per Index
SNPDat <- reshape(AllDomest, idvar = c("Chromosome", "Pos", "Index", "Gene", "Annot1", "Annot2"), timevar="Trait", direction = "wide")


SNPDat$TotTraits <- "none"
SNPDat[is.na(SNPDat)] <- 0

#if Effect.Domesticated and Effect.Wild and Effect.DmWoD != NA, call ALL
#going to call DmWoD "S" for Sensitivity

SNPDat$TotTraits <- ifelse(abs(SNPDat$Effect.Domesticated) > 0 & abs(SNPDat$Effect.Wild) >0 & abs(SNPDat$Effect.DmWoD) >0, "ALL", 
                            ifelse(abs(SNPDat$Effect.Domesticated) >0 & abs(SNPDat$Effect.Wild) >0, "DW",
                                   ifelse(abs(SNPDat$Effect.DmWoD) >0 & abs(SNPDat$Effect.Wild) >0, "WS",
                                          ifelse(abs(SNPDat$Effect.Domesticated) >0 & abs(SNPDat$Effect.DmWoD) >0, "DS",
                                                 ifelse(abs(SNPDat$Effect.Domesticated) >0, "D",
                                                        ifelse(abs(SNPDat$Effect.Wild) >0, "W", "S"))))))
#yay! now get table for venn diagram
#remove duplicates by Index
SNPDat <- SNPDat[!duplicated(SNPDat$Index),]
table(SNPDat$TotTraits)

library(eulerr)
#previously: plots/MultiPlot/meta/VennDia_domest_byGene.jpg
jpeg("paper/plots/ActualPaper/FigR8/FigR8_VennDia_domest_bySNP.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(D=(14+86+2+83), W=(14+83+157+15), S=(14+2+323+15), "D&W" =(14+83), "S&D" =(14+2), "S&W" =(14+15), "D&W&S" =14))
plot(fit, fill_opacity=0)
dev.off()

#remove duplicates by Gene
GeneDat <- GeneDat[!duplicated(GeneDat$Gene),]
GeneDat$TotTraits <- ifelse(abs(GeneDat$Effect.Domesticated) > 0 & abs(GeneDat$Effect.Wild) >0 & abs(GeneDat$Effect.DmWoD) >0, "ALL", 
      ifelse(abs(GeneDat$Effect.Domesticated) >0 & abs(GeneDat$Effect.Wild) >0, "DW",
      ifelse(abs(GeneDat$Effect.DmWoD) >0 & abs(GeneDat$Effect.Wild) >0, "WS",
      ifelse(abs(GeneDat$Effect.Domesticated) >0 & abs(GeneDat$Effect.DmWoD) >0, "DS",
      ifelse(abs(GeneDat$Effect.Domesticated) >0, "D",
      ifelse(abs(GeneDat$Effect.Wild) >0, "W", "S"))))))
#yay! now get table for venn diagram
table(GeneDat$TotTraits)

#domestication is first color, pale green. Wild is lilac.
myColors <- c("#2F4F4F", "#9EFA6C", "#AB82FF")

#install.packages("eulerr")
library(eulerr)
#previously: plots/MultiPlot/meta/VennDia_domest_byGene.jpg
jpeg("paper/plots/ActualPaper/FigR8/Routs/VennDiag/FigR8_VennDia_domest_byGene.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(D=(43+7+5+23), W=(43+23+14+11), S=(43+5+39+11), "D&W" =(43+23), "S&D" =(43+5), "S&W" =(43+11), "D&W&S" =43))
plot(fit, fill_opacity=0)
dev.off()
#fill=c("#2F4F4F", "#9EFA6C", "#AB82FF")

require(VennDiagram)
#previously plots/PosNeg_VennSNPS_domestGenes.emf
#myColors <- c("#1C86EE", "#EE7600", "#050505") with orange
# "#9EFA6C", "#AB82FF", "#2F4F4F" with lime green
venn.diagram(list(Domesticated=c(1:72,73:80,81:112,128:140), Wild=c(1:72, 81:112,113:127,141:168), Sensitivity=c(1:72,73:80,113:127,169:241)), fill_opacity=0.3, fill=c("#1C86EE", "#EE7600", "#050505"), cex = 2, cex.axis=2,filename="paper/plots/ActualPaper/FigR8_PosNeg_VennSNPS_domestGenes2.emf")

names(GeneDat)
#GeneDat3 <- GeneDat[,c("Chrom","Segment","Pos","Gene","Annot1","Annot2","Class")]
#GeneDat3 <- reshape(GeneDat3, idvar = c("Chrom","Segment","Gene","Annot1","Annot2","Class"), timevar = "Pos", direction = "wide")
write.csv(GeneDat, "data/GWAS_files/05_annotation/Domesticated/TopSNPs_domest_justGenes.csv")