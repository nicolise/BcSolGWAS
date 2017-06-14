#11_GeneAnnot_function
#Nicole E Soltis
#052217

#Input: SNP data/GWAS_files/05_annotation/Domest_topSNPs_genestoAnnot.csv 
#and Gene data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv
#from script 10_GeneAnnot_10NA_venns_figR8.R
#Output:
#Figures: None
#---------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
library(dplyr)
GenesDForAnnot <- read.csv("data/GWAS_files/05_annotation/Domest_TopSNPs_10NA_intoAnt.csv")
Genes12ForAnnot <- read.csv("data/GWAS_files/05_annotation/FINAL_2kbWindow/12plants_genestoANNOT.csv")
FuncAnnot <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_use.csv")

head(FuncAnnot)

Funcs <- FuncAnnot[,c(2,3,4)]
Funcs <- Funcs[!duplicated(Funcs[,1]),]

#now within each gene, merge back onto AnnotGenes
colnames(Funcs)[1] <- "geneID"
GenesDForAnnot <- GenesDForAnnot[,-c(1)]
DoGenAnt <- merge(GenesDForAnnot, Funcs, by="geneID")
#keep only the first mention of each gene
DoGenAnt <- DoGenAnt[!duplicated(DoGenAnt[,1]),]

DoGenAnt$PFAM_DESCRIPTION <- droplevels(DoGenAnt$PFAM_DESCRIPTION)

DoAntCats <- as.data.frame(table(DoGenAnt$PFAM_DESCRIPTION))

Funcs$PFAM_DESCRIPTION <- droplevels(Funcs$PFAM_DESCRIPTION)
WGAntCats <- as.data.frame(table(Funcs$PFAM_DESCRIPTION))

colnames(DoAntCats)[2] <-"DoGenFreq"
colnames(WGAntCats)[2] <-"WGGenFreq"

AntOverrep <- merge(DoAntCats, WGAntCats, by="Var1")
#remove blank rows
AntOverrep <- AntOverrep[-c(1,2),]
sum(AntOverrep[,2])
sum(AntOverrep[,3])
#AntOverrep <- rbind(AntOverrep, c(0001,1249,4485))

AntOverrep[AntOverrep$Var1=="Major Facilitator Superfamily",]
#test for overrepresentation of a function
#https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
#function my list, function whole list, not-function my list, not-function whole list

fisher.col <- NULL

for (i in 1:nrow(AntOverrep)){
  myrow <- AntOverrep[i,]
  myA <- myrow$DoGenFreq
  myB <- myrow$WGGenFreq
  myCmA <- sum(AntOverrep$DoGenFreq)-myA
  myDmB <- sum(AntOverrep$WGGenFreq)-myB
  fisher.p <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="greater")$p.value
  fisher.col <- append(fisher.col, fisher.p)
}

AntOverrep$Fisher.p <- fisher.col


#----------------------------------------------------------------------
#trying to figure this chunk out
#someone's function for fisher test across a dataframe (stackoverflow, un=vodka)
get_fisher <- function(df){
  mat <- matrix(as.numeric(df[c(2:3)]), ncol=2)
  f <- fisher.test(as.table(mat), alt="two.sided")
  return(c(df[1], f$p.value))
}

fishers <- apply(AntOverrep, 1,  get_fisher)

#another option
AntOverrep$p_value <- apply(AntOverrep,1,function(x) fisher.test(matrix(x[-1],nrow=728))$p.value)

fisher.test(matrix(AntOverrep[,-1],nrow=728))

#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html
dhyper(x, m, n, k, log = FALSE)
dhyper(39, 186, 4485-186, 1249, log=F)

#--------------------------------------------------------------------------------

# #summarize annotations
# DoGenAnt %>%
#   group_by(PFAM_DESCRIPTION) %>%
#   summarize(count_genes = count(geneID, na.rm = TRUE))

mytable <- as.data.frame(table(droplevels(DoGenAnt2$PFAM_DESCRIPTION)))
write.csv(mytable, "data/GWAS_files/05_annotation/Domest_NA10_FunctionsSummary.csv")

#Enrichment analysis
#use fisher exact test or hypergeometric test

#first, collapse DoGenAnt table so that there is only one row per annotation
#so: go long to wide format. 
#remove geneID, Chrom.Pos, Chrom, Segment, Pos, Index, Domesticated, Wild, DmWoD, PFAM_NAME
#split TotTraits into new columns (D, W, S... etc)
#keep PFAM_DESCRIPTION
DoAntW <- DoGenAnt[,c(12,10)]
#go wide
DoAntW <- table(PFAM_DESCRIPTION = DoAntW$PFAM_DESCRIPTION, DoAntW$TotTraits)
DoAntW <- as.data.frame(DoAntW)
DoAntW <- reshape(DoAntW, idvar="PFAM_DESCRIPTION", timevar="Var2", direction="wide")
DoAntW <- rename(DoAntW, c("PFAM_DESCRIPTION"= "PFAM_DESCRIPTION", "Freq.D"="D", "Freq.S"="S", "Freq.W"="W", "Freq.WS"="WS", "Freq.DW"="DW", "Freq.DS"="DS"))
DoAntW$DSW <- 

#add a column of "any traits"
DoAntW$Freq.Any <-

