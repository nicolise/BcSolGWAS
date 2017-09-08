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
HOGenesForAnnot <- read.csv("data/GWAS_files/05_annotation/FINAL_2kbWindow/12plants_HO_genesTOANNOT.csv")
FuncAnnot <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")

head(FuncAnnot)

Funcs <- FuncAnnot[,c(2,3,4)]
Funcs <- Funcs[!duplicated(Funcs[,1]),]

#now within each gene, merge back onto AnnotGenes
colnames(Funcs)[1] <- "geneID"
GenesDForAnnot <- GenesDForAnnot[,-c(1)]
DoGenAnt <- merge(GenesDForAnnot, Funcs, by="geneID")
write.csv(DoGenAnt, "paper/plots/TableS1/AllDOannots_byGene.csv")

HOGenAnt <- HOGenesForAnnot
HOGenAnt <- merge(HOGenAnt, Funcs, by="geneID")
write.csv(HOGenAnt, "paper/plots/TableS1/AllHOannots_byGene.csv")

#summary DF for WHOLE GENOME
#remove empty levels of the variable (in case function appears 0 times in whole genome)
Funcs$PFAM_DESCRIPTION <- droplevels(Funcs$PFAM_DESCRIPTION)
#how many times does a given function occur in the whole genome?
WGAntCats <- as.data.frame(table(Funcs$PFAM_DESCRIPTION))
colnames(WGAntCats)[2] <-"WGGenFreq"
colnames(WGAntCats)[1] <- "Function"
WGAntCats <- WGAntCats[-c(1,2),]
#need to add droplevels() to correctly count the ocurrences of each Function, but this doesn't work for underrepresentation
DoGenAnt$PFAM_DESCRIPTION <- droplevels(DoGenAnt$PFAM_DESCRIPTION)

HOGenAnt$PFAM_DESCRIPTION <- droplevels(HOGenAnt$PFAM_DESCRIPTION)

#Summary DFs for phenotypes: overrepresentation of functions test
#all phenos > 6
HOAntCats <- as.data.frame(table(HOGenAnt$PFAM_DESCRIPTION))
HOAntCats <- HOAntCats[-c(1),]
colnames(HOAntCats)[2] <- "AllHOGenFreq"
colnames(HOAntCats)[1] <- "Function"
HOAntOverrep <- merge(WGAntCats, HOAntCats, by="Function", all=T)

#summary DF for 12 phenotypes
AntCatsSUB <- HOGenAnt[,c(1,14,15,16)]
#remove rows when mean_Wild = 0
AntCats12 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 12) 
AntCats12$PFAM_DESCRIPTION <- droplevels(AntCats12$PFAM_DESCRIPTION)
AntCats12 <- as.data.frame(table(AntCats12$PFAM_DESCRIPTION))
colnames(AntCats12)[2] <-"GenFreq12p"
colnames(AntCats12)[1]<- "Function"

#now for 11 phenotypes
AntCats11 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 11) 
AntCats11$PFAM_DESCRIPTION <- droplevels(AntCats11$PFAM_DESCRIPTION)
AntCats11 <- as.data.frame(table(AntCats11$PFAM_DESCRIPTION))
colnames(AntCats11)[2] <-"GenFreq11p"
colnames(AntCats11)[1]<- "Function"

#10
AntCats10 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 10) 
AntCats10$PFAM_DESCRIPTION <- droplevels(AntCats10$PFAM_DESCRIPTION)
AntCats10 <- as.data.frame(table(AntCats10$PFAM_DESCRIPTION))
colnames(AntCats10)[2] <-"GenFreq10p"
colnames(AntCats10)[1]<- "Function"

#9
AntCats9 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 9) 
AntCats9$PFAM_DESCRIPTION <- droplevels(AntCats9$PFAM_DESCRIPTION)
AntCats9 <- as.data.frame(table(AntCats9$PFAM_DESCRIPTION))
colnames(AntCats9)[2] <-"GenFreq9p"
colnames(AntCats9)[1]<- "Function"

#8
AntCats8 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 8) 
AntCats8$PFAM_DESCRIPTION <- droplevels(AntCats8$PFAM_DESCRIPTION)
AntCats8 <- as.data.frame(table(AntCats8$PFAM_DESCRIPTION))
colnames(AntCats8)[2] <-"GenFreq8p"
colnames(AntCats8)[1]<- "Function"

#7
AntCats7 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 7) 
AntCats7$PFAM_DESCRIPTION <- droplevels(AntCats7$PFAM_DESCRIPTION)
AntCats7 <- as.data.frame(table(AntCats7$PFAM_DESCRIPTION))
colnames(AntCats7)[2] <-"GenFreq7p"
colnames(AntCats7)[1]<- "Function"

#6
AntCats6 <- subset(AntCatsSUB, AntCatsSUB[ ,2] == 6) 
AntCats6$PFAM_DESCRIPTION <- droplevels(AntCats6$PFAM_DESCRIPTION)
AntCats6 <- as.data.frame(table(AntCats6$PFAM_DESCRIPTION))
colnames(AntCats6)[2] <-"GenFreq6p"
colnames(AntCats6)[1]<- "Function"

#with all=T, can make underrepresentation work by filling in zeroes
HOAntOverrep <- merge(HOAntOverrep, AntCats12, by="Function", all=T)
HOAntOverrep <- merge(HOAntOverrep, AntCats11, by="Function", all=T)
HOAntOverrep <- merge(HOAntOverrep, AntCats10, by="Function", all=T)
HOAntOverrep <- merge(HOAntOverrep, AntCats9, by="Function", all=T)
HOAntOverrep <- merge(HOAntOverrep, AntCats8, by="Function", all=T)
HOAntOverrep <- merge(HOAntOverrep, AntCats7, by="Function", all=T)
HOAntOverrep <- merge(HOAntOverrep, AntCats6, by="Function", all=T)
HOAntOverrep[is.na(HOAntOverrep)] <- 0

#remove blank rows
HOAntOverrep <- HOAntOverrep[-c(1),]

#test for overrepresentation of a function
#https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
#fisher's exact test
#http://www.biostathandbook.com/fishers.html
#function my list, function whole list, not-function my list, not-function whole list

#next: modify this to run for each phenotype
for (j in c(3:10)){
  fisher.p.over <- NULL
  fisher.p.under <- NULL
  for (i in 1:nrow(HOAntOverrep)){
    myrow <- HOAntOverrep[i,]
    myA <- myrow[,j]
    myB <- myrow$WGGenFreq
    myCmA <- sum(HOAntOverrep[,j])-myA
    myDmB <- sum(HOAntOverrep$WGGenFreq)-myB
    fisher.p.up <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="greater")$p.value
    fisher.p.down <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="less")$p.value
    fisher.p.over <- append(fisher.p.over, fisher.p.up)
    fisher.p.under <- append(fisher.p.under, fisher.p.down)
  }
  assign(paste("fisher.p.over", names(HOAntOverrep)[j], sep="."),fisher.p.over)
  assign(paste("fisher.p.under", names(HOAntOverrep)[j], sep="."),fisher.p.under)
}

#write fisher values into df
names(HOAntOverrep)
HOAntOverrep$fisher.up.All <- fisher.p.over.AllHOGenFreq
HOAntOverrep$fisher.dn.All <- fisher.p.under.AllHOGenFreq
HOAntOverrep$fisher.up.12 <- fisher.p.over.GenFreq12p
HOAntOverrep$fisher.dn.12 <- fisher.p.under.GenFreq12p
HOAntOverrep$fisher.up.11 <- fisher.p.over.GenFreq11p
HOAntOverrep$fisher.dn.11 <- fisher.p.under.GenFreq11p
HOAntOverrep$fisher.up.10 <- fisher.p.over.GenFreq10p
HOAntOverrep$fisher.dn.10 <- fisher.p.under.GenFreq10p
HOAntOverrep$fisher.up.9 <- fisher.p.over.GenFreq9p
HOAntOverrep$fisher.dn.9 <- fisher.p.under.GenFreq9p
HOAntOverrep$fisher.up.8 <- fisher.p.over.GenFreq8p
HOAntOverrep$fisher.dn.8 <- fisher.p.under.GenFreq8p
HOAntOverrep$fisher.up.7 <- fisher.p.over.GenFreq7p
HOAntOverrep$fisher.dn.7 <- fisher.p.under.GenFreq7p
HOAntOverrep$fisher.up.6 <- fisher.p.over.GenFreq6p
HOAntOverrep$fisher.dn.6 <- fisher.p.under.GenFreq6p

write.csv(HOAntOverrep, "data/GWAS_files/05_annotation/FINAL_2kbWindow/HiOverlap12p_AnnotatedGenes.csv")


#-----------------------------------------------------------------------
#Domestication Phenotypes
#summary DF for ANY phenotype
#if there is a row in the table, it must be significant for at least one pheno, so can count all rows
DoAntCats <- as.data.frame(table(DoGenAnt$PFAM_DESCRIPTION))
DoAntCats <- DoAntCats[-c(1),]
colnames(DoAntCats)[2] <-"AllDoGenFreq"
colnames(DoAntCats)[1]<- "Function"
AntOverrep <- merge(WGAntCats, DoAntCats, by="Function", all=T)

#summary DF for Wild phenotype
WildAntCats <- DoGenAnt[,c(1,3,6,7)]
#remove rows when mean_Wild = 0
WildAntCats <- subset(WildAntCats, WildAntCats[ ,2] != 0) 
WildAntCats$PFAM_DESCRIPTION <- droplevels(WildAntCats$PFAM_DESCRIPTION)
WildAntCats <- as.data.frame(table(WildAntCats$PFAM_DESCRIPTION))
colnames(WildAntCats)[2] <-"WildGenFreq"
colnames(WildAntCats)[1]<- "Function"

#summary DF for Domestication phenotype
DomestAntCats <- DoGenAnt[,c(1,2,6,7)]
DomestAntCats <- subset(DomestAntCats, DomestAntCats[ ,2] != 0) 
DomestAntCats$PFAM_DESCRIPTION <- droplevels(DomestAntCats$PFAM_DESCRIPTION)
DomestAntCats <- as.data.frame(table(DomestAntCats$PFAM_DESCRIPTION))
colnames(DomestAntCats)[2] <-"DomestGenFreq"
colnames(DomestAntCats)[1]<- "Function"

#summary DF for Sensitivity phenotype
SensAntCats <- DoGenAnt[,c(1,4,6,7)]
SensAntCats <- subset(SensAntCats, SensAntCats[ ,2] != 0) 
SensAntCats$PFAM_DESCRIPTION <- droplevels(SensAntCats$PFAM_DESCRIPTION)
SensAntCats <- as.data.frame(table(SensAntCats$PFAM_DESCRIPTION))
colnames(SensAntCats)[2] <-"SensGenFreq"
colnames(SensAntCats)[1]<- "Function"

#with all=T, can make underrepresentation work by filling in zeroes
AntOverrep <- merge(AntOverrep, WildAntCats, by="Function", all=T)
AntOverrep <- merge(AntOverrep, DomestAntCats, by="Function", all=T)
AntOverrep <- merge(AntOverrep, SensAntCats, by="Function", all=T)
AntOverrep[is.na(AntOverrep)] <- 0

#remove blank rows
AntOverrep <- AntOverrep[-c(1,2),]

#test for overrepresentation of a function
#https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
#fisher's exact test
#http://www.biostathandbook.com/fishers.html
#function my list, function whole list, not-function my list, not-function whole list

#next: modify this to run for each phenotype

for (j in c(3:6)){
  fisher.p.over <- NULL
  fisher.p.under <- NULL
  for (i in 1:nrow(AntOverrep)){
    myrow <- AntOverrep[i,]
    myA <- myrow[,j]
    myB <- myrow$WGGenFreq
    myCmA <- sum(AntOverrep[,j])-myA
    myDmB <- sum(AntOverrep$WGGenFreq)-myB
    fisher.p.up <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="greater")$p.value
    fisher.p.down <- fisher.test(matrix(c(myA, myB, myCmA, myDmB), nrow=2, ncol=2), alternative="less")$p.value
    fisher.p.over <- append(fisher.p.over, fisher.p.up)
    fisher.p.under <- append(fisher.p.under, fisher.p.down)
  }
  assign(paste("fisher.p.over", names(AntOverrep)[j], sep="."),fisher.p.over)
  assign(paste("fisher.p.under", names(AntOverrep)[j], sep="."),fisher.p.under)
}

#write fisher values into df
names(AntOverrep)
AntOverrep$fisher.up.All <- fisher.p.over.AllDoGenFreq
AntOverrep$fisher.up.Do <- fisher.p.over.DomestGenFreq
AntOverrep$fisher.up.Wi <- fisher.p.over.WildGenFreq
AntOverrep$fisher.up.Se <- fisher.p.over.SensGenFreq
AntOverrep$fisher.dn.All <- fisher.p.under.AllDoGenFreq
AntOverrep$fisher.dn.Do <- fisher.p.under.DomestGenFreq
AntOverrep$fisher.dn.Wi <- fisher.p.under.WildGenFreq
AntOverrep$fisher.dn.Se <- fisher.p.under.SensGenFreq

write.csv(AntOverrep, "data/GWAS_files/05_annotation/FINAL_2kbWindow/Domestication_AnnotatedGenes.csv")
