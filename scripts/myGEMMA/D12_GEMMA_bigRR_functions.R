#Nicole E Soltis
#04/23/18

#------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
#this file indexes T4 gene names to PFAM annotation terms
FuncAnnot <- read.csv("data/BcGenome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")

#more functions: messy lookup list here. would need to reformat for R
#data/GEMMA_files/05_compMethods/allBotPort_bigRR_functionLookup.xlsx

##check which threshold
#99.9% threshold
pheno12HO.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_HOpheno12_overlap_999thr.csv")
pheno12.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_pheno12_overlap_999thr.csv")
DWS.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_DWS_overlap_999thr.csv")

#99% threshold
#pheno12HO.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_HOpheno12_overlap_99thr.csv")
#pheno12.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_pheno12_overlap_99thr.csv")
#DWS.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_DWS_overlap_99thr.csv")


head(FuncAnnot)

#keep only the columns we need
Funcs <- FuncAnnot[,c("GENE","PFAM_NAME","PFAM_DESCRIPTION")]
#could do this if I want to arbitrarily keep 1 function per gene
#Funcs <- Funcs[!duplicated(Funcs[,"GENE"]),]

names(Funcs)[1] <- "T4vanKan.BROAD"

#now within each gene, merge back onto AnnotGenes
DWS.overlap <- DWS.overlap[,-c(1)]
DoGenAnt <- merge(DWS.overlap, Funcs, by="T4vanKan.BROAD", all.x=TRUE)
#keep the list of genes with no function match in the T4 list
DWS.unmatch <- DWS.overlap[is.na(match(DWS.overlap$T4vanKan.BROAD, Funcs$T4vanKan.BROAD)),]
#could search for these in data/GEMMA_files/05_compMethods/allBotPort_bigRR_functionLookup.xlsx but for the most part these entries have no functional annotation!
##check file name
#keep only the first function for each gene
DoGenAnt <- DoGenAnt[!duplicated(DoGenAnt$BcinB0510gene),]
#write.csv(DoGenAnt, "data/GEMMA_files/D_08_results/AllDOfuncs_byGene_999thr.csv")

HOGenAnt <- pheno12HO.overlap
HOGenAnt <- merge(HOGenAnt, Funcs, by="T4vanKan.BROAD")
##check file name
#write.csv(HOGenAnt, "data/GEMMA_files/D_08_results/AllHOannots_byGene_999thr.csv")

Gen12Ant <- pheno12.overlap
Gen12Ant <- merge(Gen12Ant, Funcs, by="T4vanKan.BROAD", all.x=TRUE)
p12.unmatch <- pheno12.overlap[is.na(match(pheno12.overlap$T4vanKan.BROAD, Funcs$T4vanKan.BROAD)),]
Gen12Ant$TotPhenos.O <- Gen12Ant$TotPhenos.b + Gen12Ant$TotPhenos.G
Gen12Tops <- Gen12Ant[Gen12Ant$TotPhenos.b > 1,]
Gen12Tops <- Gen12Tops[Gen12Tops$TotPhenos.G > 1,]
#only keep first mention of each gene
Gen12Tops <- Gen12Tops[!duplicated(Gen12Tops$BcinB0510gene),]
##check file name
#write.csv(Gen12Tops, "data/GEMMA_files/D_08_results/All12annots_byGene_999thr.csv")

#summary DF for WHOLE GENOME
#remove empty levels of the variable (in case function appears 0 times in whole genome)
Funcs$PFAM_DESCRIPTION <- droplevels(Funcs$PFAM_DESCRIPTION)
#across how many genes does a given function occur in the whole (annotated) genome?
WGAntCats <- as.data.frame(table(Funcs$PFAM_DESCRIPTION))
colnames(WGAntCats)[2] <-"WGGenFreq"
colnames(WGAntCats)[1] <- "Function"
WGAntCats <- WGAntCats[-c(1),]
#need to add droplevels() to correctly count the ocurrences of each Function, but this doesn't work for underrepresentation
#same goes for occurrences of each Gene
DoGenAnt$PFAM_DESCRIPTION <- droplevels(DoGenAnt$PFAM_DESCRIPTION)
DoGenAnt$T4vanKan.BROAD<- droplevels(DoGenAnt$T4vanKan.BROAD)
Gen12Ant$PFAM_DESCRIPTION <- droplevels(Gen12Ant$PFAM_DESCRIPTION)
Gen12Ant$T4vanKan.BROAD <- droplevels(Gen12Ant$T4vanKan.BROAD)
HOGenAnt$PFAM_DESCRIPTION <- droplevels(HOGenAnt$PFAM_DESCRIPTION)
HOGenAnt$T4vanKan.BROAD <- droplevels(HOGenAnt$T4vanKan.BROAD)

Gen12Ant$TotPhenos.O <- (Gen12Ant$TotPhenos.b + Gen12Ant$TotPhenos.G)

#Summary DFs for phenotypes: overrepresentation of functions test

#all HO list phenos
#how many times does each PFAM show up in this full list?
  #this includes multiple genes per function, ignores # phenotypes associated
#this only has 1 observation, so it's fine
HOAntCats <- as.data.frame(table(HOGenAnt$PFAM_DESCRIPTION))
colnames(HOAntCats)[2] <- "AllHOGenFreq"
colnames(HOAntCats)[1] <- "Function"
HOAntOverrep <- merge(WGAntCats, HOAntCats, by="Function", all=T)

#all phenos
#store number of repeats of each gene name-- depends on how many functions associate with it, nothing else
#remove unneeded columns, remove repeated rows
Gen12Ant_2 <- Gen12Ant[,-c(2:20,22:33,35)]
Gen12Ant_2 <- unique(Gen12Ant_2)
#only keep each gene once: arbitrarily keep top annotation... no, this is a bad approach
Gen12AntCats <- as.data.frame(table(Gen12Ant_2$PFAM_DESCRIPTION))
colnames(Gen12AntCats)[2] <- "AllHOGenFreq"
colnames(Gen12AntCats)[1] <- "Function"
Gen12AntOverrep <- merge(WGAntCats, Gen12AntCats, by="Function", all=T)

#summary DF for 12 phenotypes
AntCatsSUB <- Gen12Ant_2[,c("T4vanKan.BROAD","TotPhenos.O","PFAM_DESCRIPTION")]
AntCats12 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 12) 
AntCats12$PFAM_DESCRIPTION <- droplevels(AntCats12$PFAM_DESCRIPTION)
#only keep each gene once
AntCats12 <- as.data.frame(table(AntCats12$PFAM_DESCRIPTION))
colnames(AntCats12)[2] <-"GenFreq12p"
colnames(AntCats12)[1]<- "Function"

# #now for 11 phenotypes
# AntCats11 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos"] == 11) 
# AntCats11$PFAM_DESCRIPTION <- droplevels(AntCats11$PFAM_DESCRIPTION)
# AntCats11 <- as.data.frame(table(AntCats11$PFAM_DESCRIPTION))
# colnames(AntCats11)[2] <-"GenFreq11p"
# colnames(AntCats11)[1]<- "Function"
# 
# #10
# AntCats10 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos"] == 10) 
# AntCats10$PFAM_DESCRIPTION <- droplevels(AntCats10$PFAM_DESCRIPTION)
# AntCats10 <- as.data.frame(table(AntCats10$PFAM_DESCRIPTION))
# colnames(AntCats10)[2] <-"GenFreq10p"
# colnames(AntCats10)[1]<- "Function"

#9
AntCats9 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 9) 
AntCats9$PFAM_DESCRIPTION <- droplevels(AntCats9$PFAM_DESCRIPTION)
AntCats9 <- as.data.frame(table(AntCats9$PFAM_DESCRIPTION))
colnames(AntCats9)[2] <-"GenFreq9p"
colnames(AntCats9)[1]<- "Function"

#8
AntCats8 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 8) 
AntCats8$PFAM_DESCRIPTION <- droplevels(AntCats8$PFAM_DESCRIPTION)
AntCats8 <- as.data.frame(table(AntCats8$PFAM_DESCRIPTION))
colnames(AntCats8)[2] <-"GenFreq8p"
colnames(AntCats8)[1]<- "Function"

#7
AntCats7 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 7) 
AntCats7$PFAM_DESCRIPTION <- droplevels(AntCats7$PFAM_DESCRIPTION)
AntCats7 <- as.data.frame(table(AntCats7$PFAM_DESCRIPTION))
colnames(AntCats7)[2] <-"GenFreq7p"
colnames(AntCats7)[1]<- "Function"

#6
AntCats6 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 6) 
AntCats6$PFAM_DESCRIPTION <- droplevels(AntCats6$PFAM_DESCRIPTION)
AntCats6 <- as.data.frame(table(AntCats6$PFAM_DESCRIPTION))
colnames(AntCats6)[2] <-"GenFreq6p"
colnames(AntCats6)[1]<- "Function"

#with all=T, can make underrepresentation work by filling in zeroes
HOAntOverrep <- merge(HOAntOverrep, AntCats12, by="Function", all=TRUE)
#HOAntOverrep <- merge(HOAntOverrep, AntCats11, by="Function", all=TRUE)
#HOAntOverrep <- merge(HOAntOverrep, AntCats10, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats9, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats8, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats7, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats6, by="Function", all=TRUE)
HOAntOverrep[is.na(HOAntOverrep)] <- 0

#test for overrepresentation of a function
#https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
#fisher's exact test
#http://www.biostathandbook.com/fishers.html
#function my list, function whole list, not-function my list, not-function whole list

#check which columns have overlap phenotypes = j
for (j in c(3:8)){
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
#HOAntOverrep$fisher.up.11 <- fisher.p.over.GenFreq11p
#HOAntOverrep$fisher.dn.11 <- fisher.p.under.GenFreq11p
#HOAntOverrep$fisher.up.10 <- fisher.p.over.GenFreq10p
#HOAntOverrep$fisher.dn.10 <- fisher.p.under.GenFreq10p
HOAntOverrep$fisher.up.9 <- fisher.p.over.GenFreq9p
HOAntOverrep$fisher.dn.9 <- fisher.p.under.GenFreq9p
HOAntOverrep$fisher.up.8 <- fisher.p.over.GenFreq8p
HOAntOverrep$fisher.dn.8 <- fisher.p.under.GenFreq8p
HOAntOverrep$fisher.up.7 <- fisher.p.over.GenFreq7p
HOAntOverrep$fisher.dn.7 <- fisher.p.under.GenFreq7p
HOAntOverrep$fisher.up.6 <- fisher.p.over.GenFreq6p
HOAntOverrep$fisher.dn.6 <- fisher.p.under.GenFreq6p

##check which threshold
write.csv(HOAntOverrep, "data/GEMMA_files/D_08_results/HiOverlap12p_FuncOverrep_999thr.csv")

#-----------------------------------------------------------------------
#Domestication Phenotypes
#summary DF for ANY phenotype
#don't bother with which/ how many domest traits. Just DoAntCats.
#if there is a row in the table, it must be significant for at least one pheno, so can count all rows
#remove unneeded columns, remove repeated rows
DoGenAnt_2 <- DoGenAnt[,-c(2:10,12:14,16)]
DoGenAnt_2 <- unique(DoGenAnt_2)

DoAntCats <- as.data.frame(table(DoGenAnt_2$PFAM_DESCRIPTION))
colnames(DoAntCats)[2] <-"AllDoGenFreq"
colnames(DoAntCats)[1]<- "Function"
AntOverrep <- merge(WGAntCats, DoAntCats, by="Function", all=T)
#with all=T, can make underrepresentation work by filling in zeroes
AntOverrep[is.na(AntOverrep)] <- 0

#test for overrepresentation of a function
#fisher's exact test

#determine which columns in AntOverrep are phenotype overlaps
for (j in c(3)){
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
AntOverrep$fisher.dn.All <- fisher.p.under.AllDoGenFreq

##check which threshold
write.csv(AntOverrep, "data/GEMMA_files/D_08_results/Domestication_FuncOverrep_999thr.csv")
