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
#99.9% GEMMA threshold
pheno12.overlap <- read.csv("data/GEMMA_files/D_08_results/ToAnnot_pheno12_overlap_999thr.csv")

head(FuncAnnot)

#keep only the columns we need
Funcs <- FuncAnnot[,c("GENE","PFAM_NAME","PFAM_DESCRIPTION")]
#could do this if I want to arbitrarily keep 1 function per gene
#Funcs <- Funcs[!duplicated(Funcs[,"GENE"]),]

names(Funcs)[1] <- "T4vanKan.BROAD"

#now within each gene, merge back onto AnnotGenes
Gen12Ant <- pheno12.overlap
Gen12Ant <- merge(Gen12Ant, Funcs, by="T4vanKan.BROAD")
p12.unmatch <- pheno12.overlap[is.na(match(pheno12.overlap$T4vanKan.BROAD, Funcs$T4vanKan.BROAD)),]
Gen12Ant$TotPhenos.O <- Gen12Ant$TotPhenos.b + Gen12Ant$TotPhenos.G
Gen12Tops <- Gen12Ant[Gen12Ant$TotPhenos.b > 1,]
Gen12Tops <- Gen12Tops[Gen12Tops$TotPhenos.G > 1,]
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
Gen12Ant$PFAM_DESCRIPTION <- droplevels(Gen12Ant$PFAM_DESCRIPTION)
Gen12Ant$T4vanKan.BROAD <- droplevels(Gen12Ant$T4vanKan.BROAD)

Gen12Ant$TotPhenos.O <- (Gen12Ant$TotPhenos.b + Gen12Ant$TotPhenos.G)

#Summary DFs for phenotypes: overrepresentation of functions test

#all HO list phenos
#how many times does each PFAM show up in this full list?
  #this includes multiple genes per function, ignores # phenotypes associated

#all phenos
#store number of repeats of each gene name-- depends on how many functions associate with it, nothing else
#remove unneeded columns, remove repeated rows
#only keep gene names, totphenos, PFAM_DESC
Gen12Ant_2 <- Gen12Ant[,-c(2,4:20,22:33,35)]
Gen12Ant_2 <- unique(Gen12Ant_2)
#only keep each gene once: arbitrarily keep top annotation... no, this is a bad approach
Gen12AntCats <- as.data.frame(table(Gen12Ant_2$PFAM_DESCRIPTION))
colnames(Gen12AntCats)[2] <- "AllHOGenFreq"
colnames(Gen12AntCats)[1] <- "Function"
Gen12AntOverrep <- merge(WGAntCats, Gen12AntCats, by="Function", all=T)

#which overlaps to include
table(Gen12Ant_2$TotPhenos.O)
#I'll do 16:10

#summary DF for 16 phenotypes
AntCatsSUB <- Gen12Ant_2[,c("T4vanKan.BROAD","TotPhenos.O","PFAM_DESCRIPTION")]
AntCats16 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 16) 
AntCats16$PFAM_DESCRIPTION <- droplevels(AntCats16$PFAM_DESCRIPTION)
#only keep each gene once
AntCats16 <- as.data.frame(table(AntCats16$PFAM_DESCRIPTION))
colnames(AntCats16)[2] <-"GenFreq16p"
colnames(AntCats16)[1]<- "Function"

#now for 15 phenotypes
AntCats15 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 15)
AntCats15$PFAM_DESCRIPTION <- droplevels(AntCats15$PFAM_DESCRIPTION)
AntCats15 <- as.data.frame(table(AntCats15$PFAM_DESCRIPTION))
colnames(AntCats15)[2] <-"GenFreq15p"
colnames(AntCats15)[1]<- "Function"

#14
AntCats14 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 14)
AntCats14$PFAM_DESCRIPTION <- droplevels(AntCats14$PFAM_DESCRIPTION)
AntCats14 <- as.data.frame(table(AntCats14$PFAM_DESCRIPTION))
colnames(AntCats14)[2] <-"GenFreq14p"
colnames(AntCats14)[1]<- "Function"

#13
AntCats13 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 13) 
AntCats13$PFAM_DESCRIPTION <- droplevels(AntCats13$PFAM_DESCRIPTION)
AntCats13 <- as.data.frame(table(AntCats13$PFAM_DESCRIPTION))
colnames(AntCats13)[2] <-"GenFreq13p"
colnames(AntCats13)[1]<- "Function"

#12
AntCats12 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 12) 
AntCats12$PFAM_DESCRIPTION <- droplevels(AntCats12$PFAM_DESCRIPTION)
AntCats12 <- as.data.frame(table(AntCats12$PFAM_DESCRIPTION))
colnames(AntCats12)[2] <-"GenFreq12p"
colnames(AntCats12)[1]<- "Function"

#11
AntCats11 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 11) 
AntCats11$PFAM_DESCRIPTION <- droplevels(AntCats11$PFAM_DESCRIPTION)
AntCats11 <- as.data.frame(table(AntCats11$PFAM_DESCRIPTION))
colnames(AntCats11)[2] <-"GenFreq11p"
colnames(AntCats11)[1]<- "Function"

#10
AntCats10 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotPhenos.O"] == 10) 
AntCats10$PFAM_DESCRIPTION <- droplevels(AntCats10$PFAM_DESCRIPTION)
AntCats10 <- as.data.frame(table(AntCats10$PFAM_DESCRIPTION))
colnames(AntCats10)[2] <-"GenFreq10p"
colnames(AntCats10)[1]<- "Function"

#with all=T, can make underrepresentation work by filling in zeroes
HOAntOverrep <- merge(Gen12AntOverrep, AntCats16, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats15, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats14, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats13, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats12, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats11, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats10, by="Function", all=TRUE)
HOAntOverrep[is.na(HOAntOverrep)] <- 0

#test for overrepresentation of a function
#https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
#fisher's exact test
#http://www.biostathandbook.com/fishers.html
#function my list, function whole list, not-function my list, not-function whole list

#check which columns have overlap phenotypes = j
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

HOAntOverrep$fisher.up.16 <- fisher.p.over.GenFreq16p
HOAntOverrep$fisher.dn.16 <- fisher.p.under.GenFreq16p
HOAntOverrep$fisher.up.15 <- fisher.p.over.GenFreq15p
HOAntOverrep$fisher.dn.15 <- fisher.p.under.GenFreq15p
HOAntOverrep$fisher.up.14 <- fisher.p.over.GenFreq14p
HOAntOverrep$fisher.dn.14 <- fisher.p.under.GenFreq14p
HOAntOverrep$fisher.up.13 <- fisher.p.over.GenFreq13p
HOAntOverrep$fisher.dn.13 <- fisher.p.under.GenFreq13p
HOAntOverrep$fisher.up.12 <- fisher.p.over.GenFreq12p
HOAntOverrep$fisher.dn.12 <- fisher.p.under.GenFreq12p
HOAntOverrep$fisher.up.11 <- fisher.p.over.GenFreq11p
HOAntOverrep$fisher.dn.11 <- fisher.p.under.GenFreq11p
HOAntOverrep$fisher.up.10 <- fisher.p.over.GenFreq10p
HOAntOverrep$fisher.dn.10 <- fisher.p.under.GenFreq10p

##check which threshold
write.csv(HOAntOverrep, "D_08_results/HiOverlap12p_FuncOverrep_999thr.csv")

#-----------------------------------------------------------------------
