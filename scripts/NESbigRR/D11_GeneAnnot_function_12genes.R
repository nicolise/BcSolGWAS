#11_GeneAnnot_function
#Nicole E Soltis
#052217

#Input: "data/GWAS_files/05_annotation/OverlapCount_plant12topgenes_99thr.csv"
#---------------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
library(dplyr)
Genes12ForAnnot <- read.csv("data/GWAS_files/05_annotation/OverlapCount_plant12topgenes_99thr.csv")
FuncAnnot <- read.csv("data/BcGenome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_mycleaned.csv")

head(FuncAnnot)

#remove T1 from gene names
Genes12ForAnnot$GENE <- gsub('T1$', '', Genes12ForAnnot$V12)


#keep only the columns we need
Funcs <- FuncAnnot[,c("GENE","PFAM_NAME","PFAM_DESCRIPTION")]
#arbitrarily keeping top function per gene
Funcs <- Funcs[!duplicated(Funcs[,"GENE"]),]

#now within each gene, merge back onto AnnotGenes
Gen12Ant <- Genes12ForAnnot
Gen12Ant <- merge(Gen12Ant, Funcs, by="GENE", all.x=TRUE)
#add on Bcin notation
#this file is the gene names from bigRR T4 genes in botportal -- use to match to GEMMA B05.10 genes
matchgenes_botport <- read.csv("data/GEMMA_files/05_compMethods/allBotPort_bigRR_geneLookup.csv")
names(matchgenes_botport)[1] <- "GENE"
matchgenes <- matchgenes_botport[,c(1,5)]
Gen12Ant <- merge(Gen12Ant, matchgenes, by = "GENE", all.x=TRUE)
Gen12Ant <- Gen12Ant[!duplicated(Gen12Ant$GENE),]
#this is a single function per gene
write.csv(Gen12Ant, "data/GWAS_files/05_annotation/window2kb/ReAnnot_All12annots_byGene.csv")
#-------------------------------------------------------------------------
#overrepresentation analysis

#keep only the columns we need
Funcs <- FuncAnnot[,c("GENE","PFAM_NAME","PFAM_DESCRIPTION")]
#could do this if I want to arbitrarily keep 1 function per gene
#Funcs <- Funcs[!duplicated(Funcs[,"GENE"]),]

names(Funcs)[1] <- "GENE"

#now within each gene, merge back onto AnnotGenes
Gen12Ant <- Genes12ForAnnot
Gen12Ant <- merge(Gen12Ant, Funcs, by="GENE")
#p12.unmatch <- pheno12.overlap[is.na(match(pheno12.overlap$T4vanKan.BROAD, Funcs$T4vanKan.BROAD)),]
Gen12Tops <- Gen12Ant[Gen12Ant$TotTraits > 2,]

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

#Summary DFs for phenotypes: overrepresentation of functions test

#how many times does each PFAM show up in this full list?
#this includes multiple genes per function, ignores # phenotypes associated
#this only has 1 observation, so it's fine

#all phenos
#remove unneeded columns, remove repeated rows
#need only TotTraits & PFAM_DESCRIPTION & GENE
Gen12Ant_2 <- Gen12Ant[,-c(2:15,17)]
Gen12Ant_2 <- unique(Gen12Ant_2)
#only keep each gene once: arbitrarily keep top annotation... no, this is a bad approach
Gen12AntCats <- as.data.frame(table(Gen12Ant_2$PFAM_DESCRIPTION))
colnames(Gen12AntCats)[2] <- "AllHOGenFreq"
colnames(Gen12AntCats)[1] <- "Function"
Gen12AntOverrep <- merge(WGAntCats, Gen12AntCats, by="Function", all=T)
#check if any over WG
#Gen12AntOverrep$blah <- Gen12AntOverrep$WGGenFreq - Gen12AntOverrep$AllHOGenFreq
#all fine

#which groups to include?
table(AntCatsSUB$TotTraits)

#summary DF for 12 phenotypes
AntCatsSUB <- Gen12Ant_2
AntCats12 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 12) 
AntCats12$PFAM_DESCRIPTION <- droplevels(AntCats12$PFAM_DESCRIPTION)
#only keep each gene once
AntCats12 <- as.data.frame(table(AntCats12$PFAM_DESCRIPTION))
colnames(AntCats12)[2] <-"GenFreq12p"
colnames(AntCats12)[1]<- "Function"

#now for 11 phenotypes
AntCats11 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 11)
AntCats11$PFAM_DESCRIPTION <- droplevels(AntCats11$PFAM_DESCRIPTION)
AntCats11 <- as.data.frame(table(AntCats11$PFAM_DESCRIPTION))
colnames(AntCats11)[2] <-"GenFreq11p"
colnames(AntCats11)[1]<- "Function"

#10
AntCats10 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 10)
AntCats10$PFAM_DESCRIPTION <- droplevels(AntCats10$PFAM_DESCRIPTION)
AntCats10 <- as.data.frame(table(AntCats10$PFAM_DESCRIPTION))
colnames(AntCats10)[2] <-"GenFreq10p"
colnames(AntCats10)[1]<- "Function"

#9
AntCats9 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 9) 
AntCats9$PFAM_DESCRIPTION <- droplevels(AntCats9$PFAM_DESCRIPTION)
AntCats9 <- as.data.frame(table(AntCats9$PFAM_DESCRIPTION))
colnames(AntCats9)[2] <-"GenFreq9p"
colnames(AntCats9)[1]<- "Function"

#8
AntCats8 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 8) 
AntCats8$PFAM_DESCRIPTION <- droplevels(AntCats8$PFAM_DESCRIPTION)
AntCats8 <- as.data.frame(table(AntCats8$PFAM_DESCRIPTION))
colnames(AntCats8)[2] <-"GenFreq8p"
colnames(AntCats8)[1]<- "Function"

#7
AntCats7 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 7) 
AntCats7$PFAM_DESCRIPTION <- droplevels(AntCats7$PFAM_DESCRIPTION)
AntCats7 <- as.data.frame(table(AntCats7$PFAM_DESCRIPTION))
colnames(AntCats7)[2] <-"GenFreq7p"
colnames(AntCats7)[1]<- "Function"

#6
AntCats6 <- subset(AntCatsSUB, AntCatsSUB[ ,"TotTraits"] == 6) 
AntCats6$PFAM_DESCRIPTION <- droplevels(AntCats6$PFAM_DESCRIPTION)
AntCats6 <- as.data.frame(table(AntCats6$PFAM_DESCRIPTION))
colnames(AntCats6)[2] <-"GenFreq6p"
colnames(AntCats6)[1]<- "Function"

#with all=T, can make underrepresentation work by filling in zeroes
HOAntOverrep <- merge(Gen12AntOverrep, AntCats12, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats11, by="Function", all=TRUE)
HOAntOverrep <- merge(HOAntOverrep, AntCats10, by="Function", all=TRUE)
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

##check which threshold
write.csv(HOAntOverrep, "data/GWAS_files/05_annotation/Rannot_HiOverlap12p_FuncOverrep_99thr.csv")
