#Nicole E Soltis
#041017
#local annotation of trueMAF genes using Suzi's lists

#----------------------------------------------------
#prepare work environment
remove(list=ls())
setwd("~/Projects/BcSolGWAS")

library(tidyr); library(plyr)

#read in annotation file
DomestAnt <- read.csv("results/Domestication_TopSNPs_SegWide_trueMAF.csv")
DomestAnt_oldMAF <- read.csv("data/GWAS_files/05_annotation/Domesticated/Domestication_TopSNPs_SegLong_annot.csv")
GeneAnnots <- read.csv("data/genome/annotation/botrytis_cinerea__t4__1_pfam_to_genes_use.csv")
#IndPlAnt <- read.csv("data/GWAS_files/05_annotation/IndPlants/Plants_TopSNPs_SegLong_0224_annotated.csv")

#convert OldMAF to wide format and fix chromosome annotation
DomestAnt_oldMAF <- DomestAnt_oldMAF[,c(2:7)]
DomestAnt_oldMAF2 <- reshape(DomestAnt_oldMAF, idvar = c("Chromosome", "Pos", "Index","Gene"), timevar = "Trait", direction = "wide")

#Reformat Chromosomes and Positions

#split chromosome and segment
unique(DomestAnt_oldMAF2$Chromosome)

DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome1$", replacement = "Chromosome1.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome2$", replacement = "Chromosome2.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome3$", replacement = "Chromosome3.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome4$", replacement = "Chromosome4.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome5$", replacement = "Chromosome5.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome6$", replacement = "Chromosome6.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome7$", replacement = "Chromosome7.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome8$", replacement = "Chromosome8.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome9$", replacement = "Chromosome9.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome10$", replacement = "Chromosome10.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome11$", replacement = "Chromosome11.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome12$", replacement = "Chromosome12.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome13$", replacement = "Chromosome13.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome14$", replacement = "Chromosome14.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome15$", replacement = "Chromosome15.0", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2$Chromosome <- gsub(pattern = "Chromosome16$", replacement = "Chromosome16.0", DomestAnt_oldMAF2$Chromosome)
unique(DomestAnt_oldMAF2$Chromosome)

DomestAnt_oldMAF2$Chromosome <- gsub("Chromosome", "", DomestAnt_oldMAF2$Chromosome)
DomestAnt_oldMAF2 <- separate (DomestAnt_oldMAF2, Chromosome, into = c("Chromosome", "Segment") )
#double check
unique(DomestAnt_oldMAF2$Chromosome)
unique(DomestAnt_oldMAF2$Segment)

DomestAnt_oldMAF2$Chromosome <- as.numeric(as.character(DomestAnt_oldMAF2$Chromosome))
DomestAnt_oldMAF2$Segment <- as.numeric(as.character(DomestAnt_oldMAF2$Segment))
DomestAnt_oldMAF2$Pos <- as.numeric(as.character(DomestAnt_oldMAF2$Pos))

#sort dataframe rows in order of Chrom, then Pos
DomestAnt_oldMAF2 <- DomestAnt_oldMAF2[with(DomestAnt_oldMAF2, order(Chromosome, Segment, Pos)), ]

#now make segments line up consecutively
DomestAnt_oldMAF2$Chrom.Seg <- paste(DomestAnt_oldMAF2$Chromosome, DomestAnt_oldMAF2$Segment, sep=".")
DomestAnt_oldMAF2$Chrom.Seg <- as.numeric(DomestAnt_oldMAF2$Chrom.Seg)

#let's try making the chrom.seg integers so that R isn't confused
unique(DomestAnt_oldMAF2$Chrom.Seg)
DomestAnt_oldMAF2$Chrom.Seg.F <- as.factor(DomestAnt_oldMAF2$Chrom.Seg)
unique(DomestAnt_oldMAF2$Chrom.Seg.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6,  6.1, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 13.2, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

DomestAnt_oldMAF2$Chrom.Seg.Int <- recode.vars$newvals[match(DomestAnt_oldMAF2$Chrom.Seg.F, recode.vars$OGvals)]
unique(DomestAnt_oldMAF2$Chrom.Seg.Int)

#Make plotting variables
DomestAnt_oldMAF2$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing

for (i in unique(DomestAnt_oldMAF2$Chrom.Seg.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of DomestAnt_oldMAF2 rows with Chromosome 1, set Index variable for each row to equal Pos.
    DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Index=DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(DomestAnt_oldMAF2,DomestAnt_oldMAF2$Chrom.Seg.Int==i-1)$Pos, 1)
    #and then for the subset of DomestAnt_oldMAF2 rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Index=DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Pos+lastbase
  }
  #set ticks to be a list of existing ticks, plus the current Index
  #floor rounds it down to the nearest whole number
  # ticks=c(ticks, DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Index[floor(length(DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom==i, ]$Index)/2)+1])
  
  ticks=c(ticks, DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Index[floor(length(DomestAnt_oldMAF2[DomestAnt_oldMAF2$Chrom.Seg.Int==i, ]$Index)/2)+1])
}
ticklim=c(min(DomestAnt_oldMAF2$Index),max(DomestAnt_oldMAF2$Index))


write.csv(DomestAnt_oldMAF2, "results/Domestication_TopSNPs_SegWide_annot.csv")
