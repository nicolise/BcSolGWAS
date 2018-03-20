#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
#first round: just one file at a time. Then convert to loops
#1 is Domest, 2 is Wild, 3 is Sensitivity
#original files, no accounting for pop str
#myGEMMA.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_1.assoc.txt", header=TRUE)
#myGEMMA.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_2.assoc.txt", header=TRUE)
#myGEMMA.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_3.assoc.txt", header=TRUE)
#based on Manhattan plots, kmat1 and kmat2 are identical
myGEMMA.2.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_1.assoc.txt", header=TRUE)
myGEMMA.2.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_2.assoc.txt", header=TRUE)
myGEMMA.2.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_3.assoc.txt", header=TRUE)

myGEMMA.1.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_1.assoc.txt", header=TRUE)
myGEMMA.1.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_2.assoc.txt", header=TRUE)
myGEMMA.1.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_3.assoc.txt", header=TRUE)

#run through this once for myGEMMA.2, once for myGEMMA.1
names(myGEMMA.1.D)
myGEMMA <- myGEMMA.1.D[,c(1,3,9,13)]
names(myGEMMA)[3] <- "beta.D"
names(myGEMMA)[4] <- "pscore.D"
myGEMMA <- cbind(myGEMMA,myGEMMA.1.W[,c(9,13)])
names(myGEMMA)[5] <- "beta.W"
names(myGEMMA)[6] <- "pscore.W"
myGEMMA <- cbind(myGEMMA,myGEMMA.1.S[,c(9,13)])
names(myGEMMA)[7] <- "beta.S"
names(myGEMMA)[8] <- "pscore.S"
#Make plotting variables
myGEMMA$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(myGEMMA$chr)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    myGEMMA[myGEMMA$chr==i, ]$Index=myGEMMA[myGEMMA$chr==i, ]$ps
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(myGEMMA,myGEMMA$chr==i-1)$ps, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    myGEMMA[myGEMMA$chr==i, ]$Index=myGEMMA[myGEMMA$chr==i, ]$ps+lastbase
  }
}

#add a tottraits variable
myGEMMA$TotTraits <- ifelse(myGEMMA$pscore.D < 0.01 & myGEMMA$pscore.W < 0.01 & myGEMMA$pscore.S <0.01, "ALL",
                              ifelse(myGEMMA$pscore.D < 0.01 & myGEMMA$pscore.W < 0.01, "DW",
                                     ifelse(myGEMMA$pscore.W < 0.01 & myGEMMA$pscore.S <0.01, "WS",
                                            ifelse(myGEMMA$pscore.D < 0.01 & myGEMMA$pscore.S <0.01, "DS",
                                                   ifelse(myGEMMA$pscore.W < 0.01, "W", 
                                                          ifelse(myGEMMA$pscore.S < 0.01,"S", 
                                                                 ifelse(myGEMMA$pscore.D < 0.01, "D", "none")))))))

table(myGEMMA$TotTraits)

myGEMMA.fulldat <- myGEMMA
write.csv(myGEMMA.fulldat, "data/GEMMA_files/04_analysis/GEMMA_allDWS_kmat1.csv")
#select just top SNPs for comparison to bigRR T4
#conditionally replace nonsig values with zero
hist(myGEMMA$beta.D)
myGEMMA$beta.D[myGEMMA$pscore.D > 0.01] <- 0
myGEMMA$beta.W[myGEMMA$pscore.W > 0.01] <- 0
myGEMMA$beta.S[myGEMMA$pscore.S > 0.01] <- 0
#remove rows if all 3 = 0 
myGEMMA_2 <- myGEMMA[!(myGEMMA$beta.D==0 & myGEMMA$beta.W==0 & myGEMMA$beta.S==0),]
#now add counting variable
myGEMMA_2$TotTraits <- ifelse(myGEMMA_2$beta.D != 0 & myGEMMA_2$beta.W != 0 & myGEMMA_2$beta.S != 0, "ALL",
                              ifelse(myGEMMA_2$beta.D != 0 & myGEMMA_2$beta.W != 0, "DW",
                                     ifelse(myGEMMA_2$beta.W != 0 & myGEMMA_2$beta.S != 0, "WS",
                                            ifelse(myGEMMA_2$beta.D != 0 & myGEMMA_2$beta.S != 0, "DS",
                                                   ifelse(myGEMMA_2$beta.D != 0, "D",
                                                          ifelse(myGEMMA_2$beta.W != 0, "W", "S"))))))

table(myGEMMA_2$TotTraits)

write.csv(myGEMMA_2, "data/GEMMA_files/04_analysis/GEMMA_peaksDWS_kmat1.csv")