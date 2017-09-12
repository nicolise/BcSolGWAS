#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/")
#on windows laptop
setwd("~/Projects/BcSolGWAS/")
#convert .tab SNP file to .csv
#sticking to MAF20, 10% missingness

 tab20 = read.delim("GWAS_files/01_tabfiles/Suzi_033016/Haploid_SNPS_97_dp6_maf20.tab")
  write.table(tab20, file="GEMMA_files/02_csvPrep/snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)


#convert each file to binary
SNPsMAF20 <- read.csv("data/GEMMA_files/02_csvPrep/snps_maf20.csv")
mySNPs <- SNPsMAF20

mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="."]<-NA #this is a true NA
allSNPs<- mySNPs
mySNPs_Chr16 <- mySNPs[grep("Chromosome16", mySNPs$X.CHROM), ]

#remove low MAFs!
names(mySNPs_Chr16)

mySNPs_Chr16$Freq.A <- rowSums(mySNPs_Chr16[,c(4:100)] =="A", na.rm=T)
mySNPs_Chr16$Freq.C <- rowSums(mySNPs_Chr16[,c(4:100)] =="C", na.rm=T)
mySNPs_Chr16$Freq.T <- rowSums(mySNPs_Chr16[,c(4:100)] =="T", na.rm=T)
mySNPs_Chr16$Freq.G <- rowSums(mySNPs_Chr16[,c(4:100)] =="G", na.rm=T)
mySNPs_Chr16$All <- (mySNPs_Chr16$Freq.A + mySNPs_Chr16$Freq.C + mySNPs_Chr16$Freq.G + mySNPs_Chr16$Freq.T)
#97 isolates total
#remove all loci (rows) with >10% missingness (mySNPs_Chr16$All <87.3)
mySNPs_Chr16 <- mySNPs_Chr16[mySNPs_Chr16$All>87,] 

#now remove all loci (rows) with MAF <20 (max > 77.6)
mySNPs_Chr16 <- transform(mySNPs_Chr16, max = pmax(Freq.A, Freq.C, Freq.G, Freq.T))
mySNPs_Chr16 <- mySNPs_Chr16[mySNPs_Chr16$max < 78,] #this is fine, keep all.
mySNPs_Chr16$missing <- apply(mySNPs_Chr16[,c(4:100)], 1, function(x) sum(is.na(x)))

#MAF is too low if: max + missing > 78
mySNPs_Chr16 <- mySNPs_Chr16[(mySNPs_Chr16$max + mySNPs_Chr16$missing) < 78,]

#NOW filter to just region of interest on Chromosome 16
mySNPs_Chr16 <- mySNPs_Chr16[which(mySNPs_Chr16$POS < 355042),]
mySNPs_Chr16 <- mySNPs_Chr16[which(mySNPs_Chr16$POS > 337285),]
#save a list of SNPs to match with 14_Fig8b_LDplot.R where I'm writing SNP.FILE
mySNPlist <- as.data.frame(mySNPs_Chr16$POS)

write.csv(mySNPlist, "data/genome/chr16_analysis/SNPlistfromPED.csv")

#and remove any duplicated POS here
mySNPs_Chr16 <- mySNPs_Chr16[!duplicated(mySNPs_Chr16$POS), ]
#excellent, now this matches the SNPs file

#and now for making PED format for PLINK!
  #do not need positional info: just SNP states for PED
#turn df sideways (individuals as rows, SNPs as columns)
#split each genotype into 2 identical columns (PED assumes diploid)
#keeping REF variable for now, to make binary phenotype for GENOTYPE.FILE
mySNPs_Chr16_2 <- mySNPs_Chr16[,-c(1:2, 101:107)]
#write.csv(mySNPs_Chr16_2, "GEMMA_files/02_csvPrep/hp_charMAF20_10NA_forPED.csv")
#df2[,c(1,3,2,4)]

#duplicate all the rows, to fake a diploid genome
#this is slow (but only takes ~4 seconds if 600 SNPs)
Sys.time()
mySNPs_Chr16_3 <- as.data.frame(NULL)
for(i in 1:nrow(mySNPs_Chr16_2)){
  mySNPs_Chr16_3<-rbind(mySNPs_Chr16_3, mySNPs_Chr16_2[i,])
  mySNPs_Chr16_3<-rbind(mySNPs_Chr16_3, mySNPs_Chr16_2[i,])
}
Sys.time()

#write.csv(mySNPs_Chr16_3, "data/genome/chr16_analysis/dp_charMAF20_10NA.csv")
#mySNPs_Chr16_3 <- read.csv("data/genome/chr16_analysis/dp_charMAF20_10NA.csv")

#transpose and format for PED
#mySNPs_Chr16_3 <- mySNPs_Chr16_3[,-c(1)]
mySNPs_Chr16_4 <- as.data.frame(t(mySNPs_Chr16_3))
#add binary phenotype = 1 (6)
mySNPs_Chr16_4 <- cbind("Pheno" = 1, mySNPs_Chr16_4)
#add individual sex = 1 (5)
mySNPs_Chr16_4 <- cbind("sex" = 1, mySNPs_Chr16_4)
#add Mother = 0 (4)
mySNPs_Chr16_4 <- cbind("Mother" = 0, mySNPs_Chr16_4)
#add Father = 0 (3)
mySNPs_Chr16_4 <- cbind("Father" = 0, mySNPs_Chr16_4)
#turn row names into column 2
mySNPs_Chr16_4 <- cbind(rownames(mySNPs_Chr16_4), mySNPs_Chr16_4)
colnames(mySNPs_Chr16_4)[1] <- 'Isolate'
#add the fam column (1)
mySNPs_Chr16_4 <- cbind("FAM" = "FAM1", mySNPs_Chr16_4)
#remove ref for myPED
myPED <- mySNPs_Chr16_4[-c(1),] 

#finally, make the GENOTYPE.FILE for SNPplotter
#transposed and formatted for modified PED
myGENOTYPE.FILE <- mySNPs_Chr16_4
rownames(myGENOTYPE.FILE)
# family ID, individual ID, father ID, mother ID, sex, and affection status followed by marker loci coded as binary factors, as shown in the example below. This file should not have column headers.
#make marker loci binary

mytest <- myGENOTYPE.FILE[1:10,7:10]
mynames <- myGENOTYPE.FILE[1:10,1:6]

for (i in rownames(mytest[2:10,])){
  mytest[i,][mytest[i,]!=mytest[1,]] <- "A"
  mytest[i,][mytest[i,]==mytest[1,]] <- 0
}

#this works as logical:
mytest[2,]!=mytest[1,]
#this works as logical: 
mytest[,1]!=mytest[1,1]

for (i in names(mySNPs[4:100])) {
  mySNPs[i][mySNPs[i]!=mySNPs$REF] <- 1
  mySNPs[i][mySNPs[i]==mySNPs$REF] <- 0
}


#make a MAP file for plink (need it to make the bed (binary ped) file from ped)
myMAP <- mySNPs_Chr16[,c("X.CHROM","POS")]
#remove "Chromosome" from X.CHROM
myMAP <- as.data.frame(lapply(myMAP, function(x) {
                 gsub("Chromosome", "", x)
              }))

#----------------------------------------------------------------------
#save all of the files
write.csv(myMAP, "data/genome/chr16_analysis/MAP_dpcharMAF20NA10.csv")
#add a column of "SNP identifiers" in excel and remove headers

write.csv(mySNPs_Chr16_3, "data/genome/chr16_analysis/dp_charMAF20_10NA.csv")
write.csv(mySNPs_Chr16, "data/genome/chr16_analysis/hp_charMAF20_10NA.csv")
write.csv(mySNPs_Chr16_4, "data/genome/chr16_analysis/PED_dpcharMAF20NA10.csv")
