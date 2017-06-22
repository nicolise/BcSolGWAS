#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/")

#convert .tab SNP file to .csv
#sticking to MAF20, 10% missingness

 tab20 = read.delim("GWAS_files/01_tabfiles/Suzi_033016/Haploid_SNPS_97_dp6_maf20.tab")
  write.table(tab20, file="GEMMA_files/02_csvPrep/snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)


#convert each file to binary
SNPsMAF20 <- read.csv("GEMMA_files/02_csvPrep/snps_maf20.csv")
mySNPs <- SNPsMAF20

mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="."]<-NA #this is a true NA
allSNPs<- mySNPs

#remove low MAFs!
names(mySNPs)

mySNPs$Freq.A <- rowSums(mySNPs[,c(4:100)] =="A", na.rm=T)
mySNPs$Freq.C <- rowSums(mySNPs[,c(4:100)] =="C", na.rm=T)
mySNPs$Freq.T <- rowSums(mySNPs[,c(4:100)] =="T", na.rm=T)
mySNPs$Freq.G <- rowSums(mySNPs[,c(4:100)] =="G", na.rm=T)
mySNPs$All <- (mySNPs$Freq.A + mySNPs$Freq.C + mySNPs$Freq.G + mySNPs$Freq.T)
#97 isolates total
#remove all loci (rows) with >10% missingness (mySNPs$All <87.3)
mySNPs <- mySNPs[mySNPs$All>87,] 

#now remove all loci (rows) with MAF <20 (max > 77.6)
mySNPs <- transform(mySNPs, max = pmax(Freq.A, Freq.C, Freq.G, Freq.T))
mySNPs <- mySNPs[mySNPs$max < 78,] #this is fine, keep all.
mySNPs$missing <- apply(mySNPs[,c(4:100)], 1, function(x) sum(is.na(x)))

#MAF is too low if: max + missing > 78
mySNPs <- mySNPs[(mySNPs$max + mySNPs$missing) < 78,]

#and now for making PED format for PLINK!

write.csv(mySNPs, "GEMMA_files/02_csvPrep/hp_charMAF20_10NA.csv")
