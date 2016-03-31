#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
#on laptop
# setwd("~/Projects/BcSolGWAS/data/genome")
#practice with mini file first
#miniSNPs <- read.csv("minisnps.csv")
#miniSNPs <- read.csv("miniMAF20.csv")
#mySNPs <- miniSNPs

#convert both all .tab SNP files to .csv
# tab5 = read.delim("snps_maf5.tab")
#  write.table(tab5, file="snps_maf5.csv",sep=",",col.names=T,row.names=FALSE)
# tab10 = read.delim("snps_maf10.tab")
#  write.table(tab10, file="snps_maf10.csv",sep=",",col.names=T,row.names=FALSE)
# tab20 = read.delim("snps_maf20.tab")
#  write.table(tab20, file="snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)
# tab20 = read.delim("Haploid_SNPS_97_dp6_maf20.tab")
#  write.table(tab20, file="hp_snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)


#convert each file to binary
# SNPsMAF5 <- read.csv("snps_maf5.csv")
# SNPsMAF10 <- read.csv("snps_maf10.csv")
# SNPsMAF20 <- read.csv("snps_maf20.csv")
SNPsMAF20hp <- read.csv("hp_snps_maf20.csv")
mySNPs <- SNPsMAF20hp

#17 possible states: A/A G/G T/T C/C G/A T/C A/T C/T T/A A/G G/T C/A A/C ./. G/C T/G C/G
#how to deal with heterozygosity?
#convert N/N format to N if diploid calls
# str(mySNPs)
#make these characters instead of factors
mySNPs[] <- lapply(mySNPs, as.character)
# mySNPs[mySNPs=="C/C"]<-"C"
# mySNPs[mySNPs=="T/T"]<-"T"
# mySNPs[mySNPs=="A/A"]<-"A"
# mySNPs[mySNPs=="G/G"]<-"G"
# mySNPs[mySNPs=="A/C"]<-"NA"
# mySNPs[mySNPs=="A/G"]<-"NA"
# mySNPs[mySNPs=="A/T"]<-"NA"
# mySNPs[mySNPs=="C/A"]<-"NA"
# mySNPs[mySNPs=="C/G"]<-"NA"
# mySNPs[mySNPs=="C/T"]<-"NA"
# mySNPs[mySNPs=="G/A"]<-"NA"
# mySNPs[mySNPs=="G/C"]<-"NA"
# mySNPs[mySNPs=="G/T"]<-"NA"
# mySNPs[mySNPs=="T/A"]<-"NA"
# mySNPs[mySNPs=="T/C"]<-"NA"
# mySNPs[mySNPs=="T/G"]<-"NA"
# mySNPs[mySNPs=="./."]<-"NA"
mySNPs[mySNPs=="."]<-"NA"
allSNPs<- mySNPs

#unlist(unique(mySNPs$X305))
#table(mySNPs$X305)

#mySNPsMAF20 <- mySNPs
#miniMAF20 <- mySNPsMAF20[1:20,1:20]
#mySNPs <- miniMAF20

#replace base with 1 if match
#loop through it
#automatically skips NAs
mySNPs <- allSNPs
for (i in names(mySNPs[4:100])) {
  mySNPs[i][mySNPs[i]!=mySNPs$REF]  <- 1
  mySNPs[i][mySNPs[i]==mySNPs$REF] <- 0
}

write.csv(allSNPs, "hp_charMAF20.csv")
write.csv(mySNPs, "hp_binaryMAF20.csv")
#-------------------------------------------------------------
#junk code
#other ideas
within(miniSNPs, X305 <- 1*(X305 == REF))
#alternate way
miniSNPs$X305b[miniSNPs$X305!=miniSNPs$REF] <- 1
miniSNPs$X305b[miniSNPs$X305==miniSNPs$REF] <- 0
#work in progress
apply( df[,grep("abc", colnames(df))] , 2 , f )
apply(within(miniSNPs, colnames(miniSNPs)))
testing <- within(miniSNPs, newvar <- 1*(i == REF))
mytest <- with(miniSNPs, ifelse(i == REF, 0, 1))
tapply(summaryvar, groupvar, function)

#doesn't work
is.na
#if else
for (i in names(mySNPs[4:18])) {
  if (mySNPs[i]==NA){
    [mySNPs[i] <- NA}else 
      if (mySNPs[i]!=mySNPs$REF) {
        mySNPs[i]<- 1}else 
          if (mySNPs[i]==mySNPs$REF) {
            mySNPs[i]<- 0}
}
