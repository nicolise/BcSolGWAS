#Nicole E Soltis
#convert .tab SNP data to binary .csv

#------------------------------------------------------
setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
#practice with mini file first
#miniSNPs <- read.csv("minisnps.csv")

#convert both .tab SNP files to .csv
#tab10 = read.delim("snps_maf10.tab")
#write.table(tab10, file="snps_maf10.csv",sep=",",col.names=T,row.names=FALSE)

#tab5 = read.delim("snps_maf5.tab")
#write.table(tab5, file="snps_maf5.csv",sep=",",col.names=T,row.names=FALSE)

#tab20 = read.delim("snps_maf20.tab")
#write.table(tab20, file="snps_maf20.csv",sep=",",col.names=T,row.names=FALSE)

SNPsMAF20 <- read.csv("snps_maf20.csv")
mySNPs <- SNPsMAF20
#convert N/N format to N
str(mySNPs)
#make these characters instead of factors
mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="C/C"]<-"C"
mySNPs[mySNPs=="T/T"]<-"T"
mySNPs[mySNPs=="A/A"]<-"A"
mySNPs[mySNPs=="G/G"]<-"G"
mySNPs[mySNPs=="./."]<-"NA"

#17 possible states: A/A G/G T/T C/C G/A T/C A/T C/T T/A A/G G/T C/A A/C ./. G/C T/G C/G
#how to deal with heterozygosity?

#replace base with 1 if match
#loop through it
for (i in names(mySNPs[4:18])) {
  mySNPs[i][mySNPs[i]!=mySNPs$REF] <- 1
  mySNPs[i][mySNPs[i]==mySNPs$REF] <- 0
}

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
