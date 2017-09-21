#Nicole E Soltis
#091317
#plot of SNPs along gene of interest
#and now: haplotype plots

#14_B_bimForHaploView
#---------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
myBIM <- read.table("data/genome/chr16_analysis/plink/myCHR16_A.bim", sep = "\t")
names(myBIM)[1] <- "Chr"
names(myBIM)[2] <- "SNP"
names(myBIM)[3] <- "Cm"
names(myBIM)[4] <- "BP"
names(myBIM)[5] <- "minAll"
names(myBIM)[6] <- "majAll"

write.table(myBIM, file="data/genome/chr16_analysis/plink/myCHR16_A_bim.txt",quote=FALSE, sep=' ',col.names=TRUE,row.names=FALSE)