#Nicole E Soltis
#091317
#plot of SNPs along gene of interest
#and now: haplotype plots

#14_D_myfile.hlist.R
#-------------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#bring Chr16_A.assoc.hap into excel or libreoffice calc. Do fixed-width import -- tab separated didn't work in R directly. Save as .csv
myHlist <-  read.csv("data/genome/chr16_analysis/haps/chr16A_assoc.hap.csv")

myHlist$first <- "*"
myHlist <- myHlist[,c("first","SNPS")]
#only keep one record per 3-SNP haplotype
myHlist <- unique(myHlist)

write.table(myHlist, file="data/genome/chr16_analysis/haps/myfile.hlist", quote=FALSE, sep=' ', col.names=FALSE, row.names=FALSE)
