#Nicole E Soltis
#05/07/18

#-----------------------------------------------------------------------------------------------
#calculate pairwise distances in R between tomato accessions based on SNP data
#https://popgen.nescent.org/2015-05-18-Dist-SNP.html
#non-evolutionary genetic distances
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
mydata <- read.csv("SlGenome/Sanger/solcap_tomato_opa_genotype_NES.csv")

library("poppr")
library("pegas")

#create a matrix with only genotype data
mygenos <- mydata[,c(4:51)]
ind <- as.character(mydata$Variety) # individual labels 
population <- rep(1,96) # dummy population labels-- can add this into myGenind or not
myGenind <- df2genind(mygenos, ploidy = 1, ind.names = ind, sep="")
myGenind2 <- genind2loci(myGenind)

#Euclidean distance in {adegenet}
distgenEUCL <- dist(myGenind, method = "euclidean", diag = FALSE, upper = FALSE, p = 1)
hist(distgenEUCL)
blah <- data.matrix(distgenEUCL)
setwd("~/Projects/BcSolGWAS")
jpeg("paper/plots/addGEMMA/SolDistance.jpg", width=10, height=10, units='in', res=600)
heatmap(blah)
dev.off()
