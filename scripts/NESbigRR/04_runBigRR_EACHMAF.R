#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/GWAS_files/")

#########################
# This makes the bigRR_update run through the GPU
## RAN first time WITH POISSON. Lesion size expected to be Gaussian
bigRR_update <- function (obj, Z, family = gaussian(link = identity), tol.err = 1e-06, 
                          tol.conv = 1e-08) 
{
  w <- as.numeric(obj$u^2/(1 - obj$leverage))
  w[w < tol.err] <- tol.err
  bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE )
}

library(bigRR) #check if version is 1.3-9

#Get genotype data
SNPs <- read.csv("03_bigRRinput/binSNP_bigRR_MAF20hp.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Base <- do.call(paste, c(SNPs[c("X.CHROM","POS")], sep="."))

rownames(SNPs) <- SNPs[,96] #set the new column of chrom.base as rownames - this could maybe be written as: rownames(SNPs) <- SNPs$Chr.Base?
SNPs <- SNPs[,4:95] #take out first three cols (X.CHROM, POS, REF) and new last col (Chr.Base)

#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}

#read in phenotype data
Phenos <- read.csv("03_bigRRinput/Sl_Pheno_bigRR_MAF20.csv", row.names = 1)
dat <- as.data.frame((Phenos[,2:13]))  #INSERT PHENOTYPE COLUMNS HERE

outpt.HEM <- colnames(SNPs)
thresh.HEM <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)

#Calculate HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { #i will be each isolate
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1)
  Pheno.BLUP.result <- bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE)
  Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  #Permute Thresholds for Phenos - this is what takes forever
  perm.u.HEM <- vector()
  for(p in 1:1000) {  
    if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
    try(temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE), silent=T)
    try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE),silent = TRUE)
    temp.HEM <- bigRR_update(temp.BLUP, SNPs) #REF change- was bigRR_update(Pheno.BLUP.result...
    perm.u.HEM <- c(perm.u.HEM, temp.HEM$u)
  }
  #write.csv(perm.u.HEM, paste("PermEffects_",colnames(dat)[i],".csv",sep=""))
  thresh.HEM$"0.95Thresh"[i] <- quantile(abs(perm.u.HEM),0.95)
  thresh.HEM$"0.975Thresh"[i] <- quantile(abs(perm.u.HEM),0.975)
  thresh.HEM$"0.99Thresh"[i] <- quantile(abs(perm.u.HEM),0.99)
  thresh.HEM$"0.999Thresh"[i] <- quantile(abs(perm.u.HEM),0.999)
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
}
#Give column names to the thresholds from the HEM list
for(j in 1:length(thresh.HEM)) {
  names(thresh.HEM[[j]]) <- colnames(dat)
}

#RF-give row names to thresh.HEM so that threshhold values will line up correctly with phenotypes, and you can see which threshold value is displayed
thresh.HEM$"0.95Thresh" <- c("0.95 Thresh", thresh.HEM$"0.95Thresh")
thresh.HEM$"0.975Thresh" <- c("0.975 Thresh", thresh.HEM$"0.975Thresh")
thresh.HEM$"0.99Thresh" <- c("0.99 Thresh", thresh.HEM$"0.99Thresh")
thresh.HEM$"0.999Thresh" <- c("0.999 Thresh", thresh.HEM$"0.999Thresh")

#Write results to output
write.csv(rbind(thresh.HEM$"0.95Thresh",thresh.HEM$"0.975Thresh",thresh.HEM$"0.99Thresh",thresh.HEM$"0.999Thresh",outpt.HEM),"04_bigRRoutput/Sl_LesionSize_MAF20.HEM.csv")

###########################################------------------------------------
#MAF10

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/GWAS_files/")

#########################
# This makes the bigRR_update run through the GPU
## RAN first time WITH POISSON. Lesion size expected to be Gaussian
bigRR_update <- function (obj, Z, family = gaussian(link = identity), tol.err = 1e-06, 
                          tol.conv = 1e-08) 
{
  w <- as.numeric(obj$u^2/(1 - obj$leverage))
  w[w < tol.err] <- tol.err
  bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE )
}

library(bigRR) #check if version is 1.3-9

#Get genotype data
SNPs <- read.csv("03_bigRRinput/binSNP_bigRR_MAF10hp.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Base <- do.call(paste, c(SNPs[c("X.CHROM","POS")], sep="."))

rownames(SNPs) <- SNPs[,96] #set the new column of chrom.base as rownames - this could maybe be written as: rownames(SNPs) <- SNPs$Chr.Base?
SNPs <- SNPs[,4:95] #take out first three cols (X.CHROM, POS, REF) and new last col (Chr.Base)

#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}

#read in phenotype data
Phenos <- read.csv("03_bigRRinput/Sl_Pheno_bigRR_MAF10.csv", row.names = 1)
dat <- as.data.frame((Phenos[,2:13]))  #INSERT PHENOTYPE COLUMNS HERE

outpt.HEM <- colnames(SNPs)
thresh.HEM <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)

#Calculate HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { #i will be each isolate
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1)
  Pheno.BLUP.result <- bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE)
  Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  #Permute Thresholds for Phenos - this is what takes forever
  perm.u.HEM <- vector()
  for(p in 1:1000) {  
    if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
    temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE)
    try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE),silent = TRUE)
    temp.HEM <- bigRR_update(temp.BLUP, SNPs) #REF change- was bigRR_update(Pheno.BLUP.result...
    perm.u.HEM <- c(perm.u.HEM, temp.HEM$u)
  }
  #write.csv(perm.u.HEM, paste("PermEffects_",colnames(dat)[i],".csv",sep=""))
  thresh.HEM$"0.95Thresh"[i] <- quantile(abs(perm.u.HEM),0.95)
  thresh.HEM$"0.975Thresh"[i] <- quantile(abs(perm.u.HEM),0.975)
  thresh.HEM$"0.99Thresh"[i] <- quantile(abs(perm.u.HEM),0.99)
  thresh.HEM$"0.999Thresh"[i] <- quantile(abs(perm.u.HEM),0.999)
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
}
#Give column names to the thresholds from the HEM list
for(j in 1:length(thresh.HEM)) {
  names(thresh.HEM[[j]]) <- colnames(dat)
}

#RF-give row names to thresh.HEM so that threshhold values will line up correctly with phenotypes, and you can see which threshold value is displayed
thresh.HEM$"0.95Thresh" <- c("0.95 Thresh", thresh.HEM$"0.95Thresh")
thresh.HEM$"0.975Thresh" <- c("0.975 Thresh", thresh.HEM$"0.975Thresh")
thresh.HEM$"0.99Thresh" <- c("0.99 Thresh", thresh.HEM$"0.99Thresh")
thresh.HEM$"0.999Thresh" <- c("0.999 Thresh", thresh.HEM$"0.999Thresh")

#Write results to output
write.csv(rbind(thresh.HEM$"0.95Thresh",thresh.HEM$"0.975Thresh",thresh.HEM$"0.99Thresh",thresh.HEM$"0.999Thresh",outpt.HEM),"04_bigRRoutput/Sl_LesionSize_MAF10.HEM.csv")

###########################################------------------------------------
#MAF5

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/data/GWAS_files/")

#########################
# This makes the bigRR_update run through the GPU
## RAN first time WITH POISSON. Lesion size expected to be Gaussian
bigRR_update <- function (obj, Z, family = gaussian(link = identity), tol.err = 1e-06, 
                          tol.conv = 1e-08) 
{
  w <- as.numeric(obj$u^2/(1 - obj$leverage))
  w[w < tol.err] <- tol.err
  bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE )
}

library(bigRR) #check if version is 1.3-9

#Get genotype data
SNPs <- read.csv("03_bigRRinput/binSNP_bigRR_MAF5hp.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Base <- do.call(paste, c(SNPs[c("X.CHROM","POS")], sep="."))

rownames(SNPs) <- SNPs[,96] #set the new column of chrom.base as rownames - this could maybe be written as: rownames(SNPs) <- SNPs$Chr.Base?
SNPs <- SNPs[,4:95] #take out first three cols (X.CHROM, POS, REF) and new last col (Chr.Base)

#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}

#read in phenotype data
Phenos <- read.csv("03_bigRRinput/Sl_Pheno_bigRR_MAF5.csv", row.names = 1)
dat <- as.data.frame((Phenos[,2:13]))  #INSERT PHENOTYPE COLUMNS HERE

outpt.HEM <- colnames(SNPs)
thresh.HEM <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)

#Calculate HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { #i will be each isolate
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1)
  Pheno.BLUP.result <- bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE)
  Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  #Permute Thresholds for Phenos - this is what takes forever
  perm.u.HEM <- vector()
  for(p in 1:1000) {  
    if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
    temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE)
    try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE),silent = TRUE)
    temp.HEM <- bigRR_update(temp.BLUP, SNPs) #REF change- was bigRR_update(Pheno.BLUP.result...
    perm.u.HEM <- c(perm.u.HEM, temp.HEM$u)
  }
  #write.csv(perm.u.HEM, paste("PermEffects_",colnames(dat)[i],".csv",sep=""))
  thresh.HEM$"0.95Thresh"[i] <- quantile(abs(perm.u.HEM),0.95)
  thresh.HEM$"0.975Thresh"[i] <- quantile(abs(perm.u.HEM),0.975)
  thresh.HEM$"0.99Thresh"[i] <- quantile(abs(perm.u.HEM),0.99)
  thresh.HEM$"0.999Thresh"[i] <- quantile(abs(perm.u.HEM),0.999)
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
}
#Give column names to the thresholds from the HEM list
for(j in 1:length(thresh.HEM)) {
  names(thresh.HEM[[j]]) <- colnames(dat)
}

#RF-give row names to thresh.HEM so that threshhold values will line up correctly with phenotypes, and you can see which threshold value is displayed
thresh.HEM$"0.95Thresh" <- c("0.95 Thresh", thresh.HEM$"0.95Thresh")
thresh.HEM$"0.975Thresh" <- c("0.975 Thresh", thresh.HEM$"0.975Thresh")
thresh.HEM$"0.99Thresh" <- c("0.99 Thresh", thresh.HEM$"0.99Thresh")
thresh.HEM$"0.999Thresh" <- c("0.999 Thresh", thresh.HEM$"0.999Thresh")

#Write results to output
write.csv(rbind(thresh.HEM$"0.95Thresh",thresh.HEM$"0.975Thresh",thresh.HEM$"0.99Thresh",thresh.HEM$"0.999Thresh",outpt.HEM),"Sl_LesionSize_MAF5.HEM.csv")