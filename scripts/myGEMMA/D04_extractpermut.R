#Nicole E Soltis
#04/05/18

#-------------------------------------------------------------------------
rm(list = ls())
#list all files in /media/ outputs
#loop through all 1000 files
# for each of 12 phenotypes, create new dataframe that keeps only selected SNPs:
  #df 1: ordered on SNP p-value top 1% = 2720 SNPs ish
  #df 2: ordered on SNP p-value small -> large, take 1st, 2nd, 25th, 250th, 2500th SNP
#save df 1 to media
#save df 2 to ~/Projects
#loop through all df2 to build a new df across all 1000 permutations using rbind
  #columns: chromosome, position, p-value on phenotype 1...12, which SNP (1 / 2 / 25 / 250 / 2500)
  #
#setwd("~/Projects/BcSolGWAS/data/GEMMA_files/D_05_bigrand")
setwd("/media/nesoltis/Data/Kliebenstein/Soltis/BcSolGWAS/data/GEMMA_files/")

#read in files from each of 1000 rand folders
#each randomization
#for (i in 1:1000){
#done: i = 1:3
#9am tomorrow is 16 hours = 480 phenos
for (i in 4:500){
  newdir <- paste0("D_06_randOUT/topouts/rand1k_",i)
  dir.create(newdir)
  #each phenotype
  for (j in 1:15){
    Sys.time()
    my_gemma <- read.table(paste("D_06_randOUT/output/rand1k_",i,"/pheno",j,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 8 seconds to read 1 phenotype
    #times 15 times 1000 = 33 hours
    ntop1 <- nrow(my_gemma)*0.01
    my_gemma_top1 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top1 <- my_gemma_top1[1:2500,]
    my_gemma_points <- my_gemma_top1[c(1,2,25,250,2500),]
    my_gemma_points$SNPnum <- c(1,2,25,250,2500)
    my_gemma_points$run <- i
    my_gemma_points$pheno <- j
    #this gives an error but it's fine
    try(ifelse( i == 1 & j ==1, write.table(my_gemma_points, "D_07_randSUMM/GEMMA_1krand_SNPsample.csv", sep = ",", col.names = T), write.table(my_gemma_points, "D_07_randSUMM/GEMMA_1krand_SNPsample.csv", sep = ",", col.names = F, append = T)))
    write.csv(my_gemma_top1, paste("D_06_randOUT/topouts/rand1k_",i,"/pheno",j,".csv", sep=""))
    Sys.time()
  }
}

#now cp D_07_randSUMM to ~/Documents/GitRepos/BcSolGWAS/
#and make sure nesoltis user account has rwx permission for the new files

      