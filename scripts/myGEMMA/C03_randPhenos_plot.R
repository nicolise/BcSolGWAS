#Nicole E Soltis
#03/29/18

#------------------------------------------------------------------------------
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")
#10: Pheno order is "V2"           "LA0410"       "LA0480"       "LA1547"       "LA1589"
#"LA1684"       "LA2093"       "LA2176"       "LA2706"       "LA3008"      
#"LA3475"       "LA4345"       "LA4355"       "Domesticated" "Wild"   "DmWoD"
my.LA0410 <- read.table("03_GEMMAouts/randtest/binMAF20NA10_randtest_PLINK_kmat1_pheno1.assoc.txt", header=TRUE)
all.Phenos <- read.csv("03_GEMMAouts/GEMMA_lmm_ind12plants/kmat1/12Plants_allSNPs_MAF20NA10_GEMMA_kmat1.csv")

#combine rand outputs
setwd("~/Projects/BcSolGWAS/data/GEMMA_files/03_GEMMAouts/randtest/kmat1/")
#read in each of 12 phenotypes
my.files <- list.files(pattern = c("assoc"))
#2_LA1589 fails with both k-matrix options.
#didn't run yet for DWS, only 1:12
my.names <- c("1_LA0410", "10_LA3475", "11_LA4345", "12_LA4355", "2_LA0480", "3_LA1547", "4_LA1589", "5_LA1684", "6_LA2093", "7_LA2176", "8_LA2706", "9_LA3008")
              #, "13_Domesticated", "14_Wild", "15_DmWoD")
#rename all files
 for(i in 1:length(my.files)) {
   my.file <- read.csv(my.files[i])
   file.rename(from=file.path(my.files[i]), to=file.path(paste(file.path(my.files[i]),my.names[i],".txt",sep="")))
 }

#now read in files and combine
my.files <- list.files(pattern = c("LA"))
for (i in 1:length(my.files)){
  my.file <- read.table(my.files[i], header=TRUE)
  ifelse(i == 1, full.file <- my.file, full.file <- cbind(full.file, my.file[,c("beta","se","p_score")]))
  ifelse(i == 1, names(full.file)[9] <- paste(my.names[i], "_beta", sep=""), names(full.file)[(ncol(full.file)-2)] <- paste(my.names[i], "_beta", sep=""))
  ifelse(i == 1, names(full.file)[10] <- paste(my.names[i], "_se", sep=""), names(full.file)[(ncol(full.file)-1)] <- paste(my.names[i], "_se", sep=""))
  ifelse(i == 1, names(full.file)[13] <- paste(my.names[i], "_pscore", sep=""), names(full.file)[(ncol(full.file))] <- paste(my.names[i], "_pscore", sep=""))
}
full.file <- full.file[,-c(2,11,12)]
pheno.bin <- full.file
#now, PosPhenos = 1 if p < 0.01, 0 if p > 0.01
for (i in c(1:12)){
  mybeta = 5 + (3*i)
  mypscore = 7 + (3*i)
  pheno.bin[,paste(colnames(pheno.bin[mybeta]),"bin", sep="")] <- ifelse(pheno.bin[,mypscore] < 0.01, 1, 0)
}

names(pheno.bin)
pheno.bin$SUMM <- rowSums(pheno.bin[,c(44:55)], na.rm=T)

rand.Phenos <- pheno.bin
all.Phenos <- all.Phenos[,-c(1)]

names(all.Phenos)[56] <- "X.SUMM"
compare.Phenos <- cbind(all.Phenos, rand.Phenos[,c(8:55)])

#reorder columns cause this is a mess
compare.Phenos <- compare.Phenos[,c("chr", "ps", "n_mis", "n_obs", "allele1", "allele0", "af", "X9_LA0410_beta", "X9_LA0410_se", "X9_LA0410_pscore","X12_LA0480_beta", "X12_LA0480_se", "X12_LA0480_pscore", "X1_LA1547_beta", "X1_LA1547_se", "X1_LA1547_pscore", "X2_LA1589_beta", "X2_LA1589_se", "X2_LA1589_pscore", "X3_LA1684_beta", "X3_LA1684_se", "X3_LA1684_pscore", "X4_LA2093_beta", "X4_LA2093_se", "X4_LA2093_pscore", "X5_LA2176_beta", "X5_LA2176_se", "X5_LA2176_pscore", "X6_LA2706_beta", "X6_LA2706_se", "X6_LA2706_pscore", "X7_LA3008_beta", "X7_LA3008_se", "X7_LA3008_pscore", "X8_LA3475_beta", "X8_LA3475_se", "X8_LA3475_pscore",  "X10_LA4345_beta", "X10_LA4345_se", "X10_LA4345_pscore", "X11_LA4355_beta", "X11_LA4355_se", "X11_LA4355_pscore", "1_LA0410_beta", "1_LA0410_se", "1_LA0410_pscore", "2_LA0480_beta", "2_LA0480_se", "2_LA0480_pscore", "3_LA1547_beta", "3_LA1547_se", "3_LA1547_pscore", "4_LA1589_beta", "4_LA1589_se", "4_LA1589_pscore", "5_LA1684_beta", "5_LA1684_se", "5_LA1684_pscore", "6_LA2093_beta", "6_LA2093_se", "6_LA2093_pscore", "7_LA2176_beta", "7_LA2176_se", "7_LA2176_pscore", "8_LA2706_beta", "8_LA2706_se", "8_LA2706_pscore", "9_LA3008_beta", "9_LA3008_se", "9_LA3008_pscore", "10_LA3475_beta", "10_LA3475_se", "10_LA3475_pscore", "11_LA4345_beta", "11_LA4345_se", "11_LA4345_pscore", "12_LA4355_beta", "12_LA4355_se", "12_LA4355_pscore", "X1_LA1547_betabin", "X2_LA1589_betabin", "X3_LA1684_betabin", "X4_LA2093_betabin", "X5_LA2176_betabin", "X6_LA2706_betabin", "X7_LA3008_betabin", "X8_LA3475_betabin", "X9_LA0410_betabin", "X10_LA4345_betabin", "X11_LA4355_betabin", "X12_LA0480_betabin",  "X.SUMM", "1_LA0410_betabin", "10_LA3475_betabin", "11_LA4345_betabin", "12_LA4355_betabin", "2_LA0480_betabin", "3_LA1547_betabin", "4_LA1589_betabin", "5_LA1684_betabin", "6_LA2093_betabin", "7_LA2176_betabin", "8_LA2706_betabin", "9_LA3008_betabin")]


#no p-values for 4 (failed)
#observed phenotypes are X1 etc. format -> column 1 -> x on qq plot
#randomized phenos are 1 etc. format -> column 2 -> y on qq plot
setwd("~/Projects/BcSolGWAS/data/GEMMA_files/04_analysis")
for (j in c(1:3,5:12)){
  mypheno <- 7 + 3*j
  myrand <- 43 + 3*j
  #hist(compare.Phenos[,mypheno])
  #hist(compare.Phenos[,myrand])
  mypheno.var <- as.data.frame(compare.Phenos[,c(mypheno, myrand)])
  jpeg(paste("randQQ/qqplot", names(mypheno.var)[1], ".jpg", sep=""), width=7.5, height=5, units='in', res=600)
  qqplot(mypheno.var[,1], mypheno.var[,2])
  dev.off()
}
#   mypheno.var <- mypheno.var[order(mypheno.var[,1]),]
#   mypheno.var$obsrank <- 1:nrow(mypheno.var)
#   mypheno.var <- mypheno.var[order(mypheno.var[,2]),]
#   mypheno.var$randrank <- 1:nrow(mypheno.var)
#   #blah <- quantile(compare.Phenos[,myrand], probs = seq(0, 0.01, by= 0.001)) 
#   #print(names(compare.Phenos)[myrand])
#   #print(blah)
#   plot(mypheno.var$obsrank ~ mypheno.var$randrank)
# }

#top 1%: 
mean(0.012, 0.014, 0.019, 0.007, 0.013, 0.016, 0.010, 0.016, 0.018, 0.02)
#top 0.1%: 
mean(0.004, 0.003, 0.001, 0.002, 0.0008, 0.004, 0.003, 0.002)
