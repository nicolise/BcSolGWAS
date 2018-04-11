#Nicole E Soltis
#041118
#D08_violin_permut
#-------------------------------------------------------------------
#draw violin plots from permutation data
rm(list=ls())
library("ggplot2")
setwd("~/Documents/GitRepos/BcSolGWAS/data/GEMMA_files")
mydat <- read.csv("D_07_randSUMM/GEMMA_1krand_SNPsample.csv")

setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#replace phenotype numbers with names
Phenos_match <- read.csv("data/GEMMA_files/D_02_randGEMMA/binMAF20NA10_fam.csv")
names(Phenos_match)

phenos_list <- names(Phenos_match)[7:21]
for (j in 1:15){
mydat$pheno[mydat$pheno == j ] <- phenos_list[j]
}

#split mydat into individual phenotypes - each a unique plot
#
eachpheno <- split( mydat , f = mydat$pheno )

for ( i in 1:15 ){
  plotdat <- eachpheno[[i]]
  plotdat$SNPset <- as.factor(plotdat$SNPnum)
  p <- ggplot(plotdat, aes(x=SNPset, y=-log10(p_score))) + 
    geom_violin(trim=TRUE, fill="grey")+ #trim tails
    theme_bw()+
    labs(list(y=expression('-log'[10]*'p'), title=plotdat$pheno))
    #scale_y_continuous(limits=c(0, 0.01))
  jpeg(paste0("paper/plots/addGEMMA/randpermuts/RandThr_", unique(plotdat$pheno), ".jpg"), width=8, height=5, units='in', res=600)
  print(p + geom_boxplot(width=0.1))
  dev.off()
}

#chosen Threshold: mean of 25th SNP
mythrmeans <- NA
for ( i in 1:15 ){
  plotdat <- eachpheno[[i]]
  plotdat$SNPset <- as.factor(plotdat$SNPnum)
  mypval <- plotdat[,c("SNPnum","p_score")]
mynewthr <- as.data.frame(aggregate(.~SNPnum, data=mypval, mean))
mynewthr$pheno <- unique(plotdat$pheno)
ifelse( i == 1, mythrmeans <- mynewthr, mythrmeans <- rbind(mythrmeans, mynewthr))
}

write.csv(mythrmeans, "data/GEMMA_files/D_07_randSUMM/GEMMA_1krand_thresholds.csv")
