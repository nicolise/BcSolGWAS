#Nicole E Soltis
#041118
#D08_violin_permut
#-------------------------------------------------------------------
#draw violin plots from permutation data
rm(list=ls())
library("ggplot2")
setwd("~/Documents/GitRepos/BcSolGWAS/data/GEMMA_files")
mydat <- read.csv("D_07_randSUMM/GEMMA_1krand_SNPsample.csv")

#split mydat into individual phenotypes - each a unique plot
#
eachpheno <- split( mydat , f = mydat$pheno )

for ( i in 1:15 ){
i <- 1
  plotdat <- eachpheno[[i]]
  plotdat$SNPset <- as.factor(plotdat$SNPnum)
  p <- ggplot(plotdat, aes(x=SNPset, y=p_score)) + 
    geom_violin(trim=TRUE)+ #trim tails
    scale_y_continuous(limits=c(0, 0.01))
  p
}