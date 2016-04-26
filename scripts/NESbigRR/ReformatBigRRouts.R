#reformat bigRR output data
#Nicole E Soltis

#--------------------------------------------------------
rm(list=ls())
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcSolGWAS/data/bigRRout")
#Import data
HEMdat <- read.csv("Sl_LesionSizePoisson_2.HEM.csv")
#RF-this is just a reorganized file of LesionSize.HEM.csv.
#RF-Basically it has two columns containing Chrom and pos separately instead of just one column with eg "III.57894". This is easiest to do with code below:
library(tidyr)
names(HEMdat)
HEMdat2 <- separate (HEMdat, X.1, into = c("Chrom", "Pos") )
#library(tidyr)
#separate(old data, V1, into = c("Chrom", "Pos"))

write.csv(HEMdat2, "Sl_LesionSizePoisson.HEM.PlotFormat.csv") 
#read in to RunningBigRRwithRF.R, line 172