#Nicole E Soltis
#042517

#-------------------------------------------------------------
rm(list=ls())
library(tidyr); library(plyr)
#setwd("~/Documents/GitRepos/BcSolGWAS/data/SNP_files")
setwd("~/Projects/BcSolGWAS/data/GWAS_files/04_bigRRoutput/trueMAF_10NA/")
#Import data
#reorganize file Sl_LesionSize.HEM.csv
HEMdat <- read.csv("SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.fix.csv")
LA2176.plotdata <- read.csv("SlBc_LA2176_trueMAF20_10NA.HEM.PlotFormat.csv")
LA0480.plotdata <- read.csv("SlBc_LA0480_c_trueMAF20_10NA.HEM.PlotFormat.csv")
LA2176.plotdata <- LA2176.plotdata[,-c(1,6,7,8,9)]
LA0480.plotdata <- LA0480.plotdata[,-c(1,6,7,8,9)]
HEMdat <- HEMdat[,-c(1, 17, 18, 19)]
names(HEMdat)
LA2176.plotdata$Chrom.Seg.Pos <- paste(LA2176.plotdata$Chrom, LA2176.plotdata$Segment, LA2176.plotdata$Pos, sep=".")
LA0480.plotdata$Chrom.Seg.Pos <- paste(LA0480.plotdata$Chrom, LA0480.plotdata$Segment, LA0480.plotdata$Pos, sep=".")
HEMdat$Chrom.Seg.Pos <- paste(HEMdat$Chrom, HEMdat$Segment, HEMdat$Pos, sep=".")
LA2176.plotdata <- LA2176.plotdata[,c(4,5)]
LA0480.plotdata <- LA0480.plotdata[,c(4,5)]
HEMdat2 <- merge(HEMdat, LA2176.plotdata, by="Chrom.Seg.Pos")
HEMdat2 <- merge(HEMdat2, LA0480.plotdata, by="Chrom.Seg.Pos")
HEMdat2 <- HEMdat2[,-c(1,9, 16)]
HEMdat2 <- rename(HEMdat2, c("X.Phenos...8.."="LA2176", "Phenos...13."= "LA0480"))

write.csv(HEMdat2, "SlBc_12plants_trueMAF20_10NA.HEM.PlotFormat.final.csv")

#and merge Thresholds
HEMthresh <- read.csv("SlBc_12plants_trueMAF20_10NA.HEM.Thresh.fix.csv")
LA2176.thresh <- read.csv("SlBc_LA2176_trueMAF20_10NA.HEM.Thresh.csv")
LA0480.thresh <- read.csv("SlBc_LA0480_c_trueMAF20_10NA.HEM.Thresh.csv")
HEMthresh <- HEMthresh[,-c(1)]
LA2176.thresh <- LA2176.thresh[,-c(1)]
LA0480.thresh <- LA0480.thresh[,-c(1)]
all.thresh <- merge(HEMthresh, LA2176.thresh, on="X.1")
all.thresh <- merge(all.thresh, LA0480.thresh, on="X.1")
all.thresh <- all.thresh[,-c(6, 13)]
all.thresh <- rename(all.thresh, c("X.Phenos...8.."="LA2176", "Phenos...13."= "LA0480"))

write.csv(all.thresh, "SlBc_12plants_trueMAF20_10NA.HEM.Thresh.final.csv")
