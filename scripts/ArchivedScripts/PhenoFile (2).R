rm(list=ls())
setwd("~/Projects/BcSolGWAS/data/genome")
myData <- read.csv("LSMforbigRR_all.csv")

names(myData)
myData <- myData[,c("Pgeno", "Igeno", "Estimate")]
names(myData)

library(reshape2)
myDat<- dcast(myData, Igeno ~ Pgeno)
write.csv(myDat, "LSMforbigRR_est.csv")
