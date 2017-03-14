#Nicole E Soltis
#Phenotype plotting Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
FigDat3 <- read.csv("data/BcSolGWAS_PhenotypePlotData_FigDat3.csv")
ModDat <- read.csv("data/BcSolGWAS_PhenotypePlotData_ModDat.csv")
FigDat4 <- read.csv("data/BcSolGWAS_PhenotypePlotData_FigDat4.csv")
#-------------------------------------------------
#-------------------------------------------------------
library(ggplot2)

#plot Domesticated vs. Wild only
FigDatD <- FigDat4[which(FigDat4$Species=="Dm"),]
FigDatW <- FigDat4[which(FigDat4$Species =="Wl"),]

wilcox.test(FigDatD$cvLS, FigDatW$cvLS, paired=T)

#----------------------------------------------------
#rank test of isolates
FigDat5_m <- FigDat4[,c("Igeno", "SpLabs", "IsoColor", "mLS", "Group")]
FigDat5_c <- FigDat4[,c("Igeno", "SpLabs", "cvLS")]
library(tidyr)
FigDat5_w <- spread(FigDat5_m, "SpLabs", "mLS")
library(plyr)
names(FigDat5_w)
FigDat5_w <- plyr::rename(FigDat5_w, c("Domesticated" = "m_Domest", "Wild" = "m_Wild"))
FigDat5_c <- spread(FigDat5_c, "SpLabs", "cvLS")
FigDat5_c <- plyr::rename(FigDat5_c, c("Domesticated" = "cv_Domest", "Wild" = "cv_Wild"))
FigDat5_w <- merge(FigDat5_w, FigDat5_c, by="Igeno")
FigDat5_w$dm <- FigDat5_w$m_Domest - FigDat5_w$m_Wild
FigDat5_w$dcv <- FigDat5_w$cv_Domest - FigDat5_w$cv_Wild

wilcox.test(FigDat5_w$m_Domest, FigDat5_w$m_Wild, paired=TRUE) 
wilcox.test(FigDat5_w$cv_Domest, FigDat5_w$cv_Wild, paired=TRUE) 
#no change for CV

mean(FigDat5_w$m_Domest)
mean(FigDat5_w$m_Wild)
(mean(FigDat5_w$m_Domest) - mean(FigDat5_w$m_Wild))/(mean(FigDat5_w$m_Wild))
