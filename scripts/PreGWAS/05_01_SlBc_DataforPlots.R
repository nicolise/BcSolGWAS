#Nicole E Soltis
#Phenotype plotting Solanum spp. x Botrytis cinerea lesions:: prepare dataframe
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
ModDat <- read.csv("data/SlBcDATAFRAME.csv")
#-------------------------------------------------


#plot my data: scatterplot!
names(ModDat)
attach(ModDat)
names(ModDat)

ModDat$SpLabs <- factor(ModDat$Species, labels = c("Domesticated", "Wild"))
ModDat$SpLabs <- factor(ModDat$SpLabs, levels =c("Wild", "Domesticated"))
#add PlantNum as an integer sorted by mean lesion size
library(plyr)
FigDat3 <- ddply(ModDat, c("PlGenoNm", "Igeno", "Species", "IsoColor"), summarise,
                 mLS   = mean(Scale.LS))
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlantNum <- c(1,2,3,4,5,6,7,8,9,10,11,12)
#MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6)
FigDat3 <- merge(FigDat3, MDmeans, by="PlGenoNm")

names(MDmeans)
MDmeans <- MDmeans[,c("PlGenoNm","PlantNum")]
ModDat <- merge(ModDat,MDmeans, by="PlGenoNm" )
ModDat$PlantNum <- as.numeric(ModDat$PlantNum)

FigDat3 <- dplyr::select(FigDat3, Plnum = mean, matches("."))
FigDat3$PlantNum <- as.numeric(FigDat3$PlantNum)
library(ggplot2)
attach(FigDat3)
names(FigDat3)
FigDat3$SpLabs <- factor(FigDat3$Species.x, labels = c("Domesticated", "Wild"))
#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$Igeno)
attach(FigDat3)
FigDat3 <- FigDat3[order(mmLS),]

#make x axis labels that actually work
library(plyr)
FigDat3$Plant.Label <- mapvalues(FigDat3$PlantNum, 
                                 c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                 c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))
FigDat3$Plant.Lab.Ord <- factor(FigDat3$Plant.Label, levels = c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))

#plant order is correct at this stage

#get list of isolate groups
IsoGroups <- read.csv("data/IsolateGroups.csv")
names(IsoGroups)
IsoGroups <- IsoGroups[,1:9]
FigDat3 <- merge(FigDat3, IsoGroups, by="Igeno")

FigDat3$SpLabs <- factor(FigDat3$SpLabs, labels = c("Domesticated", "Wild"))
FigDat3$SpLabs <- factor(FigDat3$SpLabs, levels = c("Wild", "Domesticated"))

write.csv(FigDat3, "data/BcSolGWAS_PhenotypePlotData_FigDat3.csv")
write.csv(ModDat, "data/BcSolGWAS_PhenotypePlotData_ModDat.csv")

library(plyr)
ModDat$Spec.Iso <- paste(ModDat$Species, ModDat$Igeno, sep="")
FigDat4 <- ddply(ModDat, c("Igeno", "Species", "Spec.Iso", "IsoColor"), summarise,
                 mLS   = mean(Scale.LS),
                 sdLS = sd(Scale.LS))
FigDat4$cvLS <- FigDat4$sdLS / FigDat4$mLS
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
FigDat4$SpLabs <- factor(FigDat4$Species, labels = c("Domesticated", "Wild"))
FigDat4$SpLabs <- factor(FigDat4$SpLabs, levels = c("Wild", "Domesticated"))

#get list of isolate groups
IsoGroups <- read.csv("data/IsolateGroups.csv")
IsoGroups <- IsoGroups[,1:2]
FigDat4 <- merge(FigDat4, IsoGroups, by="Igeno")
write.csv(FigDat4, "data/BcSolGWAS_PhenotypePlotData_FigDat4.csv")
