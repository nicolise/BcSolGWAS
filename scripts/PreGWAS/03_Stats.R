#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
ModDat <- read.csv("data/SlBcDATAFRAME.csv")
#get list of isolate groups
IsoGroups <- read.csv("data/IsolateGroups.csv")
IsoGroups <- IsoGroups[,1:4]
#-------------------------------------------------

names(ModDat)
attach(ModDat)
ModDat <- merge(ModDat, IsoGroups, by="Igeno")

#t-test of lesion size by species
t.test(ModDat$Scale.LS ~ ModDat$Species)

#quantiles of lesion size range by species
# function to get quantiles
qfun <- function(x, q = 5) {
  quantile <- cut(x, breaks = quantile(x, probs = 0:q/q), 
                  include.lowest = TRUE, labels = 1:q)
  quantile
}
ModDat[, qq := qfun(Scale.LS), by = Species]
tmp1 <- with(ModDat, split(Scale.LS, Species))
quantile(tmp1$Dm, 0.95)
quantile(tmp1$Dm, 0.05)
quantile(tmp1$Wl, 0.95)
quantile(tmp1$Wl, 0.05)

#t-test of tomato isolates vs. other isolates
#first, take summary data: mean of isolate per species
library(plyr)
ModDat2 <- ddply(ModDat, c("Igeno", "Species", "Tomato"), summarise,
                 meanLS = mean(Scale.LS))
ModDat2D <- ModDat2[ModDat2$Species=="Dm",]
ModDat2W <- ModDat2[ModDat2$Species=="Wl",]
t.test(ModDat2D$meanLS ~ ModDat2D$Tomato)
t.test(ModDat2W$meanLS ~ ModDat2W$Tomato)
#can't do multiple observations per isolate (pool D and W directly?)
#does no better when 1 mean per isolate
ModDat3 <- ddply(ModDat, c("Igeno", "Tomato"), summarise,
                 meanLS = mean(Scale.LS))
t.test(ModDat3$meanLS ~ ModDat3$Tomato)

#CV of lesion size between wild vs. domesticated tomato 
#wilcoxon signed-rank test
ModDat.cv <- ddply(ModDat, c("Igeno", "Species"), summarise, 
                   mLS = mean(Scale.LS),
                   sdLS = sd(Scale.LS))
ModDat.cv$cvLS <- ModDat.cv$sdLS/ModDat.cv$mLS
ModDat.cv$Species.Num <- as.numeric(ModDat.cv$Species)
t.test(ModDat.cv$cvLS, ModDat.cv$Species.Num)
wilcox.test(ModDat.cv$cvLS, ModDat.cv$Species.Num)

ModDat.cv <- ddply(ModDat, c("PlGenoNm", "Species"), summarise, 
                   mLS = mean(Scale.LS),
                   sdLS = sd(Scale.LS))

#----------------------------------------------------
#stuff from plotting I may not actually need here
#add PlantNum as an integer sorted by mean lesion size
library(plyr)
FigDat3 <- ddply(ModDat, c("PlGenoNm", "Igeno", "Species", "IsoColor"), summarise,
                 mLS   = mean(Scale.LS))
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlantNum <- c(1,2,3,4,5,6,7,8,9,10,11,12)
#MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6)
FigDat3 <- merge(FigDat3, MDmeans, by="PlGenoNm")
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


