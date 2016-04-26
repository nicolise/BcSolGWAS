#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#this time for PCA
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
ModDat <- read.csv("SlBcDATAFRAME.csv")
Isolate <- read.csv("isoPLANTFXpcs.csv")
#-------------------------------------------------
names(ModDat)
attach(ModDat)
library(plyr)

#isolateA, 10 isolates
IsolateA <- Isolate[,c("Isolate","IsoSrtA", "myColorsA")]
IsolateA <- dplyr::select(Isolate, myColor = myColorsA, matches("."))
ModDatA <- dplyr::select(ModDat, Isolate = Igeno, matches("."))
ModDatA <- merge(ModDatA, IsolateA, by = "Isolate")
FigDatA <- ddply(ModDatA, c("PlGenoNm", "IsoSrtA", "Species", "myColor"), summarise, mLS   = mean(Scale.LS))
FigDatA <- dplyr::select(FigDatA, IsoSrtName = IsoSrtA, matches("."))
FigDat <- FigDatA
ModDat2 <- ModDatA
library(RColorBrewer)
mypalette <- (c("grey68",brewer.pal(10, "Set3")))

#isolateB, 4 isolates
IsolateB <- Isolate[,c("Isolate","IsoSrtB", "myColorsB")]
IsolateB <- dplyr::select(Isolate, myColor = myColorsB, matches("."))
ModDatB <- dplyr::select(ModDat, Isolate = Igeno, matches("."))
ModDatB <- merge(ModDatB, IsolateB, by = "Isolate")
FigDatB <- ddply(ModDatB, c("PlGenoNm", "IsoSrtB", "Species", "myColor"), summarise, mLS   = mean(Scale.LS))
FigDatB <- dplyr::select(FigDatB, IsoSrtName = IsoSrtB, matches("."))
FigDat <- FigDatB
ModDat2 <- ModDatB
library(RColorBrewer)
mypalette <- (c("grey68",brewer.pal(4, "Set2")))

#isolateC, 72 isolates
IsolateC <- Isolate[,c("Isolate","IsoSrtC", "myColorsC")]
IsolateC <- dplyr::select(Isolate, myColor = myColorsC, matches("."))
ModDatC <- dplyr::select(ModDat, Isolate = Igeno, matches("."))
ModDatC <- merge(ModDatC, IsolateC, by = "Isolate")
FigDatC <- ddply(ModDatC, c("PlGenoNm", "IsoSrtC", "Species", "myColor"), summarise, mLS   = mean(Scale.LS))
FigDatC <- dplyr::select(FigDatC, IsoSrtName = IsoSrtC, matches("."))
FigDat <- FigDatC
ModDat2 <- ModDatC
mypalette <- (c("grey68","turquoise1"))

#isolateD, 1 isolate
IsolateD <- Isolate[,c("Isolate","IsoSrtD", "myColorsD")]
IsolateD <- dplyr::select(Isolate, myColor = myColorsD, matches("."))
ModDatD <- dplyr::select(ModDat, Isolate = Igeno, matches("."))
ModDatD <- merge(ModDatD, IsolateD, by = "Isolate")
FigDatD <- ddply(ModDatD, c("PlGenoNm", "IsoSrtD", "Species", "myColor"), summarise, mLS   = mean(Scale.LS))
FigDatD <- dplyr::select(FigDatD, IsoSrtName = IsoSrtD, matches("."))
FigDat <- FigDatD
ModDat2 <- ModDatD
library(RColorBrewer)
mypalette <- (c("grey68","turquoise1"))

#isolateE, 8 isolates
IsolateE <- Isolate[,c("Isolate","IsoSrtE", "myColorsE")]
IsolateE <- dplyr::select(Isolate, myColor = myColorsE, matches("."))
ModDatE <- dplyr::select(ModDat, Isolate = Igeno, matches("."))
ModDatE <- merge(ModDatE, IsolateE, by = "Isolate")
FigDatE <- ddply(ModDatE, c("PlGenoNm", "IsoSrtE", "Species", "myColor"), summarise, mLS   = mean(Scale.LS))
FigDatE <- dplyr::select(FigDatE, IsoSrtName = IsoSrtE, matches("."))
FigDat <- FigDatE
ModDat2 <- ModDatE
library(RColorBrewer)
mypalette <- (c("grey68",brewer.pal(8, "Set2")))

#run this generically for each group
MDmeans<- ddply(ModDat2, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6)
FigDat3 <- merge(FigDat, MDmeans, by="PlGenoNm")
FigDat3 <- dplyr::select(FigDat3, Plnum = mean, matches("."))
FigDat3$PlantNum <- as.numeric(FigDat3$PlNum)
FigDat3$SpLabs <- factor(FigDat3$Species.x, labels = c("Domesticated", "Wild"))
FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$IsoSrtName)
FigDat3 <- FigDat3[order(FigDat3$mmLS),]

#the plot
library(ggplot2)
library(grid)

attach(FigDat3)
ggplot(FigDat3, aes(x = PlNum, y = mLS))+
  geom_point(color="grey68")+
  theme_bw()+
  geom_line(size=1, aes(group=factor(IsoSrtName), color=factor(myColor)), show.legend=F)+
  scale_x_discrete(limit=c("1","2","3", "4", "5", "6", "1", "2", "3", "4", "5", "6"))+
  facet_grid(~SpLabs, scales="free_x", space="free_x")+
  scale_colour_manual(values=mypalette)+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 0, hjust = 0.5), panel.margin=unit(2, "lines"))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())

#another way to pick lots of colors
#library(colorRamps)
#rgb.tables(72)