#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/PhD/Research/Eudicots/Solanum/Analysis/R")
MyDat <- read.csv("AllResultsSlBc96.csv")
IsoNm <- read.csv("IsoIDs.csv")
PlantNm <- read.csv("PlantIDs.csv")
AddUnits <- read.csv("AddUnits.csv")
?anova

names(MyDat)
#subset to include only lesion size
LsDat <- MyDat[,c(1:11,152)]
#Add pixels per cm conversion
names(LsDat)
names(AddUnits)
LsDat$Sort <- paste(LsDat$PExpRep, LsDat$PImage, sep='') 
AddUnits$Sort <- paste(AddUnits$PExpRep, AddUnits$Image, sep='')
#unique(unlist(LsDat$Sort))
#unique(unlist(AddUnits$Sort))
LsDat2 <- merge(LsDat, AddUnits, by="Sort")
LsDat <- LsDat2[,c(2:13,19)]
#Add column of lesion size in cm squared
LsDat <- transform(LsDat, Scale.LS=(Lesion.Size/(pixelsPcm^2)))
library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr)
library(dplyr)

#add actual isolate names
unique(unlist(LsDat$Pexp))
LsDat96a <- filter(LsDat, Pexp == "96a") 
LsDat96b <- filter(LsDat, Pexp == "96b") 
IsoNm96a <- IsoNm
colnames(IsoNm96a)[1] <- "Piso"
IsoNm96b <- IsoNm
colnames(IsoNm96b)[3] <- "Piso"
#why would I gain data at this step?
LsDat96a2 <- merge(LsDat96a, IsoNm96a, by="Piso")
LsDat96b2 <- merge(LsDat96b, IsoNm96b, by="Piso")
LsDat96a2 <-LsDat96a2[,c(1:14,18,20)]
LsDat96b2 <-LsDat96b2[,c(1:14,18,20)]
SrtDat <- rbind(LsDat96a2, LsDat96b2)

#add a column for plant*iso interaction
SrtDat$PbyI <- paste(SrtDat$PPlant, SrtDat$Isolate, sep='') 

#remove any duplicate entries
SrtDat <- unique(SrtDat)

#prior to building MyPlot, eliminate all P*I with <2 replicates measured
#because can't estimate variance
#SrtDatpl <- SrtDat[SrtDat$n >= 3,]

#----------------------------------------------------------------------
#linear model
#nesting: B within A as A/B or A + A:B
#fixed effects: PInLflt
#random effects: PPlant, Isolate, PInPlant, PInLeaf, Pexp

#add a domestication term
names(SrtDat)
unique(SrtDat$PPlant)
#FL NC TX MA KS OR are Domest
#IL UT PA SD NY CA are Wild
SrtDat$Domest <- ifelse(SrtDat$PPlant %in% c("FL","NC","TX","MA","KS","OR"),"Dm", "Wl")

#remove Control Isolates
unique(SrtDat$Isolate)
SrtDat <- SrtDat[SrtDat$Isolate!="Control",]

#Anova for F-stats
ModDat <- SrtDat
names(ModDat)
ModDat <- dplyr::select(ModDat, ExpBlock = Pexp, Igeno = Isolate, Pgeno = PPlant, AorB = PInLflt, Leaf = PInLeaf, Plant = PInPlant, AgFlat = PImage, Species = Domest, IsoColor= OrderNJ, matches("."))
#remove any rows with NA for plant or isolate
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
ModDat <- completeFun(ModDat, c("Pgeno", "Igeno"))

#add a column of correct plant geno names
head(ModDat)
PlGenNum <- dplyr::select(PlantNm, Pgeno = PlantID, matches("."))
PlGenNum$PlGenoNm <- paste("LA", PlGenNum$PlantGeno, sep='') 
PlGenNum <- PlGenNum[c(1:12),c("Pgeno","PlGenoNm")]
ModDat <- merge(ModDat, PlGenNum, by="Pgeno")

#-------------------------------------------------
#plot my data: scatterplot!
names(ModDat)
attach(ModDat)
#add PlantNum as an integer sorted by mean lesion size
FigDat3 <- ddply(ModDat, c("PlGenoNm", "Igeno", "Species", "IsoColor"), summarise,
                 mLS   = mean(Scale.LS))
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6)
FigDat3 <- merge(FigDat3, MDmeans, by="PlGenoNm")
FigDat3 <- dplyr::select(FigDat3, Plnum = mean, matches("."))
FigDat3$PlantNum <- as.numeric(FigDat3$PlNum)
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
ggplot(FigDat3, aes(x = PlNum, y = mLS))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(mmLS), group=factor(Igeno)), show_guide=F)+
  scale_x_discrete(breaks=c("1","2","3", "4", "5", "6", "1", "2", "3", "4", "5", "6"),
                   labels=c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", 
                    "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))+
  facet_grid(~SpLabs, scales="free_x", space="free_x")+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
#+geom_smooth(aes(group = 2), size = 2, method = "lm", se = T)

#-------------------------------------------------------
#barplot with average lesion size by genotype
#and SE bars
FigDat <- ModDat
names(FigDat)
FigDat2 <- ddply(FigDat, "PlGenoNm", summarise, mLS = mean(Scale.LS))
FigDat3 <- ddply(FigDat, "PlGenoNm", summarise, seLS = se(Scale.LS))

FigDat2 <- ddply(FigDat, c("PlGenoNm", "Species"), summarise,
               N    = length(Scale.LS),
               mean = mean(Scale.LS),
               sd   = sd(Scale.LS),
               se   = sd / sqrt(N))
FigDat2$SpLabs <- factor(FigDat2$Species, labels = c("Domesticated", "Wild"))
limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(FigDat2, aes(x = factor(PlGenoNm), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())+
  geom_errorbar(limits, width=0.25)+
  facet_grid(.~SpLabs, scales="free")

#---------------------------------------------------------
#scatterplot: mean of each isolate on each host // domestication level
ModDat$SbyI <- paste(ModDat$Species, ModDat$Igeno, sep='') 
FigDat4 <- ddply(ModDat, c("Pgeno", "Igeno", "Species", "SbyI","PbyI"), summarise,PlMean   = mean(Scale.LS))
FDmeans <- ddply(ModDat, c("SbyI"), summarise, SpMean=mean(Scale.LS))
FigDat4 <- merge(FigDat4, FDmeans, by="SbyI")
attach(FigDat4)
ggplot(FigDat4, aes(x = PlMean, y = SpMean, color=Species, group=Pgeno))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(Igeno)), show_guide=F)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x="Plant Genotype")+
  theme(axis.text.x = element_blank())

#or more simply...
# mean lesion size per domestication status
names(FigDat4)
FigDat5 <- FigDat4[!duplicated(FigDat4$PbyI), ]
#color lines by positive or negative
unique(unlist(FigDat5$Species))
FigDat5w <- filter(FigDat5, Species == "Wl") 
FigDat5d <- filter(FigDat5, Species == "Dm") 
FigDat5d <-FigDat5d[,c(1,3,4,7)]
FigDat5w <-FigDat5w[,c(1,3,4,7)]
unique(FigDat5w$Igeno)
FigDat5d <- FigDat5d[!duplicated(FigDat5d$Igeno), ]
FigDat5w <- FigDat5w[!duplicated(FigDat5w$Igeno), ]
FigDat6 <- merge(FigDat5w, FigDat5d, by="Igeno")
#get mean difference (Dm - Wl)
head(FigDat6)
FigDat6 <- transform(FigDat6, MeanDiff=(SpMean.y - SpMean.x))
FigDat6$Direction <- ifelse(FigDat6$MeanDiff > 0,1,0)
FigDat7 <- FigDat6[,c(1,9)]
FigDat5 <- merge(FigDat5, FigDat7, by="Igeno")

names(FigDat5)
ggplot(FigDat5, aes(x = Species, y = SpMean, group=Igeno))+
  geom_point()+
  geom_line(size=1, aes(color=factor(Direction)), show_guide=F)+
  theme_bw()+
  labs(y=expression(Mean ~ Lesion ~ Area ~ per ~ Isolate ~ (cm^{2})), x=element_blank())+
  scale_x_discrete(labels=c("Domesticated","Wild"))

#---------------------------------------------------------------------------------------------

library(ggplot2)
ggplot (data = ModDat, 
  aes(x=GenFx, y=bigH))+
  geom_violin(adjust = 0.5, scale = "width", fill="#E6F598")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  stat_summary(fun.y="median", geom="point")+
  labs(y=expression(Heritability~(H^2)), x=NULL)+
  ylim(c(0,NA))+ #because anything less than 0 is uninteresting
  #and just means that variance for some genos > total variance
  scale_fill_manual(
    values = c("#E6F598","#E6F598", "#E6F598"))+
#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)
geom_jitter(shape=16, position=position_jitter(0.2))


library(vioplot)
#use complete cases only
MyPlotcc <- MyPlot[complete.cases(MyPlot),]
#set H<0 to 0
MyPlotcc$bigH[MyPlotcc$bigH<0] <- 0
x1 <- MyPlotcc$bigH[MyPlotcc$GenFx=="Isolate"]
x2 <- MyPlotcc$bigH[MyPlotcc$GenFx=="Plant"]
x3 <- MyPlotcc$bigH[MyPlotcc$GenFx=="PbyI"]
vioplot(x1, x2, x3, names=c("Isolate", "Plant", "PbyI"), 
        col="green")