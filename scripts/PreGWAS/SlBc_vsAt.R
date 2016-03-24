#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
MyDat <- read.csv("AllResultsSlBc96.csv")
IsoNm <- read.csv("IsoIDs.csv")
PlantNm <- read.csv("PlantIDs.csv")
AddUnits <- read.csv("AddUnits.csv")
AtDat <- read.csv("96BcRAWlesiondata.csv")
CiDat <- read.csv("CiModelData.csv")
GmDat <- read.csv("GmModelData.csv")

names(AtDat)
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

#subset to include only lesion size
AtDat <- AtDat[,c(1:10,151)]

#add a column of lesion size in cm squared for AtDat
AtDat <- transform(AtDat, Scale.LS=(Lesion.Size/(Scale..Pix.cm.^2)))

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
ModDat <- dplyr::select(ModDat, ExpBlock = Pexp, Igeno = Isolate, Pgeno = PPlant, AorB = PInLflt, Leaf = PInLeaf, Plant = PInPlant, AgFlat = PImage, Species = Domest, matches("."))
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
#-----------------------------------------------------------------
#Now to combine data sets
names(ModDat)
List1 <- unique(ModDat$Igeno)
List2 <- unique(AtDat$Isolate)
List3 <- unique(CiDat$IsolateID)
List4 <- unique(GmDat$isolatename)
AtDat <- AtDat[AtDat$Isolate!="Ctrl",]
df1 <- data.frame(matrix(unlist(List1), nrow=97, byrow=T),stringsAsFactors=FALSE)
df2 <- data.frame(matrix(unlist(List2), nrow=98, byrow=T),stringsAsFactors=FALSE)
df3 <- data.frame(matrix(unlist(List3), nrow=98, byrow=T),stringsAsFactors=FALSE)
df4 <- data.frame(matrix(unlist(List4), nrow=96, byrow=T),stringsAsFactors=FALSE)

#remove unmatched from ModDat
#AtDatOG <- AtDat
AtDat<- AtDatOG
#ModDat <- ModDat[ModDat$Igeno!=c("01.02.05","94.1"),]
ModDat <- ModDat[ModDat$Igeno!=c("Gallo3"),]
#AtDat <- AtDat[AtDat$Isolate!="01.02.13", "01.02.19", ,]
ModDat$Igeno <- revalue(ModDat$Igeno, c("MEAP6G" = "MEAPGG"))
AtDat$Isolate <- revalue(AtDat$Isolate, c("Apple 404" = "Apple404", "Triple 3 (T3)" = "Triple3", "Davis Navel" = "DavisNavel", "Gallo 2" = "Gallo2", "Kern A2" = "KernA2", "Kern B2" = "KernB2", "Kern B1" = "KernB1", "Apple 517" = "Apple517", "UK razz" = "UKRazz", "Triple 7 (T7)" = "Triple7", "Philo Menlo" = "PhiloMenlo", "Pepper sub" = "PepperSub", "Noble Rot" = "NobleRot", "Mex-03" = "Mex03", "MEA PGG" = "MEAPGG", "Katie tomato" = "KatieTomato", "Gallo 1" = "Gallo1", "Gallo 2" = "Gallo2", "Fresa S.D." = "FresaSD", "Fresa 525" = "Fresa525", "Esparato Fresa" = "EsparatoFresa", "BPA 1" = "BPA1"))
CiDat$IsolateID <- revalue(CiDat$IsolateID, c("FresaS.D." = "FresaSD", "Katietomato" = "KatieTomato", "MEAP6G" = "MEAPGG", "Mex-03" = "Mex03", "Peppersub" = "PepperSub", "Tripple3(T3)" = "Triple3", "Tripple7(T7)" = "Triple7"))
GmDat$isolatename <- revalue(GmDat$isolatename, c("Apple 404" = "Apple404", "Triple 3 (T3)" = "Triple3", "Davis Navel" = "DavisNavel", "Gallo 2" = "Gallo2", "Kern A2" = "KernA2", "Kern B2" = "KernB2", "Kern B1" = "KernB1", "Apple 517" = "Apple517", "UK razz" = "UKRazz", "Triple 7 (T7)" = "Triple7", "Philo Menlo" = "PhiloMenlo", "Pepper sub" = "PepperSub", "Noble Rot" = "NobleRot", "MEA P6G" = "MEAPGG", "Katie Tomato" = "KatieTomato", "Gallo 1" = "Gallo1", "Gallo 2" = "Gallo2", "Fresa S.D." = "FresaSD", "Fresa 525" = "Fresa525", "Esparato Fresa" = "EsparatoFresa", "BPA 1 " = "BPA1"))
#keep only: Experiment, Plant, Isolate, Scale.LS
#and ExpBlock, Igeno, Pgeno, Species, Scale.LS
names(AtDat)
names(ModDat)
names(GmDat)
AtDat2 <- AtDat[,c("Isolate", "Plant", "Scale.LS")]
SlDat2 <- ModDat[,c("Igeno", "Pgeno", "Scale.LS", "Species")]
GmDat2 <- GmDat[,c("isolatename", "plantnumber", "Scale.LS")]
CiDat2 <- CiDat[,c("IsolateID", "PlantGeno", "Scale.LS")]
SlDat2 <- dplyr::select(SlDat2, Isolate = Igeno, Plant = Pgeno, matches("."))
GmDat2 <- dplyr::select(GmDat2, Isolate = isolatename, Plant = plantnumber, matches("."))
CiDat2 <- dplyr::select(CiDat2, Isolate = IsolateID, Plant = PlantGeno, matches("."))

AtDat2$Species <- "At"
CiDat2$Species <- "Ci"
GmDat2$Species <- "Gm"

FullDat <- rbind(SlDat2, AtDat2, CiDat2, GmDat2)

#means by isolate*plant
#violin plot for Wild, Domest, At
#add "Plant" to list of Les.means factors if want to look within plant genos
names(FullDat)
Les.means <- ddply(FullDat, c("Isolate", "Species"), summarise, mean=mean(Scale.LS))
Les.means$NumSp <- ifelse(Les.means$Species == "At", 5, 
                          ifelse (Les.means$Species == "Ci", 4,
                                  ifelse (Les.means$Species == "Gm", 3,
                                     ifelse (Les.means$Species == "Dm", 1, 2))))

#add Les.means scaled to mean of each
tapply(Les.means$mean, Les.means$Species, mean)
Les.means$GpMean <- ifelse(Les.means$Species == "At", 0.2649255, 
                           ifelse(Les.means$Species == "Gm", 0.3300405,
                                  ifelse(Les.means$Species == "Ci", 1.1561966,
                           ifelse (Les.means$Species == "Dm", 0.7309409, 0.6008727))))
Les.means <- transform(Les.means, RelLes=(mean/GpMean))

#add a factor for species x exp
#Les.means2$ExpbySp <- paste(Les.means2$NumSp, Les.means2$Experiment, sep='')


#-------------------------------------------------------------------------
#scatterplot
myitalic <- element_text(face="italic")
names(Les.means)
#add a column of mRL (mean of Relative Lesion Size) per isolate
#sort dataframe by mRL
#then color by the new factor mRL
Les.means$mRL<- ave(Les.means$RelLes, Les.means$Isolate)
attach(Les.means)
Les.means <- Les.means[order(mRL),]
attach(Les.means)
ggplot(Les.means, aes(x = NumSp, y = RelLes))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(mRL), group=factor(Isolate)), show_guide=F)+
  scale_x_discrete(breaks=c("1","2","3", "4","5"),
                   labels=c(expression(paste(italic("S. lycopersicum"))), expression(paste(italic("S. pimpinellifolium"))), expression(paste(italic("G. max"))), expression(paste(italic("C. intybus"))), expression(paste(italic("A. thaliana")))))+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 15, hjust = 1, vjust=1))+
  labs(y=expression(Relative ~ Lesion ~ Area), x=element_blank())
#-----------------------------------------------
#scatterplot by experiment
attach(Les.means2)
Les.means3 <- Les.means2[Les.means2$Species!="At",]
names(Les.means3)
Les.means3$ExpbyIso <- paste(Les.means3$Experiment, Les.means3$Isolate, sep='')
file <- ddply(Les.means3, c("Experiment", "Isolate"), summarise, Newmean = mean(mean))
Les.means4 <- Les.means3[Les.means3$Isolate!="EsparatoFresa",]

ggplot(Les.means4, aes(x = ExpbySp, y = mean))+
  geom_point()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(Isolate), group=factor(Isolate)), show_guide=F)+
  scale_x_discrete(breaks=c("196a","196b","296a","296b"),
                   labels=c(expression(paste(italic("S. lycopersicum"))~ E1), expression(paste(italic("S. lycopersicum"))~ E2), expression(paste(italic("S. pimpinellifolium"))~E1), expression(paste(italic("S. pimpinellifolium"))~E2)))+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 15, hjust = 1, vjust=1))+
  labs(y=expression(Mean ~ Lesion ~ Area~(cm^2)), x=element_blank())
#--------------------------------------------------------------------------------------
#violin plot

library(ggplot2)
attach(Les.means)
ggplot (data = Les.means, 
  aes(x=factor(NumSp), y=mean))+
  geom_violin(adjust = 0.7, scale = "width", fill="mediumspringgreen")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 15, hjust = 1, vjust=1),
        text = element_text(size=24))+
  scale_x_discrete(breaks=c("1","2","3", "4","5"),
                   labels=c(expression(paste(italic("S. lycopersicum"))), expression(paste(italic("S. pimpinellifolium"))), expression(paste(italic("G. max"))), expression(paste(italic("C. intybus"))), expression(paste(italic("A. thaliana")))))+
  geom_boxplot(width=0.1)+
  labs(y=expression(Mean~Lesion~Area~(cm^2)), x=NULL)
