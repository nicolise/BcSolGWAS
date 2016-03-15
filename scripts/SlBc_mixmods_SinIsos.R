#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
MyDat <- read.csv("AllResultsSlBc96.csv")
IsoNm <- read.csv("IsoIDs.csv")
PlantNm <- read.csv("PlantIDs.csv")
AddUnits <- read.csv("AddUnits.csv")

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
library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr); library(dplyr)

#add actual isolate names
unique(unlist(LsDat$Pexp))
LsDat96a <- filter(LsDat, Pexp == "96a") 
LsDat96b <- filter(LsDat, Pexp == "96b") 
IsoNm96a <- IsoNm
colnames(IsoNm96a)[1] <- "Piso"
IsoNm96b <- IsoNm
colnames(IsoNm96b)[3] <- "Piso"
LsDat96a2 <- merge(LsDat96a, IsoNm96a, by="Piso")
LsDat96b2 <- merge(LsDat96b, IsoNm96b, by="Piso")
LsDat96a2 <-LsDat96a2[,c(1:14,18)] #18 should be Isolate
LsDat96b2 <-LsDat96b2[,c(1:14,18)]
SrtDat <- rbind(LsDat96a2, LsDat96b2)

#add a column for plant*iso interaction
SrtDat$PbyI <- paste(SrtDat$PPlant, SrtDat$Isolate, sep='') 

#remove any duplicate entries
SrtDat <- unique(SrtDat)

#add a domestication term
names(SrtDat)
unique(SrtDat$PPlant)
#FL NC TX MA KS OR are Domest
#IL UT PA SD NY CA are Wild
SrtDat$Domest <- ifelse(SrtDat$PPlant %in% c("FL","NC","TX","MA","KS","OR"),"Dm", "Wl")

#remove Control Isolates
unique(SrtDat$Isolate)
SrtDat <- SrtDat[SrtDat$Isolate!="Control",]

#rename columns as needed
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
names(ModDat)

#Plant is coded as numeric and nested within Pgeno
ModDat$IndPlant <- paste(ModDat$PlGenoNm, ModDat$Plant, sep='.') 
library(lme4); library(car); library(lmerTest)

#--------------------------------------------------------
#subset columns of interest 
ModDat <- ModDat[,c(2,3,4,5,6,7,8,10,13,14,15,16,17,18)]
head(ModDat)
#split dataset by isolate
out <- split( ModDat , f = ModDat$Igeno )
head(out[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file='ModelsBYISO_030816rand.txt')
#skip 13: blank, 59: 94.1, 68: blank, 77: Gallo3, 99: blank
for (i in c(1:12,14:58,60:67,69:76, 78:98, 100)) {
  print(unique(out[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ Species/PlGenoNm + (1|ExpBlock), data=out[[i]])
  result <- anova(Mod)
  random <- rand(Mod)
  print(result)
  print(random)
}
sink()

sink(file='ModelsBYISO_030916.txt')
for (i in c(1:12,14:58,60:67,69:76, 78:98, 100)) {
  print(unique(out[[i]]$Igeno))
  #adding AgFlat makes model worse
  #PExpRep does better than AgFlat as single term
  #don't always have enough reps for AorB
  #not meaningful (p=1) to add 1|ExpBlock/PlGenoNm
  Mod <- lmer(Scale.LS ~ Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock/PExpRep.x), data=out[[i]])
  result <- anova(Mod)
  random <- rand(Mod)
  print(result)
  print(random)
}
sink()

#-----------------------------------------------------------
#FDR cutoff
#p is a vector of p values
MyPvals <- read.csv("isoANOVApvals.csv")
names(MyPvals)
p <- MyPvals$pSpecies
MyPvals$pSpFDR <- p.adjust(p, method = "fdr", n = length(p))

pSpPl <- MyPvals$pSpPlant
MyPvals$pSpPlFDR <- p.adjust(pSpPl, method = "fdr", n = length(pSpPl))

bpSp <- MyPvals$BpSpecies
MyPvals$bpSpFDR <- p.adjust(bpSp, method = "fdr", n = length(bpSp))

bpSpPl <- MyPvals$BpSpPlant
MyPvals$bpSpPlFDR <- p.adjust(bpSpPl, method = "fdr", n = length(bpSpPl))

write.csv(MyPvals, "isoANOVAfdr.csv")
#final model: 
#lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock/Igeno) + (1|ExpBlock/PlGenoNm), data = ModDat)
#Igeno, PlGenoNm, Species, ExpBlock, AgFlat, IndPlant, AorB