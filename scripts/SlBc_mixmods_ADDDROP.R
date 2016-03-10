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
#building up a minimal model to check overfitting problem

#final model: 
#lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock|Igeno) + (1|ExpBlock|PlGenoNm), data = ModDat)
#Igeno, PlGenoNm, Species, ExpBlock, AgFlat, IndPlant, AorB

attach(ModDat)
Lesion.lm.minA1 <- lmer(Scale.LS ~ Igeno + (1|ExpBlock))
LesIso.lsm.minA1 <- lsmeans(Lesion.lm.minA1, "Igeno")

Lesion.lm.minA2 <- lmer(Scale.LS ~ PlGenoNm + (1|ExpBlock))
LesIso.lsm.minA2 <- lsmeans(Lesion.lm.minA2, "PlGenoNm")

Lesion.lm.minA3 <- lmer(Scale.LS ~ Species + (1|ExpBlock))
LesIso.lsm.minA3 <- lsmeans(Lesion.lm.minA3, "Species")

Lesion.lm.minA4 <- lmer(Scale.LS ~ AorB + (1|ExpBlock))
LesIso.lsm.minA4 <- lsmeans(Lesion.lm.minA4, "AorB")

Lesion.lm.minA5 <- lmer(Scale.LS ~ PlGenoNm + Igeno + (1|ExpBlock))

Lesion.lm.minA6 <- lmer(Scale.LS ~ Species + Igeno + (1|ExpBlock))

attach(ModDat)
Lesion.lm.minXX<- lmer(Scale.LS ~ Igeno + Species + (1|IndPlant) + (1|ExpBlock/AgFlat)) 
#first issue detected: PlGenoNm + Species

Lesion.lm.minA14 <- lmer(Scale.LS ~ PlGenoNm + Igeno + AorB +  (1|ExpBlock))

Lesion.lm.minF14 <- lmer(Scale.LS ~ PlGenoNm + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat))

Lesion.lm.min11 <- lmer(Scale.LS ~ Species + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat))

Lesion.lm.min12 <- lmer(Scale.LS ~ PlGenoNm + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat) + (1|IndPlant))

Lesion.lm.min13 <- lmer(Scale.LS ~ Species + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat) + (1|IndPlant))

Lesion.lm.min14 <- lmer(Scale.LS ~ PlGenoNm + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min15 <- lmer(Scale.LS ~ PlGenoNm + Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min16 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|ExpBlock))

Lesion.lm.min20 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#current fullest model that works
Lesion.lm.min21 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min22 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + AorB + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#Warning message: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,: Model is nearly unidentifiable: large eigenvalue ratio - Rescale variables?
Lesion.lm.min17 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock:Igeno) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min18 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock:PlGenoNm) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min23 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + AorB + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#Warning messages: 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,: unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,: Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

Lesion.lm.min23 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|IndPlant) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min24 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|IndPlant) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min24 <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock))

#fixed-effect model matrix is rank deficient so dropping 1 columns / coefficients
Lesion.lm.min07 <- lmer(Scale.LS ~ PlGenoNm + Species + (1|ExpBlock))
Lesion.lm.min09 <- lmer(Scale.LS ~ PlGenoNm + Species + (1|AgFlat))
Lesion.lm.min09 <- lmer(Scale.LS ~ PlGenoNm + Species + (1|IndPlant))

#drop 2 
Lesion.lm.min17 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (ExpBlock:Igeno))

#drop 11
Lesion.lm.min19 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat) + IndPlant)

#drop 12
Lesion.lm.min010 <- lmer(Scale.LS ~ (Species/PlGenoNm) + Igeno + (1|ExpBlock))
Lesion.lm.min09 <- lmer(Scale.LS ~ Species/PlGenoNm + Species + (1|ExpBlock))

#drop 72
Lesion.lm.min22 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat) + Species/IndPlant)

#drop 792
Lesion.lm.min19 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (PlGenoNm/IndPlant) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#myx <- lsmeans(Lesion.lm.min1, trt.vs.ctrl~Igeno, adjust="none")
#trt.vs.ctrl~PlGenoNm|Igeno

#Error in forceSymmetric(2 * solve(IE2)) : 
#error in evaluating the argument 'x' in selecting a method for function #'forceSymmetric': Error in solve.default(IE2) : 
#  system is computationally singular: reciprocal condition number = 8.36832e-17


#non mixed-models
Lesion.lm.min4 <- lmer(Scale.LS ~ (1|IndPlant) + (1|ExpBlock))
Lesion.lm.min6 <- lmer(Scale.LS ~ (1|AgFlat) + (1|ExpBlock))


myx2 <- lsmeans(lsmMod0, pairwise ~ PlGenoNm + Igeno)
write.csv(myx2, "lsmeans_firsttry.csv")

  

Sys.time()
sink(file='output021716.txt')
Sys.time()
#summary(fullmod) # the code generating output
#Sys.time()
rand(fullmod)
Anova(fullmod, type=2)
anova(fullmod)
Sys.time()
sink()
