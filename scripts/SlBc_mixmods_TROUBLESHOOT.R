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
#include or exclude Exp:Isolate and Exp:Plant
Sys.time()
Test1.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB , data = ModDat)
Sys.time()
Test2.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock:Igeno), data = ModDat)
Sys.time()
Test3.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock:PlGenoNm), data = ModDat)
Sys.time()
Test4.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock:PlGenoNm) + (1|ExpBlock:Igeno), data = ModDat)
Sys.time()
anova(Test1.lm, Test2.lm)
Sys.time()
anova(Test1.lm, Test3.lm)
Sys.time()
anova(Test2.lm, Test4.lm)
Sys.time()
anova(Test3.lm, Test4.lm)
Sys.time()
#include or exclude Exp/BLOCK/AgFlat
#include or exclude iso:Pgeno
Sys.time()
Test4.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock:PlGenoNm) + (1|ExpBlock:Igeno), data = ModDat)
Sys.time()
Test5.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock:PlGenoNm) + (1|ExpBlock:Igeno), data = ModDat)
Sys.time()
anova(Test4.lm, Test5.lm)
Sys.time()
#isolate as random or fixed effect

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

#------------------------------------------------------------------------------
#lsmeans calculations
#check for multicollinearity of variables
library("caret")
ModDat2 <- ModDat[,c("Scale.LS", "Igeno", "Species", "PlGenoNm", "ExpBlock","AgFlat","IndPlant", "AorB")]

#findLinearCombos can only deal with numeric variables, so convert each column
ModDat2$Igeno.f <- as.numeric(factor(ModDat2$Igeno))
ModDat2$Species.f <- as.numeric(factor(ModDat2$Species))
ModDat2$PlGenoNm.f <- as.numeric(factor(ModDat2$PlGenoNm))
ModDat2$ExpBlock.f <- as.numeric(factor(ModDat2$ExpBlock))
ModDat2$AgFlat.f <- as.numeric(factor(ModDat2$AgFlat))
ModDat2$IndPlant.f <- as.numeric(factor(ModDat2$IndPlant))
ModDat2$AorB.f <- as.numeric(factor(ModDat2$AorB))
ModDat3 <- ModDat2[,c("Scale.LS", "Igeno.f", "Species.f", "PlGenoNm.f", "ExpBlock.f","AgFlat.f","IndPlant.f", "AorB.f")]
findLinearCombos(ModDat3) #looks fine

#removed just igeno * plant interactions
#ideally: Igeno, Species, PlGenoNm, AorB are fixed
#ExpBlock, AgFlat, and all their interactions are random
#Igeno with PlGenoNm interaction, PlGenoNm nested within Species, AgFlat nested within ExpBlock, interaction between ExpBlock and PlGenoNm as well as ExpBlock and Igeno
Lesion.lm.01 <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock:PlGenoNm) + (1|ExpBlock:Igeno), data = ModDat)
#error: fixed-effect model matrix is rank deficient so dropping XXX columns / coefficients

#try calculating lsm separately within each plant genotype
Lesion.lm <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#library("lsmeans"); 
#also in lmerTest
#example(lsmeans)
lsm.options(pbkrtest.limit = 6285)
lsm.options(disable.pbkrtest=TRUE)
LesIso.lsm <- lsmeans(Lesion.lm, "Igeno")
s <- summary(LesIso.lsm)
class(s)
s[c("lsmean", "SE")]
LesIso.cn <- contrast(LesIso.lsm, "trt.vs.ctrlk") #pairwise comparisons vs. UKRazz
#not working for PlGenoNm

lsmeans(fm17, "Treatment")
pairs(.Last.value)

summary(LesIso.co)
lsmeans (mymod.lsm, 'Igeno', contr = "trt.vs.ctrlk")
#this is incorrect, replace 'Igeno' with what? 'specs' argument
pairs(mymod.lsm) #pairs of each isolate contrasted

mymod.lsm.Pl <- lsmeans(lsmMod0, "PlGenoNm")
contrast(mymod.lsm.Pl, "trt.vs.ctrlk") #pairwise comparisons
lsmeans (mymod.lsm, 'PlGenoNm', contr = "trt.vs.ctrlk")
warp.lsm <- lsmeans(warp.lm, ~ tension | wool)
#this is incorrect, replace 'Igeno' with what? 'specs' argument
pairs(mymod.lsm) #pairs of each isolate contrasted
myx <- lsmeans(lsmMod0, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
#Error in forceSymmetric(2 * solve(IE2)) : 
#error in evaluating the argument 'x' in selecting a method for function #'forceSymmetric': Error in solve.default(IE2) : 
#  system is computationally singular: reciprocal condition number = 8.36832e-17
myx2 <- lsmeans(lsmMod0, pairwise ~ Igeno) # + PlGenoNm?
write.csv(myx2, "lsmeans_firsttry.csv")

#more things to try
# warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
# contrast(warp.lsm, method = "poly", by = "wool")
# ( warp.lsm <- lsmeans (warp.lm,  ~ wool | tension, options = list(estName = "pred.breaks")) )
# (Oats.anal <- lsmeans(Oats.lme, list(poly ~ nitro, pairwise ~ Variety)))
# test(Oats.anal)
# lsmeans(warp.lm, eff ~ wool, adjust = "bonf") 
# lsmeans(warp.lm, eff ~ wool, options = list(adjust = "bonf"))
# summary(wool.lsm, adjust = "bonf")
# contrast(wool.lsm, "eff", adjust = "sidak")
# lsmeans(Oats.lme, ~ nitro, cov.reduce = FALSE, weights = "show.levels")

myx <- lsmeans(lsmMod0, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
#Error in forceSymmetric(2 * solve(IE2)) : 
#error in evaluating the argument 'x' in selecting a method for function #'forceSymmetric': Error in solve.default(IE2) : 
#  system is computationally singular: reciprocal condition number = 8.36832e-17
myx2 <- lsmeans(lsmMod0, pairwise ~ PlGenoNm + Igeno)
write.csv(myx2, "lsmeans_firsttry.csv")

library("lsmeans")
mod3.rg1 <- ref.grid(lsmMod3)
summary(mod3.rg1)
lsmeans(mod3.rg1, "Species")
lsmeans(mod3.rg1, "AorB")

myx <- lsmeans(lsmMod3, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
myx <- lsmeans(fullmod, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
#Error in lsmeans.character.ref.grid(object = <S4 object of class "ref.grid">,  : No variable named PlGenoNm in the reference grid
myx2 <- lsmeans(fullmod, pairwise ~ PlGenoNm + Igeno)
write.csv(myx, "lsmeans_firsttry.csv")
write.csv(myx2, "lsmeans_firsttry2.csv")
save.image(file = "0121.RData")

#writing the whole thing out to better understand nesting:
#Does it make sense to include Igeno:PlGenoNm separately? I think not.
Sys.time()
fullmodALT <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + Igeno:Species + Igeno:Species:PlGenoNm + (1|ExpBlock) + (1|ExpBlock:AgFlat) + (1|IndPlant) + AorB , data = ModDat)
save.image(file = "1119.RData") #model is saved as .RData
#use ?load to reload it
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output1117.txt')
Sys.time()
summary(fullmodALT) # the code generating output
Sys.time()
rand(fullmodALT)
Anova(fullmodALT, type=2)
anova(fullmodALT)
Sys.time()
sink()

library(sjPlot)
#
#minimized model to avoid overfitting
Sys.time()
minimod1 <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:AgFlat) + + AorB , data = ModDat)
save.image(file = "1122.RData") #model is saved as .RData
#use ?load to reload it
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output1122a.txt')
Sys.time()
summary(minimod1) # the code generating output
Sys.time()
rand(minimod1)
Sys.time()
Anova(minimod1, type=2)
Sys.time()
anova(minimod1)
Sys.time()
sink()

sjp.glmer(minimod1, type = "fe.cor")
sjp.glmer(minimod1, type = "re.qq")

#and remove all iso:plant terms
minimod2 <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + (1|ExpBlock) + (1|ExpBlock:AgFlat) + + AorB , data = ModDat)
save.image(file = "1122b.RData") #model is saved as .RData
#use ?load to reload it
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output1122b.txt')
Sys.time()
summary(minimod2) # the code generating output
Sys.time()
rand(minimod2)
Sys.time()
Anova(minimod2, type=2)
Sys.time()
anova(minimod2)
Sys.time()
sink()

sjp.glmer(minimod2, type = "fe.cor")
sjp.glmer(minimod2, type = "re.qq")

#fullmod including Isolate*Exp and Plant*Exp
#add a factor by isolate * experiment
names(ModDat)
ModDat$IbyX <- paste(ModDat$Igeno, ModDat$ExpBlock, sep='') 
sort(unique(ModDat$IbyX))
#remove isolates with no match across experiments:
#94.1, Gallo3
ModDat2 <- ModDat[ModDat$Igeno!="94.1",]
ModDat2 <- ModDat2[ModDat2$Igeno!="Gallo3",]

library(lme4); library(lmerTest)
library(car)
fullmod5 <- lmer(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant) + AorB + (1|ExpBlock) + Igeno:ExpBlock + Species:ExpBlock + ExpBlock:Species/Pgeno + ExpBlock:Igeno:Species/Pgeno, data = ModDat2)
#copy output to a text file
sink(file='output1a.txt')
summary(fullmod5) # the code generating output
sink()

sink(file='output1b.txt')
anova(fullmod5) # the code generating output
sink()
sink(file='output1c.txt')
Anova(fullmod5, type=2)
sink()
sink(file='output1d.txt')
rand(fullmod5)
sink()


fullmod2 <- lmer(Scale.LS ~ Species/Pgeno + Igeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = ModDat)
fullmod
summary(fullmod)
rand(fullmod)


anova(fullmod)
anova(fullmod2)
library(car)
Anova(fullmod, type=2)
qqnorm(residuals(fullmod))
qqline(residuals(fullmod))
plot(fullmod)
print(fullmod)
summary(fullmod)

#try recoding random effects to infer nesting and random
ModDat2 <- ModDat
ModDat2$LfinPlant <- paste(ModDat$Plant, ModDat$Leaf, sep='.')
ModDat2$FlatinBlock <- paste(ModDat$ExpBlock, ModDat$AgFlat, sep='')
names(ModDat2)
#ModDatTmp <- ModDat2[,c("Scale.LS","ExpBlock","Igeno","Pgeno","AorB","Leaf","Plant","AgFlat","Species","LfinPlant","FlatinBlock")]
#write.csv(ModDatTmp,"SoltisData.csv")
unique(ModDat2$FlatinBlock)
unique(ModDat2$LfinPlant)
#and now recode random effects as integers
ModDat2$LfinPlant <- as.numeric(ModDat2$LfinPlant)
ModDat2$FlatinBlock <- as.numeric(factor(ModDat2$FlatinBlock))
ModDat2$Plant <- as.numeric(ModDat2$Plant))
ModDat2$ExpBlock <- as.numeric(factor(ModDat2$ExpBlock))

#and model again
fullmod3 <- lm(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + ExpBlock + FlatinBlock + Plant + LfinPlant + AorB, data = ModDat2)
fullmod
anova(fullmod3)
redmod3 <- lm(Scale.LS ~ Igeno + Species/Pgeno + Species:Igeno + ExpBlock + FlatinBlock + Plant + LfinPlant + AorB, data = ModDat2)
summary(redmod3)
redmod4 <- lm(Scale.LS ~ Igeno + Species/Pgeno + ExpBlock + FlatinBlock + Plant + LfinPlant + AorB, data = ModDat2)
summary(redmod4)
Anova(redmod4, type=2)
redmod5 <- lm(Scale.LS ~ Igeno + Species/Pgeno + ExpBlock + FlatinBlock + Plant + AorB, data = ModDat2)
summary(redmod5)
Anova(redmod5, type=2)

library(car)
Anova(fullmod3, type=2)

aovmod1 <- aov(lm(Lesion.Size ~ Pexp + Pexp/PImage + PPlant*Isolate, data=SrtDat3))
summary(aovmod1)

#an alternate model
mod1 <- lm(Lesion.Size ~ PPlant*Isolate + PPlant/PInPlant/PInLeaf + PInLflt + Pexp , data=SrtDat3)
#lme for lsmeans
library(lme4)
mod1 <- lmer(Lesion.Size ~ PPlant*Isolate + (1|PImage) + (1|Pexp), data=SrtDat3)
summary(mod1) 
#lsmeans for GWAS
#doesn't work yet
library(lsmeans)
#myLSmeans <- lsmeans(mod1, Lesion.Size ~ PPlant|Isolate, adjust="none")
