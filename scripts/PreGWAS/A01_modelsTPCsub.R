#Nicole E Soltis
#mixed models, linear models, single-isolate models for TPC resub.

#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
ModDat <- read.csv("data/preGWAS/SlBc_ModelData.csv")

#remove isolates missing from one experiment
SlSumm <- as.data.frame(with(ModDat, table(Igeno,ExpBlock)))
#missing Exps: Gallo3 (96a), 94.1 (96b)
OgDat <- ModDat
ModDat <- subset(ModDat, Igeno != "Gallo3")
ModDat <- subset(ModDat, Igeno != "94.1")
ModDat <- ModDat[,-c(1)]

library(lme4); library(car); library(lmerTest)

#models I used for the paper
#linear model, fixed effects only
#this is the one used for the paper table 1: no random effects
#and type 2 ANOVA
Sys.time()
fullmod2 <- lm(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + ExpBlock/PExpRep.x + ExpBlock:Igeno + ExpBlock:Species/PlGenoNm, data = ModDat)
Sys.time()
#sink(file='output/newANOVA/fullmod_012417_BuildingModels.txt')
#print("Model: Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + ExpBlock/PExpRep.x + ExpBlock:Igeno + ExpBlock:Species/PlGenoNm")
#Sys.time()
#summary(fullmod) # the code generating output
#Anova(fullmod, type=2)
#anova(fullmod)
#Sys.time()
#sink()

#--------------------------------------------------------------------------
#testing mixed model here
#missing terms: 1|ExpBlock/PExpRep.x 1|ExpBlock:Species/PlGenoNm
mystarttime <- Sys.time()
rownames(ModDat) = make.names(rownames(ModDat), unique=TRUE)
mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock:Species), data = ModDat)
sink(file='results/output/modtest_100118.txt')
print(mystarttime)
print(Sys.time())
print("mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock/Species/PlGenoNm), data = ModDat)")
#rand{lmerTest} is deprecated. now use ranova{lmerTest}
ranova(mymmod)
Anova(mymmod, type=2)
anova(mymmod)
print(Sys.time())
sink()

#an additional mixed model- try adding 1|ExpBlock/PExpRep.x
#this one fails at ranova step
#missing terms: 1|ExpBlock:Species 1|ExpBlock:Species/PlGenoNm
mystarttime <- Sys.time()
rownames(ModDat) = make.names(rownames(ModDat), unique=TRUE)
mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock/PExpRep.x), data = ModDat)
#trying this to deal with ranova error of duplicate row.names
rownames(mymmod) = make.names(c(1:1000), unique=TRUE)
sink(file='results/output/modtest_100218.txt')
print(mystarttime)
print(Sys.time())
print("mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock/Species/PlGenoNm), data = ModDat)")
#rand{lmerTest} is deprecated. now use ranova{lmerTest}
ranova(mymmod)
Anova(mymmod, type=2)
anova(mymmod)
print(Sys.time())
sink()

#this one errors out, won't run model, let alone run anova/ ranova
#missing terms: 1|ExpBlock/PExpRep.x 
mystarttime <- Sys.time()
rownames(ModDat) = make.names(rownames(ModDat), unique=TRUE)
mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock:Species) + (1|ExpBlock:Species/PlGenoNm), data = ModDat)

#testing one more mixed model to match lsmeans
mystarttime <- Sys.time()
rownames(ModDat) = make.names(rownames(ModDat), unique=TRUE)
mymmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|IndPlant/Leaf/AorB), data = ModDat)
sink(file='results/output/modtest_lsmterms_101718.txt')
print(mystarttime)
print(Sys.time())
print("mymmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|IndPlant/Leaf/AorB), data = ModDat)")
#rand{lmerTest} is deprecated. now use ranova{lmerTest}
ranova(mymmod)
Anova(mymmod, type=2)
anova(mymmod) #won't do
print(Sys.time())
sink()

anova(mymmod, type=2)
Anova(mymmod, type=3)
#-------------------------------------------------------------------------------
#new model: try dropping 2 "domestication sensitive" isolates, rerun fixed fx model
#from the text, domestication sensitive isolates are: Fd2, Rose
ModDat.rmD <- subset(ModDat, Igeno != c("Fd2","Rose"))


library(lme4); library(car); library(lmerTest)

#models I used for the paper
#linear model, fixed effects only
#this is the one used for the paper table 1: no random effects
#and type 2 ANOVA
Sys.time()
fullmod2.rmD <- lm(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + ExpBlock/PExpRep.x + ExpBlock:Igeno + ExpBlock:Species/PlGenoNm, data = ModDat.rmD)
Sys.time()
sink(file='results/output/fixmod_100218_rmDomestIsos.txt')
print("Model: fullmod2.rmD <- lm(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + ExpBlock/PExpRep.x + ExpBlock:Igeno + ExpBlock:Species/PlGenoNm, data = ModDat.rmD)")
Sys.time()
summary(fullmod2.rmD) # the code generating output
Anova(fullmod2.rmD, type=2)
anova(fullmod2.rmD)
Sys.time()
sink()

#-------------------------------------------------------------------------------
#experiment with additional models here. not using these for the paper.
#below here: playing with reconfiguring the models
#mixed model
mystarttime <- Sys.time()
mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno*Species + Igeno*Species/PlGenoNm + (1|PExpRep.x) , data = ModDat)

mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno*Species + Igeno*Species/PlGenoNm + (1|PExpRep.x) + (1|PExpRep.x*Igeno) + (1|PExpRep.x*Species) + (1|PExpRep.x*Species/PlGenoNm) + (1|Species/PlGenoNm/IndPlant), data = ModDat)

randmodtest <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/PExpRep.x) + (1|ExpBlock:Igeno) + (1|ExpBlock:Species/PlGenoNm), data = ModDat)

randmodtest <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock:Species), data = ModDat)
#fails with 1|exp + 1|exp/rep + 1|exp:igeno + 1|exp:species/pgeno 
#fails with 1|exp + 1|exp:igeno + 1|exp:species/pgeno 
#may work with 1|exp + 1|exp:igeno + 1|exp:species

#here's the fixed effect model again, for reference
#Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + ExpBlock/PExpRep.x + ExpBlock:Igeno + ExpBlock:Species/PlGenoNm


#fixed effect model (linear model)
mystarttime <- Sys.time()
mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno*Species + Igeno*Species/PlGenoNm + (PExpRep.x) , data = ModDat)

mymmod <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + Igeno*Species + Igeno*Species/PlGenoNm + (PExpRep.x) + (PExpRep.x*Igeno) + (PExpRep.x*Species) + (PExpRep.x*Species/PlGenoNm) + (Species/PlGenoNm/IndPlant), data = ModDat)

sink(file='results/output/modtest_092818.txt')
print(mystarttime)
print(Sys.time())
print("")
#rand{lmerTest} is deprecated. now use ranova{lmerTest}
ranova(mymmod)
Anova(mymmod, type=2)
anova(mymmod)
print(Sys.time())
sink()