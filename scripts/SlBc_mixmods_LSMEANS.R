#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
#read in file from SlBc_mixmods.csv
ModDat <- read.csv("SlBcDATAFRAME.csv")

#----------------------------------------------------------------------

library(lme4); library(car); library(lmerTest)

#rename some things
ModDat <- ModDat[,c("ExpBlock", "Igeno", "AorB", "Leaf", "Plant", "AgFlat", "Species", "IsoColor", "PExpRep.x", "Scale.LS", "PlGenoNm")]
names(ModDat)
ModDat <- dplyr::select(ModDat, Exp = ExpBlock, Block = PExpRep.x, Pgeno = PlGenoNm, matches("."))

#Anova(Mod, type=2)
#anova(Mod)
#Variance output of summary(Mod) gives you SS for the random factors
#rand(Mod) gives Chi-sq  and P values for random factors in packages lmerTest
#summary(fullmod) # the code generating output
#Sys.time()
#sink()

#linear model
#nesting: B within A as A/B or A + A:B
#fixed effects: PInLflt
#random effects: PPlant, Isolate, PInPlant, PInLeaf, Pexp
#ExpBlock and AgFlat as random effects
#but maybe include a random term for "bench" = PExpRep.x??
#ExpBlock/Bench/AgFlat
#AgFlat is nested within ExpBlock
#and both are random
#Leaf is nested within Plant
#and both are random
#And we can consider Pgeno to be nested within Species
#Igeno, Species, Pgeno, AorB are fixed
#nesting terms are already included- don't need to add Species as a separate term BUT for random effects (ExpBlock alone) do need a separate term
#PlGenoNm is a term nested within Species (but not CODED as if nested within Species = it's not an implicitly nested factor)
#AgFlat IS implicitly nested within ExpBlock -- let's fix this
#non-numeric factors to use: Igeno, PlGenoNm, Species, ExpBlock, AgFlat
#optional to fix: coding of AgFlat so that it is not implicitly nested

#fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#------------------------------------------------------------------------------
#lsmeans calculations

#removed just igeno * plant interactions
#ideally: Igeno, Species, PlGenoNm, AorB are fixed
#ExpBlock, AgFlat, and all their interactions are random
#Igeno with PlGenoNm interaction, PlGenoNm nested within Species, AgFlat nested within ExpBlock, interaction between ExpBlock and PlGenoNm as well as ExpBlock and Igeno



#run model per isolate WITHIN each plant genotype
#so include no species terms or plant genotype terms
attach(ModDat)
out <- split( ModDat , f = ModDat$Pgeno)
head(out[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
#sink(file='ModelsBYISO_030816rand.txt')
#adding AorB or Plant or Leaf: error
sink(file="LSMeans032116.txt")
for (i in c(1:12)) {
  print(unique(out[[i]]$Pgeno))
  Lesion.lm <- lmer(Scale.LS ~ Igeno + Plant/Leaf/AorB + (1|Exp), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "Igeno")
  print(Lesion.lsm)
}
sink()

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