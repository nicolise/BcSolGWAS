#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
#read in file from SlBc_mixmods.csv
ModDat <- read.csv("data/preGWAS/SlBc_ModelData.csv")

#----------------------------------------------------------------------

library(lme4); library(car); library(lmerTest)

#rename some things
ModDat <- ModDat[,c("ExpBlock", "Igeno", "AorB", "Leaf", "Plant", "AgFlat", "Species", "IsoColor", "PExpRep.x", "Scale.LS", "PlGenoNm")]
names(ModDat)
ModDat <- dplyr::select(ModDat, Exp = ExpBlock, Block = PExpRep.x, Pgeno = PlGenoNm, matches("."))
#names(ModDat)[4]<-"Pgeno"
#names(ModDat)[8]<-"Exp"
#names(ModDat)[9]<-"Block"

#remove the isolates missing from one experiment
SlSumm <- as.data.frame(with(ModDat, table(Igeno,Exp)))
#missing Exps: Gallo3 (96a), 94.1 (96b)
OgDat <- ModDat
ModDat <- subset(ModDat, Igeno != "Gallo3")
ModDat <- subset(ModDat, Igeno != "94.1")

#lsmeans calculations

#removed just igeno * plant interactions
#ideally: Igeno, Species, PlGenoNm, AorB are fixed
#ExpBlock, AgFlat, and all their interactions are random
#Igeno with PlGenoNm interaction, PlGenoNm nested within Species, AgFlat nested within ExpBlock, interaction between ExpBlock and PlGenoNm as well as ExpBlock and Igeno

#run model per isolate WITHIN each plant genotype
#so include no species terms or plant genotype terms
attach(ModDat)
out <- split( ModDat , f = ModDat$Pgeno)

#try for one genotype alone, no loop
# mydat = out[[1]]
# print(unique(mydat$Pgeno))
# Lesion.lm <- lmer(Scale.LS ~ Igeno + (1|Exp) + (1|IndPlant/Leaf/AorB) + (1|Exp:Igeno), data=mydat)
# stats::anova(Lesion.lm)
# car::Anova(Lesion.lm, type=3)
# Lesion.lsm <- lsmeans(Lesion.lm, "Igeno")
# df <- as.data.frame(print(Lesion.lsm))

head(out[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
#sink(file='ModelsBYISO_030816rand.txt')
#fails to converge with igeno, plant/leaf/aorb, exp, exp/igeno
#cannot include exp/igeno
#also fails if drop aorb
#and if drop aorb and leaf

#version LSMeans061316.txt has fixed fx for AorB
d=NULL
library(data.table)
for (i in c(1:12)) {
  print(unique(out[[i]]$Pgeno))
  Lesion.lm <- lmer(Scale.LS ~ Igeno + (1|Exp) + (1|IndPlant/Leaf/AorB) + (1|Exp:Igeno), data=out[[i]])
  #Lesion.lm.2 <- lmer(Scale.LS ~ Igeno + (1|Exp) + (IndPlant/Leaf/AorB), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "Igeno")
  df <- as.data.frame(print(Lesion.lsm))
  setDT(df, keep.rownames = T)[]
  df$Plant <- unique(out[[i]]$Pgeno)
  d = rbind(d, df)
}
#add a column to d for domesticated/ wild
d$Species <- ifelse(d$Plant == "LA410", "Domesticated",
             ifelse(d$Plant == "LA4355", "Domesticated",
             ifelse(d$Plant == "LA2706", "Domesticated",
             ifelse(d$Plant == "LA4345", "Domesticated",
             ifelse(d$Plant == "LA3475", "Domesticated",
             ifelse(d$Plant == "LA3008", "Domesticated",
                    "Wild"))))))
write.csv(d, "output/lsmeans/BcSlGWAS_lsmeans.csv")

#make a file ready for bigRR
lsmdat <- read.csv("output/lsmeans/BcSlGWAS_lsmeans.csv")
names(lsmdat)
head(lsmdat)
lsmdat <- lsmdat[,c("Igeno", "Estimate", "Plant")]
library(tidyr)
data_wide <- spread(lsmdat, "Plant", "Estimate")
write.csv(data_wide, "output/lsmeans/BcSl_lsmeans_forbigRR.csv")

#now: levene's test within domesticated vs. within wild plant genos
#split dataset by isolate
lsmeans <- read.csv("output/lsmeans/BcSlGWAS_lsmeans.csv")
attach(lsmeans)
out <- split( lsmeans , f = lsmeans$Igeno)
head(out[[1]])
d2 = NULL
#95 levels, one for each isolate
#errors out for 2 isolates in a row after running for 4???
for (i in c(1:95)) {
  print(unique(out[[i]]$Igeno))
  Ltest <- leveneTest(out[[i]]$Estimate, out[[i]]$Species)
  df2 <- as.data.frame(print(Ltest))
  setDT(df2, keep.rownames = T)[]
  df2$Isolate <- unique(out[[i]]$Igeno)
  d2 = rbind(d2, df2)
}
write.csv(d2, "output/lsmeans/BcSlGWAS_LeveneTest.csv")
