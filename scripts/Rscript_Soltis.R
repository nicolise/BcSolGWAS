#R script for Nicole's mixed model analysis

BcLesions <- read.csv("SoltisData.csv")
names(BcLesions)
#Scale.LS is my lesion phenotype (area in cm^2)
#Igeno is isolate genotype
#Pgeno is plant genotype (nested within Species)
#Species is plant species
#ExpBlock is which experiment out of 2
#AgFlat is which block (nested within ExpBlock)
#Plant is individual plant
#Leaf is individual leaf (nested within Plant)
#LfinPlant is a renamed variable Leaf (without nested labeling)
#FlatinBlock is a renamed variable AgFlat (without nested labeling)

#fixed effects are Igeno, Pgeno, Species
#random effects are ExpBlock, AgFlat, Plant, Leaf

library(lme4)
fullmod <- lmer(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = BcLesions)
anova(fullmod)
library(car)
Anova(fullmod, type=2)