#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
#read in file from SlBc_mixmods.csv
ModDat <- read.csv("data/SlBcDATAFRAME.csv")

library(lme4); library(car); library(lmerTest)

#--------------------------------------------------------
#subset columns of interest 
ModDat <- ModDat[,c(2,3,4,5,6,7,8,9,18,20)]
head(ModDat)

#split dataset by isolate
out <- split( ModDat , f = ModDat$Igeno )
head(out[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
#sink(file='output/IsoSpecific/ModelsBYISO_072716rand.txt')
sink(file='output/IsoSpecific/ModelsBYISO_012717FIX.txt')
#skip 13: blank, 59: 94.1, 68: blank, 77: Gallo3, 99: blank
#adding AgFlat makes model worse
#PExpRep does better than AgFlat as single term
#don't always have enough reps for AorB
#not meaningful (p=1) to add 1|ExpBlock/PlGenoNm
# or to add 1|ExpBlock/Species
for (i in c(1:12,14:58,60:67,69:76, 78:98, 100)) {
  print(unique(out[[i]]$Igeno))
  #Mod <- lmer(Scale.LS ~ Species/PlGenoNm + (1|ExpBlock), data=out[[i]])
  Mod <- lm(Scale.LS ~ Species/PlGenoNm + ExpBlock, data=out[[i]])
  result <- anova(Mod)
  #random <- rand(Mod)
  print(result)
  #print(random)
}
sink()

head(out[[75]])
#try it as a tabular output
d = NULL
r = NULL
library(data.table)
#skip 58: 94.1, 75: Gallo3, 99: blank
for (i in c(1:57, 59:74, 76:97)) {
  print(unique(out[[i]]$Igeno))
#  Mod <- lmer(Scale.LS ~ Species + Species/PlGenoNm + (1|ExpBlock), data=out[[i]])
  Mod <- lm(Scale.LS ~ Species/PlGenoNm + ExpBlock, data=out[[i]]) 
  result <- anova(Mod)
  #random <- rand(Mod)
  df <- as.data.frame(print(result))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(out[[i]]$Igeno)
  d = rbind(d, df)
 # rf <- as.data.frame(print(random))
 # setDT(rf, keep.rownames=T)[]
 # rf$Isolate <- unique(out[[i]]$Igeno)
 # r=rbind(r, rf)
}
d.Species <- d[which(d$rn=='Species'),]
d.Plant <- d[which(d$rn=='Species:PlGenoNm'),]
p.adjust(d.Species$`Pr(>F)`, method="fdr") #None with fdr

#write.csv(d, "output/IsoSpecific/fixANOVA_072716.csv")
#write.csv(r, "output/IsoSpecific/randfx_072716.csv")
write.csv(d, "output/IsoSpecific/fixANOVA_012717.csv")

d = NULL
library(data.table)
for (i in c(1:12)) {
  print(unique(out[[i]]$PlantGeno))
  #this one works
  Lesion.lm <- lmer(Scale.LS ~ IsolateID + (1|Exp) + (1|IndPlant) + (1|Exp:IsolateID), data=out[[i]])
  Lesion.lsm <- lsmeans(Lesion.lm, "IsolateID")
  df <- as.data.frame(print(Lesion.lsm))
  setDT(df, keep.rownames = TRUE)[]
  df$Plant <- unique(out[[i]]$PlantGeno)
  d = rbind(d, df)
}
write.csv(d, "output/ModelOutputs/HaLSMeans_062016.csv")

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
MyPvals <- read.csv("data/isoANOVApvals.csv")
names(MyPvals)
p <- MyPvals$pSpecies
MyPvals$pSpFDR <- p.adjust(p, method = "fdr", n = length(p))
#2 sig by FDR
#1 sig, 1 marginally sig by bonf

pSpPl <- MyPvals$pSpPlant
MyPvals$pSpPlFDR <- p.adjust(pSpPl, method = "bonferroni", n = length(pSpPl))
#0 sig by FDR
#0 sig by bonferroni

#what is bp?
bpSp <- MyPvals$BpSpecies
MyPvals$bpSpFDR <- p.adjust(bpSp, method = "fdr", n = length(bpSp))

bpSpPl <- MyPvals$BpSpPlant
MyPvals$bpSpPlFDR <- p.adjust(bpSpPl, method = "fdr", n = length(bpSpPl))

write.csv(MyPvals, "isoANOVAfdr.csv")
#final model: 
#lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock/Igeno) + (1|ExpBlock/PlGenoNm), data = ModDat)
#Igeno, PlGenoNm, Species, ExpBlock, AgFlat, IndPlant, AorB

#-----------------------------------------------------------------------
#single isolates WITHIN species
#--------------------------------------------------------
#subset columns of interest 
ModDat <- ModDat[,c(2,3,4,5,6,7,8,10,13,14,15,16,17,18)]
head(ModDat)

#split ModDat by species
ModDat.W <- subset(ModDat, Species=="Wl")
ModDat.D <- subset(ModDat, Species=="Dm")

#split dataset by isolate
out.W <- split( ModDat.W , f = ModDat.W$Igeno )
head(out.W[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file='ModsBYISO_031816sppW2.txt')
#skip 13: blank, 59: 94.1, 68: blank, 77: Gallo3, 99: blank
for (i in c(1:12,14:58,60:67,69:76, 78:98, 100)) {
  print(unique(out.W[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ PlGenoNm + (1|ExpBlock) + (1|ExpBlock/AgFlat), data=out.W[[i]])
  result <- anova(Mod)
  random <- rand(Mod)
  print(result)
  print(random)
}
sink()

#for D
out.D <- split( ModDat.D , f = ModDat.D$Igeno )
head(out.D[[1]]) #100 elements, max. 69 obs per isolate

#Using a for loop, iterate over the list of data frames in out[[]]
sink(file='ModsBYISO_031816sppD2.txt')
#skip 13: blank, 59: 94.1, 68: blank, 77: Gallo3, 99: blank
for (i in c(1:12,14:58,60:67,69:76, 78:98, 100)) {
  print(unique(out.D[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ PlGenoNm + (1|ExpBlock) + (1|ExpBlock/AgFlat), data=out.D[[i]])
  result <- anova(Mod)
  random <- rand(Mod)
  print(result)
  print(random)
}
sink()

#final model: 
#lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock/Igeno) + (1|ExpBlock/PlGenoNm), data = ModDat)
#Igeno, PlGenoNm, Species, ExpBlock, AgFlat, IndPlant, AorB