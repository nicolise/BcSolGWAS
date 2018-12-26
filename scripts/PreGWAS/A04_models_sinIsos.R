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
ModDat$Igeno <- droplevels(ModDat$Igeno)

library(lme4); library(car); library(lmerTest)

#split dataset by isolate
datIso <- split(ModDat, f = ModDat$Igeno )
head(datIso[[1]]) #100 elements, max. 69 obs per isolate

#check all isolates
for (i in c(1:95)){
  print(unique(datIso[[i]]$Igeno))
}
  
#Using a for loop, iterate over the list of data frames in out[[]]
library(data.table)

#actually 1:95. Test with 1:5
d = NULL
r = NULL
#sink(file='results/output/indiso_100818_fixfx.txt')
for (i in c(1:95)) {
  print(unique(datIso[[i]]$Igeno))
  Mod <- lm(Scale.LS ~ Species + Species/PlGenoNm + ExpBlock, data=datIso[[i]])
  #result <- anova(Mod)
  result2 <- Anova(Mod, type=2)
  #print("Type III ANOVA")
  #print(result)
  #print("Type II ANOVA")
  #print(result2)
  df <- as.data.frame(print(result2))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(datIso[[i]]$Igeno)
  d = rbind(d, df)
}
#sink()
#try it as a tabular output
d.Species <- d[which(d$rn=='Species'),]
d.Plant <- d[which(d$rn=='Species:PlGenoNm'),]
#add fdr-corrected p values
d$fdr.P <- 0
sp_fdr <- p.adjust(d.Species$`Pr(>F)`, method="fdr") 
d[which(d$rn=='Species'),]$fdr.P <- sp_fdr
pl_fdr <- p.adjust(d.Plant$`Pr(>F)`, method="fdr") 
d[which(d$rn=='Species:PlGenoNm'),]$fdr.P <- pl_fdr
write.csv(d, "results/output/indiso_100818_fixfx_tab.csv")

#mixed model
d = NULL
r = NULL
#sink(file='results/output/indiso_100818_mixfx.txt')
for (i in c(1:95)) {
  print(unique(datIso[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ Species + Species/PlGenoNm + (1|ExpBlock), data=datIso[[i]])
  result2 <- Anova(Mod, type=2)
  random <- ranova(Mod)
  #print("Type II ANOVA")
  #print(result2)
  #print(random)
  df <- as.data.frame(print(result2))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(datIso[[i]]$Igeno)
  d = rbind(d, df)
  rf <- as.data.frame(print(random))
  setDT(rf, keep.rownames=T)[]
  rf$Isolate <- unique(datIso[[i]]$Igeno)
  r=rbind(r, rf)
  }
#sink()
#try it as a tabular output
d.Species <- d[which(d$rn=='Species'),]
d.Plant <- d[which(d$rn=='Species:PlGenoNm'),]
#add fdr-corrected p values
#could do this for experiment, too, but we don't care
d$fdr.P <- 0
sp_fdr <- p.adjust(d.Species$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species'),]$fdr.P <- sp_fdr
pl_fdr <- p.adjust(d.Plant$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species:PlGenoNm'),]$fdr.P <- pl_fdr
write.csv(d, "results/output/indiso_100818_mixfx_tab.csv")
write.csv(r, "results/output/indiso_100818_mixfx_rand_tab.csv")

#mixed model v2: testing more terms
#actually 1:95. Test with 1:5
d=NULL
r=NULL
#sink(file='results/output/indiso_100818_mixfx_expanded.txt')
for (i in c(1:95)) {
  print(unique(datIso[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ Species + Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock:Species) + (1|ExpBlock:Species:PlGenoNm), data=datIso[[i]])
  #result <- anova(Mod)
  result2 <- Anova(Mod, type=2)
  random <- ranova(Mod)
  #print("Type II ANOVA")
  #print(result2)
  #print(random)
  df <- as.data.frame(print(result2))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(datIso[[i]]$Igeno)
  d = rbind(d, df)
  rf <- as.data.frame(print(random))
  setDT(rf, keep.rownames=T)[]
  rf$Isolate <- unique(datIso[[i]]$Igeno)
  r=rbind(r, rf)
}
#sink()
#try it as a tabular output
d.Species <- d[which(d$rn=='Species'),]
d.Plant <- d[which(d$rn=='Species:PlGenoNm'),]
#add fdr-corrected p values
d$fdr.P <- 0
sp_fdr <- p.adjust(d.Species$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species'),]$fdr.P <- sp_fdr
pl_fdr <- p.adjust(d.Plant$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species:PlGenoNm'),]$fdr.P <- pl_fdr
#and write out files
write.csv(d, "results/output/indiso_100918_mixfx2_tab.csv")
write.csv(r, "results/output/indiso_100918_mixfx2_rand_tab.csv")

#mixed model v3: testing more terms
#actually 1:95. Test with 1:5
#sink(file='results/output/indiso_100818_mixfx_expandedx2.txt')
d=NULL
r=NULL
for (i in c(1:95)) {
  #print(unique(datIso[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ Species + Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock:Species) + (1|ExpBlock:Species:PlGenoNm) + (1|ExpBlock:PExpRep.x), data=datIso[[i]])
  result2 <- Anova(Mod, type=2)
  random <- ranova(Mod)
  #print("Type II ANOVA")
  #print(result2)
  #print(random)
  df <- as.data.frame(print(result2))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(datIso[[i]]$Igeno)
  d = rbind(d, df)
  rf <- as.data.frame(print(random))
  setDT(rf, keep.rownames=T)[]
  rf$Isolate <- unique(datIso[[i]]$Igeno)
  r=rbind(r, rf)
}
#sink()
#try it as a tabular output
d.Species <- d[which(d$rn=='Species'),]
d.Plant <- d[which(d$rn=='Species:PlGenoNm'),]
#add fdr-corrected p values
d$fdr.P <- 0
sp_fdr <- p.adjust(d.Species$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species'),]$fdr.P <- sp_fdr
pl_fdr <- p.adjust(d.Plant$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species:PlGenoNm'),]$fdr.P <- pl_fdr
#and write out files
write.csv(d, "results/output/indiso_100918_mixfx3_tab.csv")
write.csv(r, "results/output/indiso_100918_mixfx3_rand_tab.csv")

#--------------------------------------------------------------
#old stuff from paper

#original fixed fx model from paper
Mod <- lm(Scale.LS ~ Species/PlGenoNm + ExpBlock, data=out[[i]])
result <- anova(Mod) #this is type III

#adding AgFlat makes model worse
#PExpRep does better than AgFlat as single term
#don't always have enough reps for AorB
#not meaningful (p=1) to add 1|ExpBlock/PlGenoNm
# or to add 1|ExpBlock/Species

head(out[[75]])
#try it as a tabular output
d = NULL
r = NULL
library(data.table)
for (i in c(1:95)) {
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

#base model: ~ Species + Species/PlGenoNm + (1|ExpBlock)
#can't handle adding term: 1|ExpBlock/PExpRep.x OR 1|ExpBlock/AgFlat
for (i in c(1:5)) {
  Mod <- lmer(Scale.LS ~ Species + Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock/AgFlat), data=datIso[[i]])
}