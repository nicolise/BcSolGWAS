#Nicole E Soltis
# single-isolate models for TPC resub.

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

#model with fixed effects only
d=NULL 
r=NULL
for (i in c(1:95)) {
  print(unique(datIso[[i]]$Igeno))
  Mod <- lm(Scale.LS ~ Species + Species/PlGenoNm + ExpBlock, data=datIso[[i]])
  result2 <- Anova(Mod, type=2)
  df <- as.data.frame(print(result2))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(datIso[[i]]$Igeno)
  d = rbind(d, df)
}

#save model results as .csv output
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
for (i in c(1:95)) {
  print(unique(datIso[[i]]$Igeno))
  Mod <- lmer(Scale.LS ~ Species + Species/PlGenoNm + (1|ExpBlock), data=datIso[[i]])
  result2 <- Anova(Mod, type=2)
  random <- ranova(Mod)
  df <- as.data.frame(print(result2))
  setDT(df, keep.rownames = T)[]
  df$Isolate <- unique(datIso[[i]]$Igeno)
  d = rbind(d, df)
  rf <- as.data.frame(print(random))
  setDT(rf, keep.rownames=T)[]
  rf$Isolate <- unique(datIso[[i]]$Igeno)
  r=rbind(r, rf)
}

#save model results as .csv output
d.Species <- d[which(d$rn=='Species'),]
d.Plant <- d[which(d$rn=='Species:PlGenoNm'),]
#add fdr-corrected p values
d$fdr.P <- 0
sp_fdr <- p.adjust(d.Species$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species'),]$fdr.P <- sp_fdr
pl_fdr <- p.adjust(d.Plant$`Pr(>Chisq)`, method="fdr") 
d[which(d$rn=='Species:PlGenoNm'),]$fdr.P <- pl_fdr
write.csv(d, "results/output/indiso_100818_mixfx_tab.csv")
write.csv(r, "results/output/indiso_100818_mixfx_rand_tab.csv")