#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
ModDat <- read.csv("data/SlBcDATAFRAME.csv")
#get list of isolate groups
IsoGroups <- read.csv("data/IsolateGroups.csv")
IsoGroups <- IsoGroups[,1:4]
#-------------------------------------------------

names(ModDat)
attach(ModDat)
ModDat <- merge(ModDat, IsoGroups, by="Igeno")

#t-test of lesion size by species
t.test(ModDat$Scale.LS ~ ModDat$Species)

#quantiles of lesion size range by species
# function to get quantiles
qfun <- function(x, q = 5) {
  quantile <- cut(x, breaks = quantile(x, probs = 0:q/q), 
                  include.lowest = TRUE, labels = 1:q)
  quantile
}
ModDat[, qq := qfun(Scale.LS), by = Species]
tmp1 <- with(ModDat, split(Scale.LS, Species))
quantile(tmp1$Dm, 0.95)
quantile(tmp1$Dm, 0.05)
quantile(tmp1$Wl, 0.95)
quantile(tmp1$Wl, 0.05)

#t-test of tomato isolates vs. other isolates
#first, take summary data: mean of isolate per species
library(plyr)
ModDat2 <- ddply(ModDat, c("Igeno", "Species", "Tomato"), summarise,
                 meanLS = mean(Scale.LS))
ModDat2D <- ModDat2[ModDat2$Species=="Dm",]
ModDat2W <- ModDat2[ModDat2$Species=="Wl",]
t.test(ModDat2D$meanLS ~ ModDat2D$Tomato)
t.test(ModDat2W$meanLS ~ ModDat2W$Tomato)
t.test(ModDat2$meanLS ~ ModDat2$Tomato)
#can't do multiple observations per isolate (pool D and W directly?)
#does no better when 1 mean per isolate
ModDat3 <- ddply(ModDat, c("Igeno", "Tomato"), summarise,
                 meanLS = mean(Scale.LS))
t.test(ModDat3$meanLS ~ ModDat3$Tomato)

#CV of lesion size between wild vs. domesticated tomato 
#wilcoxon signed-rank test
ModDat.cv <- ddply(ModDat, c("Igeno", "Species"), summarise, 
                   mLS = mean(Scale.LS),
                   sdLS = sd(Scale.LS))
ModDat.cv$cvLS <- ModDat.cv$sdLS/ModDat.cv$mLS
ModDat.cv$Species.Num <- as.numeric(ModDat.cv$Species)
t.test(ModDat.cv$cvLS, ModDat.cv$Species.Num)
wilcox.test(ModDat.cv$cvLS, ModDat.cv$Species.Num)

ModDat.cv <- ddply(ModDat, c("PlGenoNm", "Species"), summarise, 
                   mLS = mean(Scale.LS),
                   sdLS = sd(Scale.LS))

#alternate way to do this: F-test for equality of variance of lesion size data between domesticated and wild
dbottle <- lm(Scale.LS ~ Igeno * Species, data=ModDat)
anova(dbottle)
summary(anova(dbottle))
#test for homogeneity of variances 
names(ModDat)
ModDat.vt <- ddply(ModDat, c("Igeno", "Species"), summarise, 
                   mLS = mean(Scale.LS))
ModDat.vt.D <- ModDat.vt[ModDat.vt$Species=="Dm",]
ModDat.vt.W <- ModDat.vt[ModDat.vt$Species=="Wl",]
var.test(ModDat.vt.D$mLS, ModDat.vt.W$mLS)
library(car)
leveneTest(Scale.LS ~ Species, data=ModDat)
leveneTest(Scale.LS ~ Igeno * Species, data=ModDat)
