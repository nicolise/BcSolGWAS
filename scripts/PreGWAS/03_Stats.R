#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
ModDat <- read.csv("data/preGWAS/SlBcDATAFRAME.csv")
#get list of isolate groups
IsoGroups <- read.csv("data/preGWAS/IsolateGroups.csv")
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

#but: sample sizes are very unequal when looking at samples from tomato vs. others, so test is underpowered
#pwr.2p2n.test power analysis for 2 proportions, unequal n. enter three of the four quantities (effect size, sample size, significance level, power) and the fourth is calculated. 
library("pwr")
#http://www.statmethods.net/stats/power.html
#Cohen suggests that h values of 0.2, 0.5, and 0.8 represent small, medium, and large effect sizes respectively.
#for IN DOMEST
#h matters a lot an I'm not sure how to set it
#estimate h - for a t test should be difference in means over sqrt( common error variance )
mean(ModDat2D[ModDat2D$Tomato=="y",]$meanLS) - mean(ModDat2D[ModDat2D$Tomato=="n",]$meanLS) / (mean(ModDat2D[ModDat2D$Tomato=="y",]$meanLS) + mean(ModDat2D[ModDat2D$Tomato=="n",]$meanLS) / 2)
count(ModDat2D$Tomato=="y")
pwr.2p2n.test(h=0.27, n1=92, n2=5, power=, sig.level=0.05)
#IN WILD
mean(ModDat2W[ModDat2W$Tomato=="y",]$meanLS) - mean(ModDat2W[ModDat2W$Tomato=="n",]$meanLS) / (mean(ModDat2W[ModDat2W$Tomato=="y",]$meanLS) + mean(ModDat2W[ModDat2W$Tomato=="n",]$meanLS) / 2)
count(ModDat2W$Tomato=="y")
pwr.2p2n.test(h=0.16, n1=92, n2=5, power=, sig.level=0.05)
#ALL HOSTS
mean(ModDat2[ModDat2$Tomato=="y",]$meanLS) - mean(ModDat2[ModDat2$Tomato=="n",]$meanLS) / (mean(ModDat2[ModDat2$Tomato=="y",]$meanLS) + mean(ModDat2[ModDat2$Tomato=="n",]$meanLS) / 2)
count(ModDat2$Tomato=="y")
pwr.2p2n.test(h=0.21, n1=184, n2=10, power=, sig.level=0.05)
#so, ~10% chance of finding a statistically significant difference, even if that difference exists.

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
wilcox.test(ModDat.cv$cvLS ~ ModDat.cv$Species.Num)
wilcox.test(ModDat.cv$mLS ~ ModDat.cv$Species.Num)

#ModDat.cv <- ddply(ModDat, c("PlGenoNm", "Species"), summarise, 
#                   mLS = mean(Scale.LS),
#                   sdLS = sd(Scale.LS))

#wilcoxon signed-rank test between each PAIR of host genotypes
#on isolate means
ModDat.plant <- ddply(ModDat, c("Igeno", "PlGenoNm"), summarise, 
                   mLS = mean(Scale.LS),
                   sdLS = sd(Scale.LS))
ModDat.plant$Plant.Num<- as.numeric(ModDat.plant$PlGenoNm)

#now subset data frame to look at pairs of plant genotypes
outdf <- as.data.frame(NULL)
newdata <- NULL
library(data.table)
for (j in c(1:12)){
  for (i in c(1:12)){
    print(unique(newdata$Plant.Num))
    tryCatch({
      newdata <- ModDat.plant[ ModDat.plant$Plant.Num %in% c(i,j), ]
      newout <- wilcox.test(newdata$mLS ~ newdata$Plant.Num)
      #print(newout)
      df <- as.data.frame(newout$statistic)
      setDT(df, keep.rownames = T)[]
      df$Plant1 <- unique(newdata$Plant.Num)[1]
      df$Plant2 <- unique(newdata$Plant.Num)[2]
      df$sig <- newout$p.value
      outdf = rbind(outdf, df)
      if (i==j) stop("Urgh, these are the same group !")
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #print(newout$statistic)
    #print(newout$p.value)
  }
}
outdf2 <- unique(outdf)
outdf2$fdr.p <- p.adjust(outdf2$sig, method="fdr")
write.csv(outdf2, "paper/plots/ActualPaper/Supp/WilcoxBYPLANT_out.csv")

t.test(ModDat.cv$cvLS, ModDat.cv$Species.Num)
wilcox.test(ModDat.cv$cvLS, ModDat.cv$Species.Num)

#alternate way to do this: F-test for equality of variance of lesion size data between domesticated and wild
dbottle <- lm(Scale.LS ~ Igeno * Species, data=ModDat)
anova(dbottle)
summary(anova(dbottle))
#test for homogeneity of variances 
names(ModDat)
#why summarize by isolate first?
ModDat.vt <- ddply(ModDat, c("Igeno", "Species"), summarise, 
                   mLS = mean(Scale.LS))
ModDat.vt.D <- ModDat.vt[ModDat.vt$Species=="Dm",]
ModDat.vt.W <- ModDat.vt[ModDat.vt$Species=="Wl",]

#this is an F test of the variances, which assumes normality:
var.test(ModDat.vt.D$mLS, ModDat.vt.W$mLS)
#test for normality: both appear to be normally distributed
shapiro.test(ModDat.vt.D$mLS)
shapiro.test(ModDat.vt.W$mLS)

#Levene's F test - page
#Levene's test does not assume normality:
library(car)
leveneTest()
leveneTest(Scale.LS ~ Species, data=ModDat)
leveneTest(Scale.LS ~ Igeno * Species, data=ModDat)
