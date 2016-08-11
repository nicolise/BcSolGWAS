h#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
setwd("~/Documents/GitRepos/BcSolGWAS/data/preGWAS")
ModDat <- read.csv("data/SlBc_ModelData.csv")

#-----------------------------------------------------
#ASSUMPTIONS BEFORE STATISTICS!
#check data structure
xtabs(~ ExpBlock + AgFlat, ModDat)
xtabs(~ Plant + Leaf, ModDat)

#check normality of Scale.LS
require(car)
require(MASS)
ModDat$Scale.LS.t <- ModDat$Scale.LS + 1
#fairly normal
qqp(ModDat$Scale.LS.t, "norm")
#definitely not log-normal
qqp(ModDat$Scale.LS.t, "lnorm")
#data must be integers for the rest
ModDat$Scale.LS.i <- ModDat$Scale.LS*100 + 100
ModDat$Scale.LS.i <- round(ModDat$Scale.LS.i)
#negative binomial looks sort of close, but losing information in the conversion?
nbinom <- fitdistr(ModDat$Scale.LS.i, "Negative Binomial")
qqp(ModDat$Scale.LS.i, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
#poisson is definitely no better than normal
#number of successes in a sequence of iid Bernoulli trials before a specified number of failures (denoted r) occurs: NB(r,p)
#works for overdispersed data (variance >> mean) e.g. counts in ecology
#overdispersion: organisms clumped/ clustered/ aggregated
#here: clumpiness in lesion size categories?
poisson <- fitdistr(ModDat$Scale.LS.i, "Poisson")
qqp(ModDat$Scale.LS.i, "pois", poisson$estimate)

#histogram of data-- goal to compare to poisson / negative binomial
dat <- hist(ModDat$Scale.LS.i, breaks=20)
hist(ModDat$Scale.LS, breaks=20)
dat
mids <- c(110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530)
counts <- c(1813,  850,  782,  649,  553,  471,  353,  242,  184,  141,   93,   61,   30,   19,   13,    6,   10,  3,    6,    3,    1,    2)
histdat <- data.frame(cbind(mids,counts))

#overlay negative binomial
#simplest normalization
ModDat$Nmod <- ModDat$Scale.LS / sum(ModDat$Scale.LS)
#alternative normalization
#temp$Nmod <- temp$N / sqrt(sum(temp$N * temp$N))
histdat$pois <- dpois(histdat$mids, lambda = mean(histdat$counts))
histdat$nbinom <- dnbinom(histdat$counts, mu = mean(histdat$counts), size = 1)
ggplot(histdat, aes(x=mids, y=counts)) +
  geom_histogram(stat="identity", binwidth = 2.5) +
  theme(panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) + 
  geom_line(aes(y = pois), col = "red") + 
  geom_line(aes(y = nbinom), col = "blue")

#----------------------------------------------------------------------
#AgFlat is nested within ExpBlock
#and both are random
#Leaf is nested within Plant
#and both are random
#And we can consider Pgeno to be nested within Species
#Igeno, Species, Pgeno, AorB are fixed

library(lme4); library(car); library(lmerTest)

#try removing isolates missing from one experiment
SlSumm <- as.data.frame(with(ModDat, table(Igeno,ExpBlock)))
#missing Exps: Gallo3 (96a), 94.1 (96b)
OgDat <- ModDat
ModDat <- subset(ModDat, Igeno != "Gallo3")
ModDat <- subset(ModDat, Igeno != "94.1")

#-------------------------------------------------------
#test for effects of agar flat
names(ModDat)
library(plyr)
SummDat.I <- ddply(ModDat, c("Igeno", "Species", "ExpBlock"), summarise,
                 mLS   = mean(Scale.LS),
                 sdLS = sd(Scale.LS))
SummDat.I$cvLS <- SummDat.I$sdLS / SummDat.I$mLS
plot(SummDat.I$Igeno, SummDat.I$cvLS)

SummDat.P <- ddply(ModDat, c("PlGenoNm", "Species", "ExpBlock"), summarise,
                   mLS   = mean(Scale.LS),
                   sdLS = sd(Scale.LS))
SummDat.P$cvLS <- SummDat.P$sdLS / SummDat.P$mLS
plot(SummDat.P$PlGenoNm, SummDat.P$cvLS)

SummDat <- ddply(ModDat, c("PlGenoNm", "Igeno", "Species", "ExpBlock"), summarise,
                   mLS   = mean(Scale.LS),
                   sdLS = sd(Scale.LS))
SummDat$cvLS <- SummDat$sdLS / SummDat$mLS
SummDat$PbyI <- paste(SummDat$PlGenoNm, SummDat$Igeno, sep=".")


flatmod <- lmer(Scale.LS ~ (1|ExpBlock) + (1|ExpBlock/PExpRep.x) + (1|ExpBlock/PExpRep.x/AgFlat) , data = ModDat)
flatmod2 <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock/PExpRep.x) + (1|ExpBlock/PExpRep.x/AgFlat) , data = ModDat)

sink(file='output/flat_effects.txt')
print("flatmod <- lmer(Scale.LS ~ (1|ExpBlock) + (1|ExpBlock/PExpRep.x) + (1|ExpBlock/PExpRep.x/AgFlat) , data = ModDat)")
rand(flatmod)
Anova(flatmod, type=2)
anova(flatmod)

print("flatmod2 <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + (1|ExpBlock) + (1|ExpBlock/PExpRep.x) + (1|ExpBlock/PExpRep.x/AgFlat) , data = ModDat)")
rand(flatmod2)
Anova(flatmod2, type=2)
anova(flatmod2)

sink()

#-----------------------------------------------------------
#GLMS for ANOVA

#Anova(Mod, type=2)
#anova(Mod)
#Variance output of summary(Mod) gives you SS for the random factors
#rand(Mod) gives Chi-sq  and P values for random factors in packages lmerTest
#SPmodD <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock/AgFlat) + AorB + (1|ExpBlock) + (1|Plant), data = MDdomest)

#linear model
#nesting: B within A as A/B or A + A:B
#fixed effects: PInLflt
#random effects: PPlant, Isolate, PInPlant, PInLeaf, Pexp
#ExpBlock and AgFlat as random effects
#but maybe include a random term for "bench"??
#ExpBlock/Bench/AgFlat

#nesting terms are already included- don't need to add Species as a separate term BUT for random effects (ExpBlock alone) do need a separate term
#PlGenoNm is a term nested within Species (but not CODED as if nested within Species = it's not an implicitly nested factor)
#AgFlat IS implicitly nested within ExpBlock -- let's fix this
#non-numeric factors to use: Igeno, PlGenoNm, Species, ExpBlock, AgFlat
#expblock is only 6 terms so I'm going to include it as a fixed effect

#optional to fix: coding of AgFlat so that it is not implicitly nested

#this model I ran for posters/ previous results
#fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#this model is consistent with the lsmeans model:
#Lesion.lm <- lmer(Scale.LS ~ Igeno + Plant/Leaf/AorB + (1|Exp), data=out[[i]])
Sys.time()

#could use without AorB as an option
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) , data = ModDat)

#this is the model I've used
fullmod1 <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant/Leaf) + AorB , data = ModDat)

#fails when including leaf AND AorB AND both exp:I exp:P terms
fullmod2a <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + (1|ExpBlock:Igeno) + (1|ExpBlock:PlGenoNm), data = ModDat)

#trying this: remove indplant. remove agflat.
fullmod2b <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock/Species/PlGenoNm), data = ModDat)

#or this: remove agflat. this one fails.
fullmod3 <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|Species/PlGenoNm/IndPlant) + (1|ExpBlock:Igeno) + (1|ExpBlock:PlGenoNm), data = ModDat)

#this one works
Sys.time()
sink(file='output/fullmod2b_062416.txt')
print("Model: fullmod2b <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:Igeno) + (1|ExpBlock/Species/PlGenoNm), data = ModDat)")
Sys.time()
#summary(fullmod) # the code generating output
rand(fullmod2b)
Anova(fullmod2b, type=2)
anova(fullmod2b)
Sys.time()
sink()

#this one works
Sys.time()
fullmod1 <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat), data = ModDat)
Sys.time()
sink(file='output/fullmod1_072816.txt')
#dropping AorB... don't need it and model fails
#fails to converge with species/plgeno/indplant/leaf
#and with species/ plgeno/indplant
print("fullmod1 <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant), data = ModDat)")
Sys.time()
#summary(fullmod) # the code generating output
rand(fullmod1)
Anova(fullmod1, type=2)
anova(fullmod1)
Sys.time()
sink()


#drop Leaf
fullmod2<- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB , data = ModDat)
#drop IndPlant
fullmod3<- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + AorB , data = ModDat)

#this model is consistent with the lsmeans model, but IndPlant does not make sense without nesting
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#trying by species: LesionWpi.lm and LesionWOpi.lm
#LesionWpi.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
#LesionWOpi.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
