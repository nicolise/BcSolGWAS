#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
setwd("~/Documents/GitRepos/BcSolGWAS/data/preGWAS")
ModDat <- read.csv("SlBc_ModelData.csv")

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
<<<<<<< HEAD
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) , data = ModDat)
=======
fullmod1 <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant/Leaf) + AorB , data = ModDat)
>>>>>>> 4f685a14a122a6ef79dc0b780d992a88a0525041
Sys.time()
sink(file='output_fullmod1_051916.txt')
print("Model: fullmod1 <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant/Leaf) + AorB , data = ModDat)")
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

<<<<<<< HEAD
#this model is consistent with the lsmeans model, but IndPlant does not make sense without nesting
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
=======
>>>>>>> 4f685a14a122a6ef79dc0b780d992a88a0525041

#trying by species: LesionWpi.lm and LesionWOpi.lm
#LesionWpi.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
#LesionWOpi.lm <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#extract p values for each individual predictor variable from anova
anova(lm)$P
# Function to extract the overall ANOVA p-value out of a linear model object
lmp <- function (modelobject) {
  +     if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  +     f <- summary(modelobject)$fstatistic
  +     p <- pf(f[1],f[2],f[3],lower.tail=F)
  +     attributes(p) <- NULL
  +     return(p)
  + }
lmp(lm) #extracts p-value for F-test

#------------------------------------------------------------------------------
#lsmeans calculations

#removed just igeno * plant interactions
#ideally: Igeno, Species, PlGenoNm, AorB are fixed
#ExpBlock, AgFlat, and all their interactions are random
#Igeno with PlGenoNm interaction, PlGenoNm nested within Species, AgFlat nested within ExpBlock, interaction between ExpBlock and PlGenoNm as well as ExpBlock and Igeno
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