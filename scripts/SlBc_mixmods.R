#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
MyDat <- read.csv("AllResultsSlBc96.csv")
IsoNm <- read.csv("IsoIDs.csv")
PlantNm <- read.csv("PlantIDs.csv")
AddUnits <- read.csv("AddUnits.csv")

names(MyDat)
#subset to include only lesion size
LsDat <- MyDat[,c(1:11,152)]
#Add pixels per cm conversion
names(LsDat)
names(AddUnits)
LsDat$Sort <- paste(LsDat$PExpRep, LsDat$PImage, sep='') 
AddUnits$Sort <- paste(AddUnits$PExpRep, AddUnits$Image, sep='')
#unique(unlist(LsDat$Sort))
#unique(unlist(AddUnits$Sort))
LsDat2 <- merge(LsDat, AddUnits, by="Sort")
LsDat <- LsDat2[,c(2:13,19)]
#Add column of lesion size in cm squared
LsDat <- transform(LsDat, Scale.LS=(Lesion.Size/(pixelsPcm^2)))
library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr); library(dplyr)

#add actual isolate names
unique(unlist(LsDat$Pexp))
LsDat96a <- filter(LsDat, Pexp == "96a") 
LsDat96b <- filter(LsDat, Pexp == "96b") 
IsoNm96a <- IsoNm
colnames(IsoNm96a)[1] <- "Piso"
IsoNm96b <- IsoNm
colnames(IsoNm96b)[3] <- "Piso"
LsDat96a2 <- merge(LsDat96a, IsoNm96a, by="Piso")
LsDat96b2 <- merge(LsDat96b, IsoNm96b, by="Piso")
LsDat96a2 <-LsDat96a2[,c(1:14,18)] #18 should be Isolate
LsDat96b2 <-LsDat96b2[,c(1:14,18)]
SrtDat <- rbind(LsDat96a2, LsDat96b2)

#add a column for plant*iso interaction
SrtDat$PbyI <- paste(SrtDat$PPlant, SrtDat$Isolate, sep='') 

#remove any duplicate entries
SrtDat <- unique(SrtDat)

#add a domestication term
names(SrtDat)
unique(SrtDat$PPlant)
#FL NC TX MA KS OR are Domest
#IL UT PA SD NY CA are Wild
SrtDat$Domest <- ifelse(SrtDat$PPlant %in% c("FL","NC","TX","MA","KS","OR"),"Dm", "Wl")

#remove Control Isolates
unique(SrtDat$Isolate)
SrtDat <- SrtDat[SrtDat$Isolate!="Control",]

#rename columns as needed
ModDat <- SrtDat
names(ModDat)
ModDat <- dplyr::select(ModDat, ExpBlock = Pexp, Igeno = Isolate, Pgeno = PPlant, AorB = PInLflt, Leaf = PInLeaf, Plant = PInPlant, AgFlat = PImage, Species = Domest, matches("."))

#remove any rows with NA for plant or isolate
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
ModDat <- completeFun(ModDat, c("Pgeno", "Igeno"))

#add a column of correct plant geno names
head(ModDat)
PlGenNum <- dplyr::select(PlantNm, Pgeno = PlantID, matches("."))
PlGenNum$PlGenoNm <- paste("LA", PlGenNum$PlantGeno, sep='') 
PlGenNum <- PlGenNum[c(1:12),c("Pgeno","PlGenoNm")]
ModDat <- merge(ModDat, PlGenNum, by="Pgeno")
names(ModDat)

#Plant is coded as numeric and nested within Pgeno
ModDat$IndPlant <- paste(ModDat$PlGenoNm, ModDat$Plant, sep='.') 
#--------------------------------------------------------
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
Sys.time()
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
Sys.time()
sink(file='output021716.txt')
Sys.time()
#summary(fullmod) # the code generating output
#Sys.time()
rand(fullmod)
Anova(fullmod, type=2)
anova(fullmod)
Sys.time()
sink()

#------------------------------------------------------------------------------
#lsmeans calculations
#check for multicollinearity of variables
library("caret")
ModDat2 <- ModDat[,c("Scale.LS", "Igeno", "Species", "PlGenoNm", "ExpBlock","AgFlat","IndPlant", "AorB")]

#findLinearCombos can only deal with numeric variables, so convert each column
ModDat2$Igeno.f <- as.numeric(factor(ModDat2$Igeno))
ModDat2$Species.f <- as.numeric(factor(ModDat2$Species))
ModDat2$PlGenoNm.f <- as.numeric(factor(ModDat2$PlGenoNm))
ModDat2$ExpBlock.f <- as.numeric(factor(ModDat2$ExpBlock))
ModDat2$AgFlat.f <- as.numeric(factor(ModDat2$AgFlat))
ModDat2$IndPlant.f <- as.numeric(factor(ModDat2$IndPlant))
ModDat2$AorB.f <- as.numeric(factor(ModDat2$AorB))
ModDat3 <- ModDat2[,c("Scale.LS", "Igeno.f", "Species.f", "PlGenoNm.f", "ExpBlock.f","AgFlat.f","IndPlant.f", "AorB.f")]
findLinearCombos(ModDat3) #looks fine

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