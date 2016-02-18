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
library(beanplot); library(ggplot2); library(RColorBrewer); library(plyr)
library(dplyr)

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

#total variance
var(SrtDat$Scale.LS) #equals 0.401204
TotV <- 0.401204

#variance within genotype = environmental variance = Ve
#total variance = Vp = Vg + Ve
#Vg = Vp - Ve
#big H squared = Vg / Vp = (Vp - Ve)/ Vp
#Ve = LsI, LsP, or LsPbyI
#can get negative values if Ve > Vp
# SrtDat3 <- transform(SrtDat, bigHi=((TotV - LsI)/TotV))
# SrtDat3 <- transform(SrtDat3, bigHpl=((TotV - LsP)/TotV))
# SrtDat3 <- transform(SrtDat3, bigHpbyi=((TotV - LsPbyI)/TotV))
# SrtDat3 <- SrtDat3[order("bigHpl"),] 

#prior to building MyPlot, eliminate all P*I with <2 replicates measured
#because can't estimate variance
#SrtDatpl <- SrtDat[SrtDat$n >= 3,]

#----------------------------------------------------------------------
#linear model
#nesting: B within A as A/B or A + A:B
#fixed effects: PInLflt
#random effects: PPlant, Isolate, PInPlant, PInLeaf, Pexp

#add a domestication term
names(SrtDat)
unique(SrtDat$PPlant)
#FL NC TX MA KS OR are Domest
#IL UT PA SD NY CA are Wild
SrtDat$Domest <- ifelse(SrtDat$PPlant %in% c("FL","NC","TX","MA","KS","OR"),"Dm", "Wl")

#remove Control Isolates
unique(SrtDat$Isolate)
SrtDat <- SrtDat[SrtDat$Isolate!="Control",]

#Anova for F-stats
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
#STATISTICS!
#check data structure
xtabs(~ ExpBlock + AgFlat, ModDat)
xtabs(~ Plant + Leaf, ModDat)
#AgFlat is nested within ExpBlock
  #and both are random
#Leaf is nested within Plant
  #and both are random
#And we can consider Pgeno to be nested within Species
#Igeno, Species, Pgeno, AorB are fixed

#library(nlme)
#this model is not quite right
#ModDatF <- na.omit(ModDat)
#lme(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + AorB, random = list(~ 1|ExpBlock/AgFlat + Plant/Leaf), data=ModDatF)

library(lme4)

#split analysis within species
#unique(unlist(ModDat$Species))
#MDwild <- filter(ModDat, Species == "Wl") 
#MDdomest <- filter(ModDat, Species == "Dm") 
#SPmodW <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock/AgFlat) + AorB + (1|ExpBlock) + (1|Plant), data = MDwild)
library(car)
#Anova(SPmodW, type=2)
#anova(SPmodW)
#summary(SPmodW)
#SPmodD <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock/AgFlat) + AorB + (1|ExpBlock) + (1|Plant), data = MDdomest)
#Anova(SPmodD, type=2)
#anova(SPmodD)
#summary(SPmodD)

#Variance output of summary(Mod) gives you SS for the random factors
#rand(Mod) gives Chi-sq  and P values for random factors in packages lmerTest
library(lmerTest)
#rand(SPmodW)
#rand(SPmodD)

#------------------------------------------------------------------------------------------------
#first check normality of Scale.LS
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

#------------------------------------------------------------------------------------------------
#full model
#nesting terms are already included- don't need to add Species as a separate term BUT for random effects (ExpBlock alone) do need a separate term

#error messages: fixed-effect model matrix is rank deficient so dropping 1164 columns / coefficients
#Warning messages:
#1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  : unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  : Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
#Warning message:
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :Model is nearly unidentifiable: large eigenvalue ratio- Rescale variables?

#only 2 ExpBlocks so does it make more sense to include them as fixed effects?
#reasonable to consider AgFlat as a random effect -- have 16 x 3 per exp
#but maybe include a term for "bench"?? random or fixed?
#ExpBlock/Bench/AgFlat

#PlGenoNm is a term nested within Species (but not CODED as if nested within Species = it's not an implicitly nested factor)
#AgFlat IS implicitly nested within ExpBlock -- let's fix this
#non-numeric factors to use: Igeno, PlGenoNm, Species, ExpBlock, AgFlat

#optional to fix: coding of AgFlat so that it is not implicitly nested
fullmod <- lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

#model for lsmeans - Igeno removed
lsmMod <- lmer(Scale.LS ~ Species/PlGenoNm + Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
#model for lsmeans - PlGenoNm removed
lsmMod2 <- lmer(Scale.LS ~ Igeno + Species + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
#model for lsmeans - Igeno and PlGenoNm removed
# lsmMod3 <- lmer(Scale.LS ~ Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
#gives errors: trying expblock as fixed effect (below)
#keeping ExpBlock/AgFlat as random because Flat is random
lsmMod3 <- lmer(Scale.LS ~ Species + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)

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

#removed just igeno * plant interactions (and expblock as fixed)
lsmMod0 <- lmer(Scale.LS ~ Igeno + Species + Species/PlGenoNm + ExpBlock + (1|ExpBlock/AgFlat) + (1|IndPlant) + AorB , data = ModDat)
#error: fixed-effect model matrix is rank deficient so dropping 12 columns / coefficients

library("lsmeans"); 
#also in lmerTest
#example(lsmeans)
modIso.lsm <- lsmeans(lsmMod0, "Igeno")
modPlant.lsm <- lsmeans(lsmMod0, "PlGenoNm")
LesIso.cn <- contrast(modIso.lsm, "trt.vs.ctrlk") #pairwise comparisons vs. UKRazz
LesPlant.cn <- contrast(modPlant.lsm, "trt.vs.ctrlk")

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
#Error in forceSymmetric(2 * solve(IE2)) : 
#error in evaluating the argument 'x' in selecting a method for function #'forceSymmetric': Error in solve.default(IE2) : 
#  system is computationally singular: reciprocal condition number = 8.36832e-17
myx2 <- lsmeans(lsmMod0, pairwise ~ Igeno) # + PlGenoNm?
write.csv(myx2, "lsmeans_firsttry.csv")

#more things to try
# warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
# contrast(warp.lsm, method = "poly", by = "wool")
# ( warp.lsm <- lsmeans (warp.lm,  ~ wool | tension, options = list(estName = "pred.breaks")) )
# (Oats.anal <- lsmeans(Oats.lme, list(poly ~ nitro, pairwise ~ Variety)))
# test(Oats.anal)
# lsmeans(warp.lm, eff ~ wool, adjust = "bonf") 
# lsmeans(warp.lm, eff ~ wool, options = list(adjust = "bonf"))
# summary(wool.lsm, adjust = "bonf")
# contrast(wool.lsm, "eff", adjust = "sidak")
# lsmeans(Oats.lme, ~ nitro, cov.reduce = FALSE, weights = "show.levels")

myx <- lsmeans(lsmMod0, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
#Error in forceSymmetric(2 * solve(IE2)) : 
#error in evaluating the argument 'x' in selecting a method for function #'forceSymmetric': Error in solve.default(IE2) : 
#  system is computationally singular: reciprocal condition number = 8.36832e-17
myx2 <- lsmeans(lsmMod0, pairwise ~ PlGenoNm + Igeno)
write.csv(myx2, "lsmeans_firsttry.csv")

library("lsmeans")
mod3.rg1 <- ref.grid(lsmMod3)
summary(mod3.rg1)
lsmeans(mod3.rg1, "Species")
lsmeans(mod3.rg1, "AorB")

myx <- lsmeans(lsmMod3, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
myx <- lsmeans(fullmod, trt.vs.ctrl~PlGenoNm|Igeno, adjust="none")
#Error in lsmeans.character.ref.grid(object = <S4 object of class "ref.grid">,  : No variable named PlGenoNm in the reference grid
myx2 <- lsmeans(fullmod, pairwise ~ PlGenoNm + Igeno)
write.csv(myx, "lsmeans_firsttry.csv")
write.csv(myx2, "lsmeans_firsttry2.csv")
save.image(file = "0121.RData")

#writing the whole thing out to better understand nesting:
#Does it make sense to include Igeno:PlGenoNm separately? I think not.
Sys.time()
fullmodALT <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + Igeno:Species + Igeno:Species:PlGenoNm + (1|ExpBlock) + (1|ExpBlock:AgFlat) + (1|IndPlant) + AorB , data = ModDat)
save.image(file = "1119.RData") #model is saved as .RData
#use ?load to reload it
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output1117.txt')
Sys.time()
summary(fullmodALT) # the code generating output
Sys.time()
rand(fullmodALT)
Anova(fullmodALT, type=2)
anova(fullmodALT)
Sys.time()
sink()

library(sjPlot)
#
#minimized model to avoid overfitting
Sys.time()
minimod1 <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock:AgFlat) + + AorB , data = ModDat)
save.image(file = "1122.RData") #model is saved as .RData
#use ?load to reload it
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output1122a.txt')
Sys.time()
summary(minimod1) # the code generating output
Sys.time()
rand(minimod1)
Sys.time()
Anova(minimod1, type=2)
Sys.time()
anova(minimod1)
Sys.time()
sink()

sjp.glmer(minimod1, type = "fe.cor")
sjp.glmer(minimod1, type = "re.qq")

#and remove all iso:plant terms
minimod2 <- lmer(Scale.LS ~ Igeno + Species + Species:PlGenoNm + (1|ExpBlock) + (1|ExpBlock:AgFlat) + + AorB , data = ModDat)
save.image(file = "1122b.RData") #model is saved as .RData
#use ?load to reload it
Sys.time()

#copy output to a text file
#other model options- glmm.object? glmm package? nlme?
sink(file='output1122b.txt')
Sys.time()
summary(minimod2) # the code generating output
Sys.time()
rand(minimod2)
Sys.time()
Anova(minimod2, type=2)
Sys.time()
anova(minimod2)
Sys.time()
sink()

sjp.glmer(minimod2, type = "fe.cor")
sjp.glmer(minimod2, type = "re.qq")

#fullmod including Isolate*Exp and Plant*Exp
#add a factor by isolate * experiment
names(ModDat)
ModDat$IbyX <- paste(ModDat$Igeno, ModDat$ExpBlock, sep='') 
sort(unique(ModDat$IbyX))
#remove isolates with no match across experiments:
#94.1, Gallo3
ModDat2 <- ModDat[ModDat$Igeno!="94.1",]
ModDat2 <- ModDat2[ModDat2$Igeno!="Gallo3",]

library(lme4); library(lmerTest)
library(car)
fullmod5 <- lmer(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant) + AorB + (1|ExpBlock) + Igeno:ExpBlock + Species:ExpBlock + ExpBlock:Species/Pgeno + ExpBlock:Igeno:Species/Pgeno, data = ModDat2)
#copy output to a text file
sink(file='output1a.txt')
summary(fullmod5) # the code generating output
sink()

sink(file='output1b.txt')
anova(fullmod5) # the code generating output
sink()
sink(file='output1c.txt')
Anova(fullmod5, type=2)
sink()
sink(file='output1d.txt')
rand(fullmod5)
sink()


fullmod2 <- lmer(Scale.LS ~ Species/Pgeno + Igeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = ModDat)
fullmod
summary(fullmod)
rand(fullmod)


anova(fullmod)
anova(fullmod2)
library(car)
Anova(fullmod, type=2)
qqnorm(residuals(fullmod))
qqline(residuals(fullmod))
plot(fullmod)
print(fullmod)
summary(fullmod)

#try recoding random effects to infer nesting and random
ModDat2 <- ModDat
ModDat2$LfinPlant <- paste(ModDat$Plant, ModDat$Leaf, sep='.')
ModDat2$FlatinBlock <- paste(ModDat$ExpBlock, ModDat$AgFlat, sep='')
names(ModDat2)
#ModDatTmp <- ModDat2[,c("Scale.LS","ExpBlock","Igeno","Pgeno","AorB","Leaf","Plant","AgFlat","Species","LfinPlant","FlatinBlock")]
#write.csv(ModDatTmp,"SoltisData.csv")
unique(ModDat2$FlatinBlock)
unique(ModDat2$LfinPlant)
#and now recode random effects as integers
ModDat2$LfinPlant <- as.numeric(ModDat2$LfinPlant)
ModDat2$FlatinBlock <- as.numeric(factor(ModDat2$FlatinBlock))
ModDat2$Plant <- as.numeric(ModDat2$Plant))
ModDat2$ExpBlock <- as.numeric(factor(ModDat2$ExpBlock))

#and model again
fullmod3 <- lm(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + ExpBlock + FlatinBlock + Plant + LfinPlant + AorB, data = ModDat2)
fullmod
anova(fullmod3)
redmod3 <- lm(Scale.LS ~ Igeno + Species/Pgeno + Species:Igeno + ExpBlock + FlatinBlock + Plant + LfinPlant + AorB, data = ModDat2)
summary(redmod3)
redmod4 <- lm(Scale.LS ~ Igeno + Species/Pgeno + ExpBlock + FlatinBlock + Plant + LfinPlant + AorB, data = ModDat2)
summary(redmod4)
Anova(redmod4, type=2)
redmod5 <- lm(Scale.LS ~ Igeno + Species/Pgeno + ExpBlock + FlatinBlock + Plant + AorB, data = ModDat2)
summary(redmod5)
Anova(redmod5, type=2)

library(car)
Anova(fullmod3, type=2)

aovmod1 <- aov(lm(Lesion.Size ~ Pexp + Pexp/PImage + PPlant*Isolate, data=SrtDat3))
summary(aovmod1)
aovmod2 <- aov(lm(Lesion.Size ~ Pexp + PPlant*Isolate, data=SrtDat3))
summary(aovmod2)
aovmod3 <- aov(lm(Lesion.Size ~ Pexp + Pexp/PImage + PPlant*Isolate + Domest/PPlant/PInPlant, data=SrtDat3))
summary(aovmod3)
#an alternate model
mod1 <- lm(Lesion.Size ~ PPlant*Isolate + PPlant/PInPlant/PInLeaf + PInLflt + Pexp , data=SrtDat3)
#lme for lsmeans
library(lme4)
mod1 <- lmer(Lesion.Size ~ PPlant*Isolate + (1|PImage) + (1|Pexp), data=SrtDat3)
summary(mod1) 
#lsmeans for GWAS
#doesn't work yet
library(lsmeans)
#myLSmeans <- lsmeans(mod1, Lesion.Size ~ PPlant|Isolate, adjust="none")

#make new charts of just single-genotype variance
DFbyIso <- data.frame(VarE=NA, GenFx=NA, Geno=NA)[numeric(0), ]
##how to add column of correct genotype?
IsoVar <- by(SrtDatpl$Scale.LS, SrtDatpl$PIso, var)
newthing <- cbind(IsoVar)
dframe <- data.frame(newthing)
dframe <- transform(dframe, GenFx = "Isolate")
colnames(dframe)[1] <- "VarE"

PlVar <- by(SrtDatpl$Scale.LS, SrtDatpl$PPlant, var)
newtwo <- cbind(PlVar)
dframe2 <- data.frame(newtwo)
dframe2 <- transform(dframe2, GenFx = "Plant")
colnames(dframe2)[1] <- "VarE"

PIVar <- by(SrtDatpl$Scale.LS, SrtDatpl$PbyI, var)
newthree <- cbind(PIVar)
dframe3 <- data.frame(newthree)
dframe3 <- transform(dframe3, GenFx = "PbyI")
colnames(dframe3)[1] <- "VarE"
MyPlot <- rbind(dframe, dframe2, dframe3)


#add total variance column
MyPlot <- transform(MyPlot, TotV = 0.401204)

#add Hsquared column
#variance within genotype = environmental variance = Ve
#total variance = Vp = Vg + Ve
#Vg = Vp - Ve
#big H squared = Vg / Vp = (Vp - Ve)/ Vp
#Ve = LsI, LsP, or LsPbyI
#can get negative values if Ve > Vp
MyPlot <- transform(MyPlot, bigH=((TotV - VarE)/TotV))
beanplot(MyPlot$bigH~MyPlot$GenFx)

library(ggplot2)
ggplot (data = MyPlot, 
  aes(x=GenFx, y=bigH))+
  geom_violin(adjust = 0.5, scale = "width", fill="#E6F598")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  stat_summary(fun.y="median", geom="point")+
  labs(y=expression(Heritability~(H^2)), x=NULL)+
  ylim(c(0,NA))+ #because anything less than 0 is uninteresting
  #and just means that variance for some genos > total variance
  scale_fill_manual(
    values = c("#E6F598","#E6F598", "#E6F598"))+
#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)
geom_jitter(shape=16, position=position_jitter(0.2))


library(vioplot)
#use complete cases only
MyPlotcc <- MyPlot[complete.cases(MyPlot),]
#set H<0 to 0
MyPlotcc$bigH[MyPlotcc$bigH<0] <- 0
x1 <- MyPlotcc$bigH[MyPlotcc$GenFx=="Isolate"]
x2 <- MyPlotcc$bigH[MyPlotcc$GenFx=="Plant"]
x3 <- MyPlotcc$bigH[MyPlotcc$GenFx=="PbyI"]
vioplot(x1, x2, x3, names=c("Isolate", "Plant", "PbyI"), 
        col="green")