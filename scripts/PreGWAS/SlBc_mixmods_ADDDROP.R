#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")

#read in file from SlBc_mixmods.R
ModDat <- read.csv("SlBcDATAFRAME.csv")

#--------------------------------------------------------
#building up a minimal model to check overfitting problem

#final model: 
#lmer(Scale.LS ~ Igeno + Species/PlGenoNm + Igeno:Species/PlGenoNm + Igeno:Species + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (1|Species/PlGenoNm/IndPlant) + AorB + (1|ExpBlock|Igeno) + (1|ExpBlock|PlGenoNm), data = ModDat)
#Igeno, PlGenoNm, Species, ExpBlock, AgFlat, IndPlant, AorB

attach(ModDat)
Lesion.lm.minA1 <- lmer(Scale.LS ~ Igeno + (1|ExpBlock))
LesIso.lsm.minA1 <- lsmeans(Lesion.lm.minA1, "Igeno")

Lesion.lm.minA2 <- lmer(Scale.LS ~ PlGenoNm + (1|ExpBlock))
LesIso.lsm.minA2 <- lsmeans(Lesion.lm.minA2, "PlGenoNm")

Lesion.lm.minA3 <- lmer(Scale.LS ~ Species + (1|ExpBlock))
LesIso.lsm.minA3 <- lsmeans(Lesion.lm.minA3, "Species")

Lesion.lm.minA4 <- lmer(Scale.LS ~ AorB + (1|ExpBlock))
LesIso.lsm.minA4 <- lsmeans(Lesion.lm.minA4, "AorB")

Lesion.lm.minA5 <- lmer(Scale.LS ~ PlGenoNm + Igeno + (1|ExpBlock))

Lesion.lm.minA6 <- lmer(Scale.LS ~ Species + Igeno + (1|ExpBlock))

attach(ModDat)
Lesion.lm.minXX<- lmer(Scale.LS ~ Igeno + Species + (1|IndPlant) + (1|ExpBlock/AgFlat)) 
#first issue detected: PlGenoNm + Species

Lesion.lm.minA14 <- lmer(Scale.LS ~ PlGenoNm + Igeno + AorB +  (1|ExpBlock))

Lesion.lm.minF14 <- lmer(Scale.LS ~ PlGenoNm + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat))

Lesion.lm.min11 <- lmer(Scale.LS ~ Species + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat))

Lesion.lm.min12 <- lmer(Scale.LS ~ PlGenoNm + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat) + (1|IndPlant))

Lesion.lm.min13 <- lmer(Scale.LS ~ Species + Igeno + AorB +  (1|ExpBlock) + (1|AgFlat) + (1|IndPlant))

Lesion.lm.min14 <- lmer(Scale.LS ~ PlGenoNm + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min15 <- lmer(Scale.LS ~ PlGenoNm + Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min16 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|ExpBlock))

Lesion.lm.min20 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#current fullest model that works
Lesion.lm.min21 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min22 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + AorB + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#Warning message: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,: Model is nearly unidentifiable: large eigenvalue ratio - Rescale variables?
Lesion.lm.min17 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock:Igeno) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min18 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock:PlGenoNm) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min23 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + AorB + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#Warning messages: 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,: unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,: Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

Lesion.lm.min23 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|IndPlant) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min24 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|IndPlant) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

Lesion.lm.min24 <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock))

#fixed-effect model matrix is rank deficient so dropping 1 columns / coefficients
Lesion.lm.min07 <- lmer(Scale.LS ~ PlGenoNm + Species + (1|ExpBlock))
Lesion.lm.min09 <- lmer(Scale.LS ~ PlGenoNm + Species + (1|AgFlat))
Lesion.lm.min09 <- lmer(Scale.LS ~ PlGenoNm + Species + (1|IndPlant))

#drop 2 
Lesion.lm.min17 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat) + (ExpBlock:Igeno))

#drop 11
Lesion.lm.min19 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat) + IndPlant)

#drop 12
Lesion.lm.min010 <- lmer(Scale.LS ~ (Species/PlGenoNm) + Igeno + (1|ExpBlock))
Lesion.lm.min09 <- lmer(Scale.LS ~ Species/PlGenoNm + Species + (1|ExpBlock))

#drop 72
Lesion.lm.min22 <- lmer(Scale.LS ~ Species + Igeno + Species:Igeno + (1|ExpBlock) + (1|ExpBlock/AgFlat) + Species/IndPlant)

#drop 792
Lesion.lm.min19 <- lmer(Scale.LS ~ PlGenoNm + Igeno + PlGenoNm:Igeno + (PlGenoNm/IndPlant) + (1|ExpBlock) + (1|ExpBlock/AgFlat))

#myx <- lsmeans(Lesion.lm.min1, trt.vs.ctrl~Igeno, adjust="none")
#trt.vs.ctrl~PlGenoNm|Igeno

#Error in forceSymmetric(2 * solve(IE2)) : 
#error in evaluating the argument 'x' in selecting a method for function #'forceSymmetric': Error in solve.default(IE2) : 
#  system is computationally singular: reciprocal condition number = 8.36832e-17


#non mixed-models
Lesion.lm.min4 <- lmer(Scale.LS ~ (1|IndPlant) + (1|ExpBlock))
Lesion.lm.min6 <- lmer(Scale.LS ~ (1|AgFlat) + (1|ExpBlock))


myx2 <- lsmeans(lsmMod0, pairwise ~ PlGenoNm + Igeno)
write.csv(myx2, "lsmeans_firsttry.csv")

  

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
