rm(list=ls())
setwd("~/PhD/Research/Eudicots/Solanum/Analysis/R")
MyDat <- read.csv("AllResultsSlBc96.csv")
IsoNm <- read.csv("IsoIDs.csv")
PlantNm <- read.csv("PlantIDs.csv")
AddUnits <- read.csv("AddUnits.csv")
?anova

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
colnames(IsoNm96b)[4] <- "Piso"
#why would I gain data at this step?
LsDat96a2 <- merge(LsDat96a, IsoNm96a, by="Piso")
LsDat96b2 <- merge(LsDat96b, IsoNm96b, by="Piso")
LsDat96a2 <-LsDat96a2[,c(1:14,16)]
LsDat96b2 <-LsDat96b2[,c(1:14,17)]
SrtDat <- rbind(LsDat96a2, LsDat96b2)

#add a column for plant*iso interaction
SrtDat$PbyI <- paste(SrtDat$PPlant, SrtDat$Isolate, sep='') 

#First look at variance by geno
#IsoVar <- with(SrtDat, tapply(Scale.LS, Isolate, var))
#PlantVar <- with(SrtDat, tapply(Scale.LS, PPlant, var))
#PbyIVar <- with(SrtDat, tapply(Scale.LS, PbyI, var))

#add column of total phenotypic variance

#use Plyr to add columns of phenotypic variance by P, I, and PbyI
#SrtDat <- ddply(SrtDat, .(Isolate), transform, LsI = var(Scale.LS))
#SrtDat <- ddply(SrtDat, .(PPlant), transform, LsP = var(Scale.LS))
#SrtDat <- ddply(SrtDat, .(Isolate,PPlant), transform, LsPbyI = var(Scale.LS))
#PbyIreps = count(SrtDat, PbyI)
#SrtDat <- merge(SrtDat, PbyIreps, by="PbyI")
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
 #SrtDat3 <- SrtDat3[order("bigHpl"),] 
# write.csv(SrtDat3, "GenoOrders.csv")
# names(LsDat)
# beanplot(LsDat$bigHi)

beanplot(Scale.LS~PPlant, data=SrtDat)
beanplot(Scale.LS~Isolate, data=SrtDat)

#prior to building MyPlot, eliminate all P*I with <2 replicates measured
#because can't estimate variance
#SrtDatpl <- SrtDat[SrtDat$n >= 3,]

#check isolate consistency across experiments
# SrtDat96a <- filter(SrtDat, Pexp == "96a") 
# SrtDat96b <- filter(SrtDat, Pexp == "96b") 
# names(SrtDat96a)
# means96a <- ddply(SrtDat96a, ~PbyI, summarise, mean=mean(Lesion.Size))
# means96b <- ddply(SrtDat96b, ~PbyI, summarise, mean=mean(Lesion.Size))
# means96 <- merge(means96a, means96b, by="PbyI")
# plot(means96$mean.x, means96$mean.y)
# attach(means96)
# myline <- lm(mean.x~mean.y)
# par(cex=.8)
# plot(mean.x, mean.y)
# abline(myline)

#check isolate consistency across replicates
# SrtDat208a <- filter(SrtDat96a, PExpRep.x == "96a208a") 
# SrtDat208b <- filter(SrtDat96a, PExpRep.x == "96a208a") 
# names(SrtDat208a)
# means208a <- ddply(SrtDat208a, ~PbyI, summarise, mean=mean(Lesion.Size))
# means208b <- ddply(SrtDat208b, ~PbyI, summarise, mean=mean(Lesion.Size))
# means208 <- merge(means208a, means208b, by="PbyI")
# plot(means208$mean.x, means208$mean.y)

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

#-------------------------------------------------
#plot my data: scatterplot!
names(ModDat)
attach(ModDat)
#add PlantNum as an integer sorted by mean lesion size
FigDat3 <- ddply(ModDat, c("Pgeno", "Igeno", "Species"), summarise,
                 mLS   = mean(Scale.LS))
MDmeans <- ddply(ModDat, c("Pgeno","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlNum <- 1:nrow(MDmeans)
FigDat3 <- merge(FigDat3, MDmeans, by="Pgeno")
#FigDat <- dplyr::select(FigDat, Plnum = mean, matches("."))
#FigDat$PlantNum <- as.numeric(FigDat$PlNum)
library(ggplot2)
attach(FigDat3)
names(FigDat3)
FigDat3$SpLabs <- factor(FigDat3$Species.x, labels = c("Domesticated", "Wild"))
ggplot(FigDat3, aes(x = PlNum, y = mLS))+
  geom_point()+
  facet_grid(.~SpLabs, scales="free")+
  theme_bw()+
  geom_line(size=1, aes(color=factor(Igeno)), show_guide=F)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x="Plant Genotype")+
  theme(axis.text.x = element_blank())
#+geom_smooth(aes(group = 2), size = 2, method = "lm", se = T)

#-------------------------------------------------------
#barplot with average lesion size by genotype
#and SE bars
names(FigDat)
FigDat2 <- ddply(FigDat, "Pgeno", summarise, mLS = mean(Scale.LS))
FigDat3 <- ddply(FigDat, "Pgeno", summarise, seLS = se(Scale.LS))

FigDat2 <- ddply(FigDat, c("Pgeno", "Species"), summarise,
               N    = length(Scale.LS),
               mean = mean(Scale.LS),
               sd   = sd(Scale.LS),
               se   = sd / sqrt(N))
FigDat2$SpLabs <- factor(FigDat2$Species, labels = c("Domesticated", "Wild"))
limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(FigDat2, aes(x = factor(Pgeno), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  labs(y=expression(Mean ~ Lesion ~ Area ~ (cm^{2})), x="Plant Genotype")+
  facet_grid(.~SpLabs, scales="free")+
  theme(axis.text.x = element_blank())+
  geom_errorbar(limits)

#---------------------------------------------------------
#scatterplot: mean of each isolate on each host // domestication level
ModDat$SbyI <- paste(ModDat$Species, ModDat$Igeno, sep='') 
FigDat4 <- ddply(ModDat, c("Pgeno", "Igeno", "Species", "SbyI","PbyI"), summarise,PlMean   = mean(Scale.LS))
FDmeans <- ddply(ModDat, c("SbyI"), summarise, SpMean=mean(Scale.LS))
FigDat4 <- merge(FigDat4, FDmeans, by="SbyI")
attach(FigDat4)
ggplot(FigDat4, aes(x = PlMean, y = SpMean, color=Species, group=Pgeno))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_line(size=1, aes(color=factor(Igeno)), show_guide=F)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x="Plant Genotype")+
  theme(axis.text.x = element_blank())

#or more simply...
# mean lesion size per domestication status
names(FigDat4)
FigDat5 <- FigDat4[!duplicated(FigDat4$PbyI), ]
ggplot(FigDat4, aes(x = SpMean, y = Species, group=Igeno))+
  geom_point()+
  geom_line()+
  theme_bw()

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
unique(unlist(ModDat$Species))
MDwild <- filter(ModDat, Species == "Wl") 
MDdomest <- filter(ModDat, Species == "Dm") 
SPmodW <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = MDwild)
library(car)
Anova(SPmodW, type=2)
SPmodD <- lmer(Scale.LS ~ Igeno + Pgeno + Igeno:Pgeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = MDdomest)
Anova(SPmodD, type=2)
anova(SPmodD)
anova(SPmodW)
summary(SPmodD)
#Variance output of summary(Mod) gives you SS for the random factors
#rand(Mod) gives Chi-sq  and P values for random factors in packages lmerTest
library(lmerTest)
rand(SPmodW)

fullmod <- lmer(Scale.LS ~ Igeno + Species/Pgeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = ModDat)
fullmod2 <- lmer(Scale.LS ~ Species/Pgeno + Igeno + Igeno:Species/Pgeno + Species:Igeno + (1|ExpBlock/AgFlat) + (1|Plant/Leaf) + AorB, data = ModDat)
fullmod
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