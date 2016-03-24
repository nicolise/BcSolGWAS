#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#112315
#-------------------------------------------------------------
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
colnames(IsoNm96b)[3] <- "Piso"
#why would I gain data at this step?
LsDat96a2 <- merge(LsDat96a, IsoNm96a, by="Piso")
LsDat96b2 <- merge(LsDat96b, IsoNm96b, by="Piso")
LsDat96a2 <-LsDat96a2[,c(1:14,18,20)]
LsDat96b2 <-LsDat96b2[,c(1:14,18,20)]
SrtDat <- rbind(LsDat96a2, LsDat96b2)

#add rankings for min iso*p and max iso*p
IsoNm2 <- IsoNm[,c(5,8,9)]
SrtDat <- merge(SrtDat, IsoNm2, by="Isolate")

#add a column for plant*iso interaction
SrtDat$PbyI <- paste(SrtDat$PPlant, SrtDat$Isolate, sep='') 

#remove any duplicate entries
SrtDat <- unique(SrtDat)

#prior to building MyPlot, eliminate all P*I with <2 replicates measured
#because can't estimate variance
#SrtDatpl <- SrtDat[SrtDat$n >= 3,]

#----------------------------------------------------------------------

#add a domestication term
names(SrtDat)
unique(SrtDat$PPlant)
#FL NC TX MA KS OR are Domest
#IL UT PA SD NY CA are Wild
SrtDat$Domest <- ifelse(SrtDat$PPlant %in% c("FL","NC","TX","MA","KS","OR"),"Dm", "Wl")

#remove Control Isolates
unique(SrtDat$Isolate)
SrtDat <- SrtDat[SrtDat$Isolate!="Control",]

ModDat <- SrtDat
names(ModDat)
ModDat <- dplyr::select(ModDat, ExpBlock = Pexp, Igeno = Isolate, Pgeno = PPlant, AorB = PInLflt, Leaf = PInLeaf, Plant = PInPlant, AgFlat = PImage, Species = Domest, IsoColor= OrderNJ, IsoMin = minLSbyP, IsoMax = maxLSbyP, matches("."))
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

#add a column: isolates ranked by min lesion size (after avg per plant)
#and a column: isolates ranked by MAX ls (after avg per plant)
# ModDatx1 <- ModDat[complete.cases(ModDat[,21]),]
# ModDatx1 <- ddply(ModDatx1, c("PlGenoNm", "Igeno"), summarise,
#                   mLS   = mean(Scale.LS))
# ModDatx2 <- ddply(ModDatx1, c("Igeno"), summarise,
#                   min  = min(mLS))
# ModDatx3 <- ddply(ModDatx1, c("Igeno"), summarise,
#                   max   = max(mLS))

ddply(ModDat,~PExpRep.x,summarise,mean=mean(Scale.LS))
p <- ggplot(ModDat, aes(factor(PExpRep.x), Lesion.Size))
p + geom_violin() + geom_jitter(height=0)

#-----------------------------------------------------------------
#PCA averaged by plant:iso combination
names(ModDat)
#remove all rows with Scale.LS=0
ModDat$Scale.LS[which(ModDat$Scale.LS==0)] = NA
#log transform the variable
ModDat <- transform(ModDat, log.LS=(log(Scale.LS)))
#remove rows with NA
ModDat1 <- ModDat[complete.cases(ModDat[,19]),]
# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
library(reshape2)
names(ModDat1)
#add a shorter plant geno name
values <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")
ModDat1$Pcol2 <- values[ModDat1$Pgeno]
#can color isolate by IsoMin, IsoMax, or IsoColor (from phylogeny)
#can label plant as PlGenoNm or Pcol2 (short version)
ModDat_wide4 <- dcast(ModDat1, Igeno + IsoMin ~ Species+ Pcol2, value.var="log.LS", mean)
Iso.Color <- ModDat_wide4[,2]
#remove one ExpRep.x at a time... 96a208b is the weird one
logLS4 <- ModDat_wide4[,c(3:14)]

for(i in 1:12){
  logLS4[is.na(logLS4[,i]), i] <- mean(logLS4[,i], na.rm = TRUE)
}

LS.pca4 <- prcomp(logLS4,
                  center = TRUE,
                  scale = TRUE) 
print(LS.pca4)
summary(LS.pca4)
loadings <- LS.pca4$rotation
scores <- LS.pca4$x

##  Hierarchical clustering of PC 3 of scores
scores3 <- scores[,3]
scores2 <- scores[,c(1:6)]

#compare PCs to isolate min, max, mean
scores4 <- data.frame(scores2)
scores4 <- cbind(scores4, ModDat_wide4[,1])
names(scores4)[7] <- "Isolate"
PCbyIso <- merge(scores4, IsoNm, by="Isolate")
names(PCbyIso)
plot(miLSPn~PC1, data=PCbyIso)
plot(maLSPn~PC1, data=PCbyIso)
PCmod <- lm(maLSPn~PC3, data=PCbyIso)
summary(PCmod)
head(PCbyIso)
#PC1 to minLS R2 = 0.3668, p=2.998e-11, 95 df
#PC1 to meanLS R2 = 0.5086, p=2.2e-16, 95 df
#PC1 to maxLS R2= 0.3946, p=3.439e-12, 95 df
#PC2 to minLS R2 = 0, p=0.64, 95 df
#PC2 to meanLS R2 = 0, p=0.40, 95 df
#PC2 to maxLS R2 = 0, p=0.58, 95 df

#histogram of PC3
h <- hist(scores2, xlab="PC3", ylim = c(0,25), breaks=15)
x <- scores2
xfit<-seq(min(x),max(x),length=97) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)

##  Calculate the hierarchical cluster of scores2
scores_dist <- dist (scores2)

##  Calculate the hierarchical clustering of the loadings
hclust_scores <- hclust (scores_dist)
HC.labs <- ModDat_wide4[,1]
plot(hclust_scores, labels=HC.labs, hang=0.1)

##  Show the hierarchical clustering
png (filename="scores-hc.png", width=480, height=480, pointsize=12)
#png (filename="scores-hc.png", width=480, height=480, bg="transparent", pointsize=12)
plot (hclust_scores)
dev.off ()  ##  Close the image
#--------------------------------------------------
#plot it
library(ggbiplot)
#Best Pair = PCs 1 and 3
g <- ggbiplot(LS.pca4, choices = c(1,3), obs.scale = 1, var.scale = 0.1, groups = factor(Iso.Color), ellipse = F, circle = T, varname.size=4, varname.adjust=10, varname.abbrev=T) + theme_bw() 
g <- g + scale_color_discrete(name = '')
g <- g  + scale_x_continuous(limits=c(-5,7))+
  theme(legend.position = 'none')+
  scale_y_continuous(limits=c(-6,5))
print(g)

#another plot option
#library(ggfortify)
#autoplot(LS.pca4, label=T, shape=T) + theme_bw()

#-------------------------------------------------------------
#LA1589 vs LA3475
ModDatmin <- ModDat[ModDat$PlGenoNm=c("LA1589","LA3475"),]
df[df$ID %in% keep, ]
keep <- c("LA1589","LA3475")
ModDatmin <- ModDat[ ModDat$PlGenoNm %in% keep,]
ModDatm_w <- dcast(ModDatmin, Igeno + IsoColor ~ PlGenoNm + Species, value.var="Scale.LS", mean)
names(ModDatm_w)
p <- ggplot(ModDatm_w, aes(LA1589_Wl, LA3475_Dm))
p + geom_point()
p + geom_point(aes(color = factor(IsoColor))) + theme(legend.position = 'none') + geom_smooth(method="lm")

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
p <- ggplot(ModDatm_w, aes(LA1589_Wl, LA3475_Dm))
fit1 <- lm(LA1589_Wl ~ LA3475_Dm, data = ModDatm_w)
ggplotRegression(fit1)

#other pairs
levels(ModDat$PlGenoNm)
keep <- c("LA1589","LA2706")
ModDatminx <- ModDat[ ModDat$PlGenoNm %in% keep,]
ModDatx <- dcast(ModDatminx, Igeno + IsoColor ~ PlGenoNm + Species, value.var="Scale.LS", mean)
names(ModDatx)
fit1 <- lm(LA1589_Wl ~ LA2706_Dm, data = ModDatx)
summary(fit1)
#------------------------------------------------------------------
ModDat_wide5 <- dcast(ModDat1, Igeno + PExpRep.x ~ PlGenoNm + Species, value.var="log.LS", mean)
Exp.Color <- ModDat_wide5[,2]
#remove one ExpRep.x at a time... 96a208b is the weird one
logLS5 <- ModDat_wide5[,c(3:14)]

for(i in 1:12){
  logLS4[is.na(logLS4[,i]), i] <- mean(logLS4[,i], na.rm = TRUE)
}

LS.pca4 <- prcomp(logLS4,
                  center = TRUE,
                  scale = TRUE) 
print(LS.pca4)
summary(LS.pca4)

library(ggbiplot)
#Best Pair = 1,3
g <- ggbiplot(LS.pca4, choices = c(1,3), obs.scale = 1, var.scale = 0.5, groups = factor(Iso.Color), ellipse = F, circle = T, varname.size=4, varname.adjust=5, varname.abbrev=F) + theme_bw() 
g <- g + scale_color_discrete(name = '')
g <- g  + scale_x_continuous(limits=c(-5,7))+
  theme(legend.position = 'none')+
  scale_y_continuous(limits=c(-5,7))
print(g)
#-------------------------------------------------------------------
#PCA with each observation separate
names(ModDat)
#remove all rows with Scale.LS=0
ModDat$Scale.LS[which(ModDat$Scale.LS==0)] = NA
#log transform the variable
ModDat <- transform(ModDat, log.LS=(log(Scale.LS)))
#remove rows with NA
ModDat1 <- ModDat[complete.cases(ModDat[,20]),]
# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
library(reshape2)
names(ModDat1)
ModDat_wide <- dcast(ModDat1, Igeno + PExpRep.x ~ PlGenoNm + Species, value.var="log.LS")
#for 96b208a, FL and IL are both labeled FL! ... fixed
Iso.List <- ModDat_wide[,1]
#remove one ExpRep.x at a time... 96a208b is the weird one
names(ModDat_wide)
ModDat_wide3 <- ModDat_wide[ModDat_wide$PExpRep.x!="96a208b", ]
logLS <- ModDat_wide[,c(3:14)]
logLS3 <- ModDat_wide3[,c(3:14)]
Iso.List3 <- ModDat_wide3[,1]

for(i in 1:12){
  logLS3[is.na(logLS3[,i]), i] <- mean(logLS3[,i], na.rm = TRUE)
}

LS.pca3 <- prcomp(logLS3,
                 center = TRUE,
                 scale = TRUE) 
print(LS.pca3)
summary(LS.pca3)

library(devtools)
library(ggbiplot)
#Best Pair = 2,4
g <- ggbiplot(LS.pca3, choices = c(2,4), obs.scale = 1, var.scale = 1, 
              groups = Iso.List3, ellipse = F, 
              circle = T, varname.size=4, varname.adjust=5, varname.abbrev=F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.position = 'none') 
#+ scale_x_continuous(limits=c(-120,9))+
#  scale_y_continuous(limits=c(-30,40))
print(g)




#remove NAs from table
for(i in 1:12){
  logLS[is.na(logLS[,i]), i] <- mean(logLS[,i], na.rm = TRUE)
}

attach(ModDat1)
LS.pca <- prcomp(logLS,
                 center = TRUE,
                 scale = TRUE) 
print(LS.pca)
summary(LS.pca)

library(devtools)
library(ggbiplot)
#best pair: PC 4,6
g <- ggbiplot(LS.pca, choices = c(4,6), obs.scale = 1, var.scale = 1, 
              groups = Iso.List, ellipse = F, 
              circle = T, varname.size=4, varname.adjust=5, varname.abbrev=F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.position = 'none') 
#+ scale_x_continuous(limits=c(-120,9))+
#  scale_y_continuous(limits=c(-30,40))
print(g)

Exp.List <- ModDat_wide[,2]
g <- ggbiplot(LS.pca, choices = c(1,2), obs.scale = 1, var.scale = 1, 
              groups = Exp.List, ellipse = F, 
              circle = T, varname.size=4, varname.adjust=5, varname.abbrev=F)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.position = 'none') 
#+ scale_x_continuous(limits=c(-120,9))+
#  scale_y_continuous(limits=c(-30,40))
print(g)