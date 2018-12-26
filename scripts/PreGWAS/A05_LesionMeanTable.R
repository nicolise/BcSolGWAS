#Nicole E Soltis

#------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
ModDat <- read.csv("data/preGWAS/SlBc_ModelData.csv")

#no
#mylesion <- read.csv("paper/Submissions/PlantCell/TPC\ revision/Figures/Supplemental\ Figures\ &\ Tables/Supp1_meanlesion.csv")

#remove isolates missing from one experiment
SlSumm <- as.data.frame(with(ModDat, table(Igeno,ExpBlock)))
#missing Exps: Gallo3 (96a), 94.1 (96b)
OgDat <- ModDat
ModDat <- subset(ModDat, Igeno != "Gallo3")
ModDat <- subset(ModDat, Igeno != "94.1")
ModDat <- ModDat[,-c(1)]
ModDat$Igeno <- droplevels(ModDat$Igeno)

myDat.plants <- ddply(ModDat, c("Igeno", "Species", "PlGenoNm"), summarise, 
                      meanLs.cm = mean(Scale.LS, na.rm=TRUE),
                      sd.Ls.cm = sd(Scale.LS),
                      n.Ls = length(Scale.LS),
                      se.Ls.cm = (sd.Ls.cm/(n.Ls^0.5)))

myDat.domest <- ddply(ModDat, c("Igeno", "Species"), summarise, 
                      meanLs.cm = mean(Scale.LS, na.rm=TRUE),
                      sd.Ls.cm = sd(Scale.LS),
                      n.Ls = length(Scale.LS),
                      se.Ls.cm = (sd.Ls.cm/(n.Ls^0.5)))
tbl.plants <- myDat.plants[,c("Igeno","PlGenoNm","meanLs.cm","se.Ls.cm")]
tbl.plants <- reshape(tbl.plants, idvar = "Igeno", timevar = "PlGenoNm", direction = "wide")

tbl.domest <- myDat.domest[,c("Igeno","Species","meanLs.cm","se.Ls.cm")]
tbl.domest <- reshape(tbl.domest, idvar = "Igeno", timevar = "Species", direction = "wide")
setwd("~/Projects/BcSolGWAS/paper")
write.csv(tbl.domest, "Submissions/PlantCell/TPC revision/TableX1_lesiondomest.csv")
write.csv(tbl.plants, "Submissions/PlantCell/TPC revision/TableX1_lesionplants.csv")

#check mean difference D vs. W
head(tbl.domest)
mean(tbl.domest$meanLs.cm.Dm)
mean(tbl.domest$meanLs.cm.Wl)
(0.727 - 0.614) / (0.727) #16% difference
#test mean differences when dropping domestication-associated isolates
#need to come back to this when I decide which are the domest-assoc isos...
t.dom.red <- subset(tbl.domest, Igeno != "Rose")
t.dom.red <- subset(t.dom.red, Igeno != "Fd2")
mean(t.dom.red$meanLs.cm.Dm)
mean(t.dom.red$meanLs.cm.Wl)
(0.723 - 0.619) / (0.723) #14% difference
