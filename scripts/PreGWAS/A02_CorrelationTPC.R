#Nicole E Soltis
#100218

#-------------------------------------------------------------------
#correlation between Bc lesion on tomato Domest, tomato Wild, and Arabidopsis
rm(list=ls())

#get Arabidopsis data 
setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/02_csvprep")
AtPhenos <- read.csv("LSMeanCamLes4Map_FIN.csv")

#get tomato data
setwd("~/Projects/BcSolGWAS/")

mylesion <- read.csv("paper/Submissions/PlantCell/TPC\ revision/Figures/Supplemental\ Figures\ &\ Tables/Supp1_meanlesion.csv")

ModDat <- read.csv("data/preGWAS/SlBc_ModelData.csv")
library(plyr)
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

head(myDat.domest)
myDat.corr <- myDat.domest[,c(1,2,3)]
myDat.corr <- reshape(myDat.corr, idvar = "Igeno", timevar = "Species", direction = "wide")
##here, attach At data, matching on Igeno
#may need to rename isolates to match properly

# test correlation -- replace this with At later
cor(myDat.corr$meanLs.cm.Dm, myDat.corr$meanLs.cm.Wl)

# Correlations with significance levels
cor.test(myDat.corr$meanLs.cm.Dm, myDat.corr$meanLs.cm.Wl)
library(Hmisc)
myDat.corr.mx <- myDat.corr[,2:3]
rcorr(as.matrix(myDat.corr.mx), type="pearson") # type can be pearson or spearman

#scatterplot... could look up old notes on prettier regression lines from Julin's class maybe
library(ggplot2)

my_x_title <- substitute(paste("Mean Lesion Size on ",italic('A. thaliana'), " (cm)"))
ggplot(myDat.corr, aes(x=meanLs.cm.Dm, y=meanLs.cm.Wl)) + geom_point()+
  theme_bw()+
  labs(list(y="Mean Lesion Size on Tomato (cm)", x=my_x_title))+
  geom_smooth(method='lm',formula=y~x)
