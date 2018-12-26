#Nicole E Soltis
#100218

#-------------------------------------------------------------------
#correlation between Bc lesion on tomato Domest, tomato Wild, and Arabidopsis
rm(list=ls())

#get Arabidopsis data 
setwd("~/Projects/BcEudicotGWAS/data/BcAtGWAS/02_csvPrep")
AtPhenos <- read.csv("LSMeanCamLes4Map_FIN.csv")

#get tomato data
setwd("~/Projects/BcSolGWAS/")

ModDat <- read.csv("data/preGWAS/SlBc_ModelData.csv")
ModDat <- subset(ModDat, Igeno != "Gallo3")
ModDat <- subset(ModDat, Igeno != "94.1")
ModDat <- ModDat[,-c(1)]
ModDat$Igeno <- droplevels(ModDat$Igeno)
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
names(AtPhenos)[2] <- "Igeno"
AtLes <- AtPhenos[,c(2,4)]
myDat.corr <- merge(myDat.corr, AtLes, by="Igeno")

# Correlations with significance levels
cor.test(myDat.corr$meanLs.cm.Dm, myDat.corr$Col0.Les)
cor.test(myDat.corr$meanLs.cm.Wl, myDat.corr$Col0.Les)
# library(Hmisc)
# myDat.corr.mx <- myDat.corr[,2:3]
# rcorr(as.matrix(myDat.corr.mx), type="pearson") # type can be pearson or spearman

#scatterplot... could look up old notes on prettier regression lines from Julin's class maybe
library(ggplot2)

my_x_title <- substitute(paste("Mean Lesion Size on ",italic('A. thaliana'), " (cm)"))
plota <- ggplot(myDat.corr, aes(x=Col0.Les, y=meanLs.cm.Wl)) + geom_point()+
  theme_bw()+
  labs(list(y="Mean Lesion Size on Wild Tomato (cm)", x=my_x_title))+
  geom_smooth(method='lm',formula=y~x)

plotb <- ggplot(myDat.corr, aes(x=Col0.Les, y=meanLs.cm.Dm)) + geom_point()+
  theme_bw()+
  labs(list(y="Mean Lesion Size on Domesticated Tomato (cm)", x=my_x_title))+
  geom_smooth(method='lm',formula=y~x)

#define multiplot function below!
jpeg("paper/Submissions/PlantCell/TPC revision/Figures/Supplemental Figures & Tables/Supp3_AtCorrelation.jpg", width=7.5, height=10, units='in', res=600)
#
multiplot(plota, plotb, cols=1)
dev.off()

#-----------------------------------------------------------------------
#make multiplot function
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}