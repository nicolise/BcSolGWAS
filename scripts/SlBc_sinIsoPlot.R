#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
ModDat <- read.csv("SlBcDATAFRAME.csv")
Isolate <- read.csv("isoPLANTFX.csv")
#-------------------------------------------------
names(ModDat)
attach(ModDat)
library(plyr)

#first add a column for isolate coloring
names(ModDat)
names(Isolate)
Isolate <- Isolate[,c("Isolate","Sort", "IsoSrtName")]
ModDat <- dplyr::select(ModDat, Isolate = Igeno, matches("."))
ModDat <- merge(ModDat, Isolate, by = "Isolate")

#plot my data: scatterplot!
#add PlantNum as an integer sorted by mean lesion size
FigDat3 <- ddply(ModDat, c("PlGenoNm", "IsoSrtName", "Species", "IsoColor", "Sort"), summarise, mLS   = mean(Scale.LS))
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6)
##MDmeans$PlNum <- c(1,2,3,4,5,6,7,8,9,10,11,12)
FigDat3 <- merge(FigDat3, MDmeans, by="PlGenoNm")
#FigDat3 <- dplyr::select(FigDat3, mean = Plnum, matches("."))
FigDat3$PlantNum <- as.numeric(FigDat3$PlNum)
library(ggplot2)
attach(FigDat3)
names(FigDat3)
FigDat3$SpLabs <- factor(FigDat3$Species.x, labels = c("Domesticated", "Wild"))

#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
#scale_x_discrete(breaks=c("1","2","3", "4", "5", "6", "1", "2", "3", "4", "5", "6"),
#                 labels=c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))+

#make x axis labels
library(plyr)  

#make color list
mypalette <- c("grey68", "#FF6666", "#3399FF", "yellowgreen", "mediumpurple2")
FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$IsoSrtName)
attach(FigDat3)
FigDat3 <- FigDat3[order(mmLS),]
library(grid)
ggplot(FigDat3, aes(x = PlNum, y = mLS))+
  geom_point(color="grey68")+
  theme_bw()+
  geom_line(size=1, aes(group=factor(IsoSrtName), color=factor(Sort)), show.legend=F)+
  scale_x_discrete(limit=c("1","2","3", "4", "5", "6", "1", "2", "3", "4", "5", "6"))+
  facet_grid(~SpLabs, scales="free_x", space="free_x")+
  scale_colour_manual(values=mypalette)+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 0, hjust = 0.5), panel.margin=unit(2, "lines"))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())

#try drawing the figure with subsets of data
datasubW <- subset(FigDat3, (FigDat3$SpLabs == "Wild"))
#remove incomplete cases for isolate
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

datasubW <- completeFun(datasubW, "Isolate")
datasubD <- subset(FigDat3, (FigDat3$SpLabs == "Domesticated"))
#plot it myself as two panels
opar <- par(mfrow=c(1,1))
par(mfrow=c(1,2))
p1 <- ggplot(datasubD, aes(x = PlNum, y = mLS))+
  geom_point(color="grey68")+
  theme_bw()+
  geom_line(size=1, aes(group=factor(IsoSrtName), color=factor(Sort)), show.legend=F)+
  scale_x_discrete(limit=c("1","2","3", "4", "5", "6"),
                   labels=c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410"))+
  scale_colour_manual(values=mypalette)+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank(), title = "Domesticated")

p2 <- ggplot(datasubW, aes(x = PlNum, y = mLS))+
  geom_point(color="grey68")+
  theme_bw()+
  geom_line(size=1, aes(group=factor(IsoSrtName), color=factor(Sort)), show.legend=F)+
  scale_x_discrete(limit=c("1","2","3", "4", "5", "6"),
                   labels=c("LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))+
  scale_colour_manual(values=mypalette)+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 90, hjust = -10), axis.ticks=element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+ 
  labs(title = "Wild", x = "")

pdf(file="example.pdf", width=20, height=10)
library(gridExtra)
grid.arrange(p1, p2, ncol=2)
dev.off()