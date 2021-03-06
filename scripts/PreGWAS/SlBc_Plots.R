#Preliminary data analysis Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
ModDat <- read.csv("data/SlBcDATAFRAME.csv")
#-------------------------------------------------

#plot my data: scatterplot!
names(ModDat)
attach(ModDat)
names(ModDat)

ModDat$SpLabs <- factor(ModDat$Species, labels = c("Domesticated", "Wild"))
ModDat$SpLabs <- factor(ModDat$SpLabs, levels =c("Wild", "Domesticated"))
#add PlantNum as an integer sorted by mean lesion size
library(plyr)
FigDat3 <- ddply(ModDat, c("PlGenoNm", "Igeno", "Species", "IsoColor"), summarise,
                 mLS   = mean(Scale.LS))
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
MDmeans <- MDmeans[order(MDmeans$Species, MDmeans$mean),] 
MDmeans$PlantNum <- c(1,2,3,4,5,6,7,8,9,10,11,12)
#MDmeans$PlNum <- c(1,2,3,4,5,6,1,2,3,4,5,6)
FigDat3 <- merge(FigDat3, MDmeans, by="PlGenoNm")
FigDat3 <- dplyr::select(FigDat3, Plnum = mean, matches("."))
FigDat3$PlantNum <- as.numeric(FigDat3$PlantNum)
library(ggplot2)
attach(FigDat3)
names(FigDat3)
FigDat3$SpLabs <- factor(FigDat3$Species.x, labels = c("Domesticated", "Wild"))
#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$Igeno)
attach(FigDat3)
FigDat3 <- FigDat3[order(mmLS),]

#make x axis labels that actually work
library(plyr)
FigDat3$Plant.Label <- mapvalues(FigDat3$PlantNum, 
                                 c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                                 c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))
FigDat3$Plant.Lab.Ord <- factor(FigDat3$Plant.Label, levels = c("LA4345", "LA3008", "LA4355", "LA2706", "LA3475", "LA0410", "LA1547", "LA2093", "LA1684", "LA1589", "LA0480", "LA2176"))

# #plot it
# tiff("plots/Sl_LesionSize_rainbowIntx.tiff", width=10, height=6, units='in', res=600)
# ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
#   theme_bw()+
#   geom_line(size=1, aes(color=factor(mmLS), group=factor(Igeno)), show.legend=F, alpha=0.6)+
#   facet_grid(.~SpLabs, scales="free_x")+
#   theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())+
#   labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
# dev.off()
 
#----------------------------------------------------------- 
#scatter plot with subset of isolates colored in

#get list of isolate groups
IsoGroups <- read.csv("data/IsolateGroups.csv")
names(IsoGroups)
IsoGroups <- IsoGroups[,1:7]
FigDat3 <- merge(FigDat3, IsoGroups, by="Igeno")

# #plots
# tiff("plots/Sl_LesionSize_greyIntx.tiff", width=6, height=4, units='in', res=600)
# ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
#   theme_bw()+
#   #order: all, b05.10, high10, intx, low10
#   scale_color_manual(values = c("black", "black",  "black", "black", "black"))+
#   geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.4)+
#   facet_grid(.~SpLabs, scales="free_x")+
#   theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())+
#   labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
# dev.off()

#plot as panels
opar <- par(mfrow=c(1,1))
#first panel: all
p1 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, intx, low10
  scale_color_manual(values = c("black", "black",  "black", "black", "black"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.4)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), strip.background = element_blank(), aspect.ratio=1/1, axis.title.x=element_text(color=NA))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#second panel: b05.10
#mylab <- as.data.frame(c("B05.10", "B05.10"))
p2 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("grey80", "black",  "grey80", "grey80"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), axis.title.y = element_text(color=NA), strip.background = element_blank(), axis.title.x=element_text(color=NA), aspect.ratio=1/1)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#third panel: high subset
p3 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("grey80", "grey80",  "black", "grey80"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), strip.background = element_blank(), strip.text.x=element_text(color=NA), aspect.ratio=1/1, axis.title.x=element_text(color=NA))+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#fourth panel: low subset
p4 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("grey80",  "grey80", "grey80", "black"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), strip.background = element_blank(), axis.title.y = element_text(color=NA), axis.title.x=element_text(color=NA), strip.text.x= element_text(color=NA), aspect.ratio=1/1)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#fifth panel: varying subset
p5 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, Intx, tomato
  scale_color_manual(values = c("grey80",  "black", "grey80"))+
  geom_line(size=0.5, aes(color=factor(GroupsB), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x=element_text(color=NA), axis.title.x=element_text(color=NA), aspect.ratio=1/1)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#sixth panel: tomato subset
p6 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, intx, low10
  scale_color_manual(values = c("grey80",  "grey80", "black"))+
  geom_line(size=0.5, aes(color=factor(GroupsB), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x=element_text(color=NA), axis.title.y = element_text(color=NA), axis.title.x=element_text(color=NA), aspect.ratio=1/1)+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#############################
#the color version
#plot as panels
names(FigDat3)
#plot for species fx
p0 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("turquoise4", "darkorchid3","chartreuse3", "grey80"))+
  geom_line(size=0.5, aes(color=factor(FdrSpp), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1), axis.title.y = element_text(), strip.background = element_blank(), axis.title.x=element_text(color=NA), aspect.ratio=1/1)+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

tiff("plots/poster/Sl_LesionSize_SPintx.tif", width=6, height=3, units='in', res=600)
p0
dev.off()

p9 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("turquoise4", "darkorchid3","chartreuse3", "darkgoldenrod1", "orange", "indianred1", "deeppink", "black", "blue", "coral4", "grey80"))+
  geom_line(size=0.5, aes(color=factor(SpeciesFx), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1), axis.title.y = element_text(), strip.background = element_blank(), axis.title.x=element_text(color=NA), aspect.ratio=1/1)+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#color list: "darkturquoise" "darkorchid3"   "goldenrod1"    "orangered"     "chartreuse2"
#alternates: "indianred1"   "orchid3"    "darkgoldenrod1"  "orange" "chartreuse3"


tiff("plots/poster/Sl_LesionSize_SPintx_noFDR.tif", width=6, height=3, units='in', res=600)
p9
dev.off()

#first panel: all

p1 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+ 
  theme_bw()+ 
 geom_line(size=0.5, aes(color=factor(mmLS), group=factor(Igeno)), show.legend=F, alpha=0.4)+ facet_grid(.~SpLabs, scales="free_x")+ 
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1, color = NA), strip.background = element_blank(), aspect.ratio=1/1, axis.title.x=element_text(color=NA), plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2}))) 

#second panel: b05.10
#mylab <- as.data.frame(c("B05.10", "B05.10"))
p2 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("grey80", "turquoise4",  "grey80", "grey80"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), axis.title.y = element_text(color=NA), strip.background = element_blank(), axis.title.x=element_text(color=NA), aspect.ratio=1/1, plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#third panel: high subset
p3 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("grey80", "grey80",  "orange", "grey80"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), strip.background = element_blank(), strip.text.x=element_text(color=NA), aspect.ratio=1/1, axis.title.x=element_text(color=NA), plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#fourth panel: low subset
p4 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, low10
  scale_color_manual(values = c("grey80",  "grey80", "grey80", "orchid3"))+
  geom_line(size=0.5, aes(color=factor(Group), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle=45, hjust=1, color=NA), strip.background = element_blank(), axis.title.y = element_text(color=NA), axis.title.x=element_text(color=NA), strip.text.x= element_text(color=NA), aspect.ratio=1/1, plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#fifth panel: varying subset
p5 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, Intx, tomato
  scale_color_manual(values = c("grey80",  "chartreuse3", "grey80"))+
  geom_line(size=0.5, aes(color=factor(GroupsB), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x=element_text(color=NA), axis.title.x=element_text(color=NA), aspect.ratio=1/1, plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))

#sixth panel: tomato subset
p6 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
  theme_bw()+
  #order: all, b05.10, high10, intx, low10
  scale_color_manual(values = c("grey80",  "grey80", "indianred1"))+
  geom_line(size=0.5, aes(color=factor(GroupsB), group=factor(Igeno)), show.legend=F, alpha=0.8)+
  facet_grid(.~SpLabs, scales="free_x")+
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), strip.text.x=element_text(color=NA), axis.title.y = element_text(color=NA), axis.title.x=element_text(color=NA), aspect.ratio=1/1, plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
  labs(y=expression(Lesion ~ Area ~ (cm^{2})))
#color list: "darkturquoise" "darkorchid3"   "goldenrod1"    "orangered"     "chartreuse2"
#alternates: "indianred1"   "orchid3"    "darkgoldenrod1"  "orange" "chartreuse3"

tiff("plots/poster/Sl_LesionSize_Intx_a.tif", width=6, height=3, units='in', res=600)
p1 
dev.off()

tiff("plots/poster/Sl_LesionSize_Intx_b.tif", width=6, height=3, units='in', res=600)
p2
dev.off()

tiff("plots/poster/Sl_LesionSize_Intx_c.tif", width=6, height=3, units='in', res=600)
p3 
dev.off()

tiff("plots/poster/Sl_LesionSize_Intx_d.tif", width=6, height=3, units='in', res=600)
p4 
dev.off()

tiff("plots/poster/Sl_LesionSize_Intx_e.tif", width=6, height=3, units='in', res=600)
p5
dev.off()

tiff("plots/poster/Sl_LesionSize_greyIntx_f.tif", width=6, height=3, units='in', res=600)
p6
dev.off()

library(gridExtra)
tiff("plots/paper/DomestRight/Sl_LesionSize_PANELS.tif", width=12, height=8, units='in', res=600)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)
dev.off()
par(opar)
#--------------------------------------------------------------
#plot Domesticated vs. Wild only

library(plyr)
ModDat$Spec.Iso <- paste(ModDat$Species, ModDat$Igeno, sep="")
FigDat4 <- ddply(ModDat, c("Igeno", "Species", "Spec.Iso", "IsoColor"), summarise,
                 mLS   = mean(Scale.LS),
                 sdLS = sd(Scale.LS))
FigDat4$cvLS <- FigDat4$sdLS / FigDat4$mLS
MDmeans <- ddply(ModDat, c("PlGenoNm","Species"), summarise, mean=mean(Scale.LS))
FigDat4$SpLabs <- factor(FigDat4$Species, labels = c("Domesticated", "Wild"))
FigDat4$SpLabs <- factor(FigDat4$SpLabs, levels = c("Wild", "Domesticated"))

#get list of isolate groups
IsoGroups <- read.csv("data/IsolateGroups.csv")
IsoGroups <- IsoGroups[,1:2]
FigDat4 <- merge(FigDat4, IsoGroups, by="Igeno")

FigDatD <- FigDat4[which(FigDat4$Species=="Dm"),]
FigDatW <- FigDat4[which(FigDat4$Species =="Wl"),]

wilcox.test(FigDatD$cvLS, FigDatW$cvLS, paired=T)

#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
#FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$Igeno)
#attach(FigDat3)
#FigDat3 <- FigDat3[order(mmLS),]

tiff("plots/paper/DomestRight/Sl_LesionSize_IntMean_DW.tif", width=6, height=4, units='in', res=600)
ggplot(FigDat4, aes(x = SpLabs, y = mLS, group=factor(Igeno)))+
  theme_bw()+
  geom_line(size=0.5, alpha=0.4, show.legend = F)+
  ylim(0,1.7)+
  theme(text = element_text(size=14), axis.text.x = element_text(), strip.background = element_blank())+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()

#color plot
tiff("plots/poster/Sl_LesionSize_IntMean_DW.tif", width=6, height=4, units='in', res=600)
ggplot(FigDat4, aes(x = SpLabs, y = mLS, group=factor(Igeno)))+
  theme_bw()+
  geom_line(size=0.7, alpha=0.6, show.legend = F, aes(color=factor(mLS)))+
  ylim(0,1.7)+
  theme(text = element_text(size=14), axis.text.x = element_text(), strip.background = element_blank())+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()

ModDat$SpLabs <- factor(ModDat$Species, labels = c("Domesticated", "Wild"))
ModDat2 <- ddply(ModDat, c("Igeno", "SpLabs", "Pgeno"), summarise,
                 mLS   = mean(Scale.LS))
ModDat2$SpLabs <- factor(ModDat2$SpLabs, levels=c("Wild", "Domesticated"))

tiff("plots/poster/Sl_LesionSize_vio_DW.tif", width=4, height=4, units='in', res=600)
ggplot(ModDat2, aes(x = SpLabs, y = mLS))+
  theme_bw() +
  geom_violin(fill = "darkturquoise")+
  geom_boxplot(width = 0.2) +
  ylim(0,1.7) + 
  theme(text = element_text(size=14), axis.text.x = element_text(size = 14, hjust = 1), axis.text.y = element_text(size = 14), aspect.ratio=1.5/1)+
  labs(y = expression (Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_vio_DW.tif", width=4, height=4, units='in', res=600)
ggplot(ModDat2, aes(x = SpLabs, y = mLS))+
  theme_bw() +
  geom_violin(fill = "grey60")+
  geom_boxplot(width = 0.2) +
  ylim(0,1.7) + 
  theme(text = element_text(size=14), axis.text.x = element_text(size = 14, hjust = 1), axis.text.y = element_text(size = 14), aspect.ratio=1.5/1)+
  labs(y = expression (Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()

tiff("plots/Sl_LesionSize_IntCV_DW.tif", width=6, height=4, units='in', res=600)
ggplot(FigDat4, aes(x=SpLabs, y= cvLS, group=factor(Igeno)))+
  theme_bw()+
  geom_line(size=0.5, alpha=0.4, show.legend = F)+
  theme(text = element_text(size=14), axis.text.x = element_text(), strip.background = element_blank())+
  labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()

#---------------------------------------------------------
#rank test of isolates
FigDat5_m <- FigDat4[,c("Igeno", "SpLabs", "IsoColor", "mLS", "Group")]
FigDat5_c <- FigDat4[,c("Igeno", "SpLabs", "cvLS")]
library(tidyr)
FigDat5_w <- spread(FigDat5_m, "SpLabs", "mLS")
library(plyr)
names(FigDat5_w)
FigDat5_w <- plyr::rename(FigDat5_w, c("Domesticated" = "m_Domest", "Wild" = "m_Wild"))
FigDat5_c <- spread(FigDat5_c, "SpLabs", "cvLS")
FigDat5_c <- plyr::rename(FigDat5_c, c("Domesticated" = "cv_Domest", "Wild" = "cv_Wild"))
FigDat5_w <- merge(FigDat5_w, FigDat5_c, by="Igeno")
FigDat5_w$dm <- FigDat5_w$m_Domest - FigDat5_w$m_Wild
FigDat5_w$dcv <- FigDat5_w$cv_Domest - FigDat5_w$cv_Wild

wilcox.test(FigDat5_w$m_Domest, FigDat5_w$m_Wild, paired=TRUE) 
wilcox.test(FigDat5_w$cv_Domest, FigDat5_w$cv_Wild, paired=TRUE) 
#no change for CV

mean(FigDat5_w$m_Domest)
mean(FigDat5_w$m_Wild)
(mean(FigDat5_w$m_Domest) - mean(FigDat5_w$m_Wild))/(mean(FigDat5_w$m_Wild))
#-------------------------------------------------------
#barplot with average lesion size by genotype
#and SE bars
FigDat <- ModDat
names(FigDat)
FigDat2 <- ddply(FigDat, "PlGenoNm", summarise, mLS = mean(Scale.LS))
FigDat3 <- ddply(FigDat, "PlGenoNm", summarise, seLS = se(Scale.LS))

FigDat2 <- ddply(FigDat, c("PlGenoNm", "Species"), summarise,
               N    = length(Scale.LS),
               mean = mean(Scale.LS),
               sd   = sd(Scale.LS),
               se   = sd / sqrt(N))
FigDat2$SpLabs <- factor(FigDat2$Species, labels = c("Domesticated", "Wild"))
limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(FigDat2, aes(x = factor(PlGenoNm), y = mean))+
  geom_bar(stat="identity", fill="dodgerblue3")+
  theme_bw()+
  theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y=expression(Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())+
  geom_errorbar(limits, width=0.25)+
  facet_grid(.~SpLabs, scales="free")


#instead, I'll do this as violin plots
#fill = gray70 for b&w, darkturquoise for color
tiff("plots/paper/DomestRight/Sl_LesionSize_beanplots.tiff", width=6, height=4, units='in', res=600)
ggplot(ModDat, aes(x = factor(PlGenoNm), y = Scale.LS)) + 
  theme_bw() +
  geom_violin(fill = "gray70") + #aes(fill = factor(cyl))
  facet_grid(.~SpLabs, scales="free") +
  geom_boxplot(width = 0.2) +
  theme(text = element_text(size=14), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), axis.text.y = element_text(size = 14), strip.text.x = element_text(size = 14), strip.background=element_blank())+
  labs(y = expression (Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()
#---------------------------------------------------------
#scatterplot: mean of each isolate on each host // domestication level

# mean lesion size per domestication status
names(FigDat4)
FigDat5 <- FigDat4[!duplicated(FigDat4$PbyI), ]
#color lines by positive or negative
unique(unlist(FigDat5$Species))
FigDat5w <- filter(FigDat5, Species == "Wl") 
FigDat5d <- filter(FigDat5, Species == "Dm") 
FigDat5d <-FigDat5d[,c(1,3,4,7)]
FigDat5w <-FigDat5w[,c(1,3,4,7)]
unique(FigDat5w$Igeno)
FigDat5d <- FigDat5d[!duplicated(FigDat5d$Igeno), ]
FigDat5w <- FigDat5w[!duplicated(FigDat5w$Igeno), ]
FigDat6 <- merge(FigDat5w, FigDat5d, by="Igeno")
#get mean difference (Dm - Wl)
head(FigDat6)
FigDat6 <- transform(FigDat6, MeanDiff=(SpMean.y - SpMean.x))
FigDat6$Direction <- ifelse(FigDat6$MeanDiff > 0,1,0)
FigDat7 <- FigDat6[,c(1,9)]
FigDat5 <- merge(FigDat5, FigDat7, by="Igeno")

names(FigDat5)
ggplot(FigDat5, aes(x = Species, y = SpMean, group=Igeno))+
  geom_point()+
  geom_line(size=1, aes(color=factor(Direction)), show.legend=F)+
  theme_bw()+
  labs(y=expression(Mean ~ Lesion ~ Area ~ per ~ Isolate ~ (cm^{2})), x=element_blank())+
  scale_x_discrete(labels=c("Domesticated","Wild"))

#------------------------------------------------------------------
#plot sum of squares for host and for domestication
#for each isolate
library(ggplot2)
SinIsos <- read.csv("output/IsoSpecific/fixANOVA_072716.csv")
names(SinIsos)
ggplot(SinIsos, aes(x=Sum.Sq, fill=rn))+
  geom_density(alpha=0.3)+
  theme_bw()+
  labs(y=expression(Density), x=expression(Sum ~ of ~ Squares))+
  guide_legend()
