#Nicole E Soltis
#Phenotype plotting Solanum spp. x Botrytis cinerea lesions
#101315
#-------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS")
FigDat3 <- read.csv("data/BcSolGWAS_PhenotypePlotData_FigDat3.csv")
ModDat <- read.csv("data/BcSolGWAS_PhenotypePlotData_ModDat.csv")
FigDat4 <- read.csv("data/BcSolGWAS_PhenotypePlotData_FigDat4.csv")
#-------------------------------------------------
#-------------------------------------------------------
#after loading dataframe have to reorder domestication levels to plot
ModDat$SpLabs <- factor(ModDat$SpLabs, levels = c("Wild", "Domesticated"))
FigDat3$SpLabs <- factor(FigDat3$SpLabs, levels = c("Wild", "Domesticated"))
FigDat4$SpLabs <- factor(FigDat4$SpLabs, levels = c("Wild", "Domesticated"))
library(ggplot2)
#domestication 2 set
#domestication is first color, pale green. Wild is lilac.
myColors <- c("#9EFA6C", "#AB82FF")
names(myColors) <- levels(ModDat$Species)
colScale <- scale_fill_manual(name = "Species",values = myColors)

#Violin plots with mean per genotype
names(ModDat)

tiff("plots/paper/DomestRight/Sl_LesionSize_beanplots.tiff", width=7.5, height=5.4, units='in', res=600)
ggplot(ModDat, aes(x = factor(PlGenoNm), y = Scale.LS)) + 
  theme_bw() +
  guides(fill=F)+
  geom_violin(aes(fill=factor(Species)))+
  colScale+
  facet_grid(.~SpLabs, scales="free") +
  geom_boxplot(width = 0.2) +
  theme(text = element_text(size=14), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), axis.text.y = element_text(size = 14), strip.text.x = element_text(size = 14), strip.background=element_blank())+
  labs(y = expression (Mean ~ Lesion ~ Area ~ (cm^{2})), x=element_blank())
dev.off()
#--------------------------------------------------------------
#Domesticated Vs. Wild ONLY

#add a column of mmLS (mean of mean lesion size) per isolate
#sort dataframe by mmLS 
#then color by the new factor mmLS
#FigDat3$mmLS <- ave(FigDat3$mLS, FigDat3$Igeno)
#attach(FigDat3)
#FigDat3 <- FigDat3[order(mmLS),]

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
library(plyr)
ModDat2 <- ddply(ModDat, c("Igeno", "SpLabs", "Pgeno"), summarise,
                 mLS   = mean(Scale.LS))
ModDat2$SpLabs <- factor(ModDat2$SpLabs, levels=c("Wild", "Domesticated"))

#can draw half-violins if I download vioplot2, but I'm lazy
tiff("plots/poster/Sl_LesionSize_vio_DW.tif", width=6, height=4, units='in', res=600)
ggplot(ModDat2, aes(x = SpLabs, y = mLS))+
  theme_bw() +
  colScale+
  geom_violin(aes(fill=factor(SpLabs)))+
  guides(fill=F)+
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
  guide_legend(title=NULL)
#------------------------------------------------------------------------
#panel plots with subset of isolates colored in

#really can't use facet_grid with margins, so have to fix this by hand, plot.margin=unit(c(0,0,0,0), "cm")
#############################
#the color version
#plot as panels
names(FigDat3)

#color list: "darkturquoise" "darkorchid3"   "goldenrod1"    "orangered"     "chartreuse2"
#alternates: "indianred1"   "orchid3"    "darkgoldenrod1"  "orange" "chartreuse3"

#first panel: all

p1 <- ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+ 
  theme_bw()+ 
 geom_line(size=0.5, aes(color=factor(mmLS), group=factor(Igeno)), show.legend=F, alpha=0.4)+ facet_grid(.~SpLabs, scales="free_x")+ 
  theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(), aspect.ratio=1/1, axis.title.x=element_text(), plot.margin=unit(c(-1,-2,-1,-2),"lines"))+ 
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
############################
tiff("plots/poster/DomestRight/Sl_LesionSize_Intx_a.tif", width=6, height=3, units='in', res=600)
p1 
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_Intx_b.tif", width=6, height=3, units='in', res=600)
p2
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_Intx_c.tif", width=6, height=3, units='in', res=600)
p3 
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_Intx_d.tif", width=6, height=3, units='in', res=600)
p4 
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_Intx_e.tif", width=6, height=3, units='in', res=600)
p5
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_greyIntx_f.tif", width=6, height=3, units='in', res=600)
p6
dev.off()

tiff("plots/paper/DomestRight/Sl_LesionSize_greyIntx_g.tif", width=6, height=3, units='in', res=600)
p7
dev.off()

library(gridExtra)
tiff("plots/paper/DomestRight/Sl_LesionSize_PANELS.tif", width=12, height=8, units='in', res=600)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)
dev.off()
par(opar)

###################################
#extra plotting scripts
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


tiff("plots/poster/Sl_LesionSize_SPintx_noFDR.tif", width=6, height=3, units='in', res=600)
p9
dev.off()