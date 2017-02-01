

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

#lesion size interaction plot
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


#rainbow lesion size interaction plot
# #plot it
# tiff("plots/Sl_LesionSize_rainbowIntx.tiff", width=10, height=6, units='in', res=600)
# ggplot(FigDat3, aes(x = Plant.Lab.Ord, y = mLS))+
#   theme_bw()+
#   geom_line(size=1, aes(color=factor(mmLS), group=factor(Igeno)), show.legend=F, alpha=0.6)+
#   facet_grid(.~SpLabs, scales="free_x")+
#   theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())+
#   labs(y=expression(Lesion ~ Area ~ (cm^{2})), x=element_blank())
# dev.off()