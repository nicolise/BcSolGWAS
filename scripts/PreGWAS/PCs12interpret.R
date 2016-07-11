#Nicole E Soltis
#120115
#exploration of groups clustered by PC1 and PC2
#----------------------------------------------------------------------------

#these are the points in small clusters outside the center
#keep <- c(84,25,90,43,29,3,65,93,95,96,75,46,76,50,7,70,47,10,18,23,82,12,87)
#FigDat3min <- FigDat3[ !(FigDat3$IsoColor %in% keep),]
#write.csv(FigDat3min, "lesionsPCsubset2.csv")

NewDat <- read.csv("lesionsPCsubset.csv")
names(NewDat)
NewDat3 <- ddply(NewDat, c("Cluster","Species.x"), summarise, mean=mean(mLS))
NewDat2 <- ddply(NewDat, c("Cluster","Plnum"), summarise, mean=mean(mLS))
p <- ggplot(NewDat2, aes(x=Plnum, y=mean, color=Cluster))
p + geom_point() + geom_smooth(method=lm, se=T)

Mod1 <- lm(mLS~Cluster*Species.x, data=NewDat) #no effect
Mod2 <- lm(mLS~Cluster*PlGenoNm, data=NewDat)
anova(Mod2)
