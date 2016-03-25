#Make a new scaled data 

library(ggplot2)
library(grid)


#Can impement a quick loop to include each threshold sparately to get the scaled thresholds
CamLes.scaled.plot.data <- cbind(HEM.plotdata[,1:2],abs(apply(HEM.plotdata[,3:12],2,scale)))

ScaledThresh.HEM.99 <- abs(apply(rbind(thresh.HEM$'0.99Thresh',HEM.plotdata[,3:12]),2,scale)[1,])
ScaledThresh.HEM.95 <- abs(apply(rbind(thresh.HEM$'0.95Thresh',HEM.plotdata[,3:12]),2,scale)[1,])


###############################
###Make a HEM plot of a region (Cam and Lesion)
#This code has several dependent objects from the single HEM plotting code above

#PRLM3 (Chr4 Pos 9560135:9565454) in this example.  Just change the Chrom and positions to find another gene.
Reg.min <- 9560135 - 50000
Reg.max <- Reg.min + 1000000


Reg.min <- 13750000 - 50000
Reg.max <- Reg.min + 100000


Region2Plot <- CamLes.scaled.plot.data[CamLes.scaled.plot.data$Chrom == 4 & CamLes.scaled.plot.data$Pos %in% c(Reg.min:Reg.max),]
Region2Plot <- Region2Plot[,-match("pos",colnames(Region2Plot))]
Region2Plot <- melt(Region2Plot, id=c("Chrom","Pos"), variable.name = "Isolate", value.name = "Effect")
Region2Plot$Pheno <- NA
Region2Plot$Pheno[grepl("Cam",Region2Plot$Isolate)] <- "Camalexin"
Region2Plot$Pheno[grepl("Lesion",Region2Plot$Isolate)] <- "Lesion"


Region2Plot$Isolate <- gsub("Camalexin.ng.LesPer._|Lesion.Size.mm2_|.HEM","",Region2Plot$Isolate)

plot.reg <- qplot(Pos,abs(as.numeric(Effect)),data=Region2Plot[Region2Plot$Pheno == "Camalexin",], ylab="|Scaled SNP Effect Estimate|" , colour=factor(Isolate))

#Missing a gene at 13755649-13759740; At4g27550 (wouldn't print)
gene.start = c(13759841,13763464,13765661,13767763,13771185,13772819)
gene.stop = c(13761559,13765064,13766597,13769961,13772594,13777290)

plot.reg <- plot.reg + geom_segment(aes(x = gene.start, y = -2, xend = gene.stop, yend = -2), colour = "red")

#####
###For some reason, I had to use the following code at the terminal to get this to print right.
#xwd > PolyMorphicRegOnChr4.xwd
#convert PolyMorphicRegOnChr4.xwd PolyMorphicRegOnChr4.jpg

jpeg(filename = "PolyMorphicRegionOnChr4.jpg", width = 1500, height = 800)
print(plot2)
dev.off()

















#This works for just about any gene.  Just change the Chr and Position and run the code below


#PAD4
Chr <- 3
Positions <- c(19431371,19434403)

#LYK4
Chr <- 2
Positions <- c(10120242,10122080)

#RLM3
Chr <- 4
Positions <- c(9560135,9565454)

#COI1
Chr <- 2
Positions <- c(16672493,16675748)

#NPR1
Chr <- 1
Positions <- c(23852748,23855566)

#MAM1
Chr <- 5
Positions <- c(7703092,7706896)

#ATHOOK (AT5G54930)
Chr <- 5
Positions <- c(22305734,22307564)

#BIK1
Chr <- 2
Positions <- c(16531688,16533931)

#BOS1
Chr <- 3
Positions <- c(2003393,2006624)

#BOI - invariant in the population
Chr <- 4
Positions <- c(10713561,10714707)


#ANAC055 - invariant in the population
Chr <- 3
Positions <- c(5234617,5236088)

#TGA3
Chr <- 1
Positions <- c(7789348,7792114)

##Heather's loci from Genetics QTL paper
#Cam Only, MSAT3.18
Chr <- 3
Positions <- c(21388082,21388082)

#Cam and Lesion, MSAT4.15
Chr <- 4
Positions <- c(9362675,9362675)




####Run this code to get the plot
#Find positions
Reg.min <- mean(Positions) - 12500
Reg.max <- mean(Positions) + 12500



##Find data to plot
Region2Plot <- CamLes.scaled.plot.data[CamLes.scaled.plot.data$Chrom == Chr & CamLes.scaled.plot.data$Pos %in% c(as.integer(Reg.min):as.integer(Reg.max)),]
Region2Plot <- Region2Plot[,-match("Lesion.Size.mm2_Ctrl.HEM", colnames(Region2Plot))]
Region2Plot <- melt(Region2Plot, id=c("Chrom","Pos"), variable.name = "Isolate", value.name = "Effect")
Region2Plot$Pheno <- NA
Region2Plot$Pheno[grepl("Cam",Region2Plot$Isolate)] <- "Camalexin"
Region2Plot$Pheno[grepl("Lesion",Region2Plot$Isolate)] <- "Lesion"

Region2Plot$Isolate <- gsub("Camalexin.ng.LesPer._|Lesion.Size.mm2_|.HEM","",Region2Plot$Isolate)
Region2Plot$Effect <- abs(Region2Plot$Effect)



##Plot the region

#Make the initial plot
plot.reg <- qplot(Pos,as.numeric(Effect),data=Region2Plot, ylab="|Scaled SNP Effect Estimate|" , shape = Isolate, colour=Pheno)

#Plot the gene region
plot.reg <- plot.reg + geom_segment(aes(x = Positions[1], y = -0.25, xend = Positions[2], yend = -0.25), colour = "darkorange2", size = 2) 


#Plot the colors and facet the plots
my.colours <- c("royalblue3","chartreuse3")


plot.reg <- plot.reg + scale_colour_manual(values=my.colours) + scale_shape_manual(values=c(22,23,21,24,25))
plot.reg <- plot.reg + facet_wrap(~ Pheno, ncol = 1)


#Plot the thresholds
hline.data <- data.frame(z = c(mean(ScaledThresh.HEM.95[1:5]),mean(ScaledThresh.HEM.95[6:10])), Pheno = c("Camalexin","Lesion"))
plot.reg <- plot.reg + geom_hline(aes(yintercept = z), hline.data,colour="red")





#Facets and colours
Cam.subset <- subset(Region2Plot, Pheno == "Camalexin")
Lesion.subset <- subset(Region2Plot, Pheno == "Lesion")
plot.reg + geom_point(data = Cam.subset, colour = "royalblue3", fill="royalblue3",size = 3) + 
	#scale_shape_manual(values=c(22,23,21,24,25)) +
	geom_point(data = Lesion.subset, colour = "chartreuse3", fill="chartreuse3", size = 3) + 
	scale_shape_manual(values=c(22,23,21,24,25)) 



plot.reg + geom_point(data = Cam.subset, size = 1) + 
	scale_shape_manual(values=c(22,23,21,24,25)) +
	geom_point(data = Lesion.subset, size = 1) + 
	scale_shape_manual(values=c(22,23,21,24,25)) 





plot.reg

plot.reg + scale_shape_manual(values=c(22,23,21,24,25)) 
plot.reg + geom_point(fill=c(my.colours,my.colours,"royalblue3",my.colours,my.colours)) + scale_shape_manual(values=c(22,23,21,24,25)) 


plot.reg + scale_colour_discrete(name  ="Phenotype",
                            breaks=c("Apple517", "B05.10","Ctrl","Supersteak","UKRazz"),
                            labels=c("Camalexin", "Lesion")) +
      scale_shape_discrete(name  ="Phenotype",
                            breaks=c("Apple517", "B05.10","Ctrl","Supersteak","UKRazz"),
                            labels=c("Camalexin", "Lesion"))



plot.reg + scale_colour_discrete(name  ="Phenotype",
                            labels=c("Camalexin", "Lesion")) +
      scale_shape_discrete(name  ="Phenotype",
                            labels=c("Camalexin", "Lesion"))








####Genes for the region near the apple517 peak on chromIV
#Missing a gene at 13755649-13759740; At4g27550 (wouldn't print)
gene.start = c(13759841,13763464,13765661,13767763,13771185,13772819)
gene.stop = c(13761559,13765064,13766597,13769961,13772594,13777290)

plot.reg <- plot.reg + geom_segment(aes(x = gene.start, y = -0.001, xend = gene.stop, yend = -.001), colour = "red")

#####
###For some reason, I had to use the following code at the terminal to get this to print right.
#xwd > PolyMorphicRegOnChr4.xwd
#convert PolyMorphicRegOnChr4.xwd PolyMorphicRegOnChr4.jpg

jpeg(filename = "PolyMorphicRegionOnChr4.jpg", width = 1500, height = 800)
print(plot2)
dev.off()








########################################################################

##Blues
"#5986e0"
"#4060a0"
"#334c80"
"#263960"
"#131c30"

##Greens
"#86c959"
"#609040"
"#4c7333"
"#395626"
"#1c2b13"

Blues <- c("#5986e0","#4060a0","#000000","#334c80","#263960")
Greens <- c("#86c959","#609040","#4c7333","#395626")




Cam.subset <- subset(Region2Plot, Pheno == "Camalexin")
Lesion.subset <- subset(Region2Plot, Pheno == "Lesion")



plot.reg + geom_point(data = Cam.subset, size = 3) + 
	scale_shape_manual(data = Cam.subset,values=c(22,23,21,24,25)) +
	scale_colour_manual(data = Cam.subset,values=Blues) +
	scale_fill_manual(data = Cam.subset,values=Blues) +
	geom_point(data = Lesion.subset, size = 3) + 
	scale_shape_manual(data = Lesion.subset,values=c(22,23,21,24,25)) +
	scale_colour_manual(data = Lesion.subset,values=Greens) +
	scale_fill_manual(data = Lesion.subset,values=Greens) 


#########################################################################

plot.reg <- qplot(Pos,as.numeric(Effect),data=Region2Plot, ylab="|Scaled SNP Effect Estimate|", facets = Pheno ~ .)

plot.reg + geom_point(data=Cam.subset,aes(colour = "red"))




hline.data <- data.frame(z = c(mean(ScaledThresh.HEM.95[1:5]),mean(ScaledThresh.HEM.95[6:10])), Pheno = c("Camalexin","Lesion"))
Blues <- c("#007FFF","#006ADD","#000000","#004AAA","#001555")
Greens <- c("#7fff00","#6ADD00","#4AAA00","#155500")

ggplot() + geom_point(data=subset(Region2Plot, Isolate=="Ctrl" & Pheno=="Camalexin"),aes(Pos,Effect), colour = "#000000", fill = "#000000",size = 3, shape = 21) +
	geom_point(data=subset(Region2Plot, Isolate=="Apple517" & Pheno=="Camalexin"),aes(Pos,Effect), colour = Blues[1], fill = Blues[1],size = 3, shape = 22) +
	geom_point(data=subset(Region2Plot, Isolate=="B05.10" & Pheno=="Camalexin"),aes(Pos,Effect), colour = Blues[2], fill = Blues[2],size = 3, shape = 23) +
	geom_point(data=subset(Region2Plot, Isolate=="Supersteak" & Pheno=="Camalexin"),aes(Pos,Effect), colour = Blues[3], fill = Blues[3],size = 3, shape = 24) +
	geom_point(data=subset(Region2Plot, Isolate=="UKRazz" & Pheno=="Camalexin"),aes(Pos,Effect), colour = Blues[4], fill = Blues[4],size = 3, shape = 25) +
	geom_point(data=subset(Region2Plot, Isolate=="Apple517" & Pheno=="Lesion"),aes(Pos,Effect), colour = Greens[1], fill = Greens[1],size = 3, shape = 22) +
	geom_point(data=subset(Region2Plot, Isolate=="B05.10" & Pheno=="Lesion"),aes(Pos,Effect), colour = Greens[2], fill = Greens[2],size = 3, shape = 23) +
	geom_point(data=subset(Region2Plot, Isolate=="Supersteak" & Pheno=="Lesion"),aes(Pos,Effect), colour = Greens[3], fill = Greens[3],size = 3, shape = 24) +
	geom_point(data=subset(Region2Plot, Isolate=="UKRazz" & Pheno=="Lesion"),aes(Pos,Effect), colour = Greens[4], fill = Greens[4],size = 3, shape = 25) +
	geom_segment(data = data.frame(x = rep(Positions,2), y = rep(-0.25,4), Pheno = c(rep("Camalexin",2),rep("Lesion",2))),aes(x = Positions[1], y = -0.25, xend = Positions[2], yend = -0.25), colour = "darkorange2", size = 2) +
	geom_hline(aes(yintercept = z), hline.data,colour="red") +
	facet_grid(Pheno~.)



plot.reg


