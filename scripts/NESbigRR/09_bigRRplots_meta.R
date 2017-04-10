#Nicole E Soltis
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input files: 
#Output files: 
#Plots: Sl_LesionSize_MAF20_meta.ManhattanPlot.jpg AKA the Manhattan plot with Phenotypes > threshold for Domestication phenos vs. Single geno phenos vs. Both
#AND FigR8_Sl_LesionSize_MAF20_Domestmeta.ManhattanPlot.jpg AKA the Manhattan plot with DmWoD vs. Domesticated vs. Wild 
###########################################################################
#For Plant/ Domest plot
#these source files have all been revised to include positive AND NEGATIVE effect SNPs
Top50SNP.wide.PL <- read.csv("results/Plants_TopSNPs_SegWide.csv")
#combine plant with domestication data
#easiest: take wide data for plant. take wide data for domestication. Omit effects: we only need a count of phenotypes per level.
Top50SNP.wide.DM <- read.csv("results/Domestication_TopSNPs_SegWide_trueMAF.csv")
#remove X column
Top50SNP.wide.DM <- Top50SNP.wide.DM[,2:8]
#count the number of phenotypes
Top50SNP.wide.DM$PlantPhenos <- apply(Top50SNP.wide.DM, 1, function(x) sum(!is.na(x)))
Top50SNP.wide.DM$PlantPhenos <- (Top50SNP.wide.DM$PlantPhenos - 4)
#now remove effect columns
Top50SNP.wide.DM <- Top50SNP.wide.DM[,c(1:4,8)]
Top50SNP.wide.PL <- Top50SNP.wide.PL[,c(2:5,17)]
names(Top50SNP.wide.DM)
names(Top50SNP.wide.PL)
Top50SNP.wide.DM$Trait <- "Domest"
Top50SNP.wide.PL$Trait <- "Plant"
#now stack these
Top50SNP.all <- rbind(Top50SNP.wide.DM, Top50SNP.wide.PL)
#now make it wide
Top50SNP.all.w <- reshape(Top50SNP.all,
                          timevar = "Trait",
                          idvar = c("Chrom", "Segment", "Pos", "Index"),
                          direction = "wide")
#now add a column for BOTH phenos
names(Top50SNP.all.w)
Top50SNP.all.w$PlantPhenos.Both <- (Top50SNP.all.w$PlantPhenos.Domest + Top50SNP.all.w$PlantPhenos.Plant )
#and BOTH does not double-count, because 0 = NA for Plant and Domest. X + NA = NA.

#read in Domestication file for DmWoD vs. Domesticated vs. Wild. This will overwrite the Top50SNP.wide.DM in memory for making the BOTH phenos plot.
Top50SNP.wide.DM <- read.csv("results/Domestication_TopSNPs_SegWide_trueMAF.csv")

#make plots for each phenotype

#create a custom color scale for black and white version
#myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
#colorblind friendly options
#+scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
library(ggplot2)
#domestication is first color, pale green. Wild is lilac.
#colorblind-friendly 3 set
myColors <- c("#2F4F4F", "#9EFA6C", "#AB82FF")

#another 3 set (blue, orange, black)
#domesticated is blue 1C, wild is orange EE, sensitivity is black 05
myColors <- c("#050505", "#1C86EE", "#EE7600")
#names(myColors) <- levels(HEM.plotdata$Phenos)
colScale <- scale_colour_manual(name="Phenotypes", values=myColors)
#names(myColors) <- levels(HEM.plotdata$Chrom)
#colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#get midpoint positions per chromosome
#((max(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index) - min(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index))/2+min(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index))
#for full SNPs
#c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687)


#Fig R7
#plot: all traits
#need colorblind-friendly plot
names(Top50SNP.all.w)
jpeg("paper/plots/ActualPaper/FigR7/Routs/FigR7b_Sl_MAF20_meta.ManhattanPlot2_legend.jpg", width=7.5, height=5, units='in', res=600)
  ggplot(Top50SNP.all.w)+
    theme_bw()+
    colScale+ #add this to modify color palette
    geom_point(aes(x=Index, y=(PlantPhenos.Domest), color = "Domestication"))+
    geom_point(aes(x=Index, y=(PlantPhenos.Plant), color = "Genotype", alpha=1/2))+
    scale_alpha(guide = "none")+
    geom_point(aes(x=Index, y=(PlantPhenos.Both), color = "Both", alpha=1/2))+
    labs(list(y="Number of Phenotypes", x="Chromosome", title=element_blank()))+
    theme(legend.position="right")+
    scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    scale_y_continuous(breaks=c(0,2,4,6,8,10))
   expand_limits(y=0)
dev.off()

#Fig R8
#plot: just domestication
#need colorblind-friendly plot that matches domestication gene Venn Diagram
#domesticated is blue 1C, wild is orange EE, sensitivity is black 05
myColors <- c("#1C86EE", "#050505", "#EE7600")
#names(myColors) <- levels(HEM.plotdata$Phenos)
colScale <- scale_colour_manual(name="Phenotypes", values=myColors)

jpeg("paper/plots/ActualPaper/FigR8/FigR8_SlBc_trueMAF20_domest.ManhattanPlot.jpg", width=7.5, height=5, units='in', res=600)
ggplot(Top50SNP.wide.DM)+
  theme_bw()+
  #    colScale+ #add this for BW version only
  colScale+
  geom_point(aes(x=Index, y=(Effect.Domesticated), color = "Domesticated"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Effect.Wild), color = "Wild"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Effect.DmWoD), color = "Sensitivity"), alpha=1/2)+
  labs(list(y=expression(paste("Estimated Effect Size (",cm^{2},")")), x="Chromosome position", title=element_blank()))+
  #theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  guides(col = guide_legend(nrow = 8, title="Traits"))+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0)
dev.off()
