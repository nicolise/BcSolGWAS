#Nicole E Soltis
#plotting from bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#Input file: TopSNPs_alltraits_golong.csv
#Output file: TopSNPs_alltraits_Indexed.csv AND TopSNPs_alltraits_wide.csv AND TopSNPs_domest_genesummarized.csv AND TopSNPs_domest_justGenes.csv AND TopSNPs_domestONLY_wide.csv
#Plots: Sl_LesionSize_MAF20_meta.ManhattanPlot.jpg AKA the Manhattan plot with Phenotypes > threshold for Domestication phenos vs. Single geno phenos vs. Both
#AND Sl_LesionSize_MAF20_Domestmeta.ManhattanPlot.jpg AKA the Manhattan plot with DmWoD vs. Domesticated vs. Wild
############################################################################
#Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_golong.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation

# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Segment <- as.numeric(as.character(HEM.plotdata$Segment))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Segment, Pos)), ]

#now make segments line up consecutively
HEM.plotdata$Chrom.Seg <- paste(HEM.plotdata$Chrom, HEM.plotdata$Segment, sep=".")
HEM.plotdata$Chrom.Seg <- as.numeric(HEM.plotdata$Chrom.Seg)

#let's try making the chrom.seg integers so that R isn't confused
unique(HEM.plotdata$Chrom.Seg)
library(plyr)
HEM.plotdata$Chrom.Seg.F <- as.factor(HEM.plotdata$Chrom.Seg)
unique(HEM.plotdata$Chrom.Seg.F)
recode.vars <- data.frame(OGvals = c(1,1.1,2,2.1,2.2,2.3,2.4,2.5, 3, 3.1, 3.2, 4,  4.1, 5,  5.1, 6,  6.1, 6.2, 6.3, 7,  7.1, 7.2, 8,  8.1, 8.2, 9,  9.1, 10, 10.1, 11, 11.1, 12, 12.1, 13, 13.1, 13.2, 14, 14.1, 14.2, 15, 15.1, 15.2, 15.3, 15.4, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 16.11), newvals = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
print(recode.vars)

HEM.plotdata$Chrom.Seg.Int <- recode.vars$newvals[match(HEM.plotdata$Chrom.Seg.F, recode.vars$OGvals)]
unique(HEM.plotdata$Chrom.Seg.Int)

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing

for (i in unique(HEM.plotdata$Chrom.Seg.Int)) {
  print(i)
  #for chromosome 1
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    #changed lastbase+tail to lastbase+max
    lastbase=lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom.Seg.Int==i-1)$Pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Pos+lastbase
  }
  #set ticks to be a list of existing ticks, plus the current Index
  #floor rounds it down to the nearest whole number
  # ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
  
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom.Seg.Int==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

write.csv(HEM.plotdata, "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_Indexed_contigs.csv")

#now go long to wide
library(reshape2)
names(HEM.plotdata)
plotdata.wide <- reshape(HEM.plotdata, idvar=c("Chrom", "Segment", "Pos", "Index"), timevar="Trait", direction = "wide")
names(plotdata.wide)
plotdata.wide$Domest.Num <- rowSums(!is.na(plotdata.wide[,c(6,7,12)]))
plotdata.wide$Plant.Num <- rowSums(!is.na(plotdata.wide[,c(5,8,9,10,11,13,14,15,16,17,18)]))
plotdata.wide$Trait.Num <- rowSums(!is.na(plotdata.wide[,5:18]))
plotdata.wide$Phenos <- "NA"
plotdata.wide$Phenos <- ifelse(plotdata.wide$Plant.Num >= 1 & plotdata.wide$Domest.Num >= 1, "Both", ifelse(plotdata.wide$Plant.Num >= 1, "IndPlants", "Domestication"))
#plotdata.wide$Domest.Phenos <- "NA"
#plotdata.wide$Domest.Phenos <- ifelse(!is.na(plotdata.wide$Effect.Wild) & !is.na(plotdata.wide$Effect.DmWoD), "Wild", ifelse(!is.na) )

write.csv(plotdata.wide, "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_wide_contigs.csv")

#make plots for each phenotype

#create a custom color scale for black and white version
#myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
#colorblind friendly options
#+scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
myColors <- c("#CC79A7", "#009E73", "#E69F00")
names(myColors) <- levels(HEM.plotdata$Phenos)
colScale <- scale_colour_manual(name="Phenotypes", values=myColors)
#names(myColors) <- levels(HEM.plotdata$Chrom)
#colScale <- scale_colour_manual(name = "Chrom",values = myColors)
plotdata.wide <- read.csv("data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_wide_contigs.csv")
names(plotdata.wide)

#indexing for these is shorter than for the WHOLE-SNP figures. So need to redo midpoints list for labeling

#get midpoint positions per chromosome
((max(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index) - min(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index))/2+min(HEM.plotdata[ which(HEM.plotdata$Chrom=='16'),]$Index))
#for full SNPs
#c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687)
#for just this dataset
#c(1778621, 5125545, 7993251, 9994370, 12554828, 15633188, 17342790, 19012547, 20952804, 23185770, 24961413, 25042674, 26494222, 28549075, 29994155, 32411031)

#Fig R8
#plot: all traits
#need colorblind-friendly plot
jpeg("plots/MultiPlot/meta/FigR8_Sl_MAF20_meta.ManhattanPlot.jpg", width=8, height=4, units='in', res=600)
  ggplot(plotdata.wide, aes(x=Index, y=abs(Trait.Num)))+
    theme_bw()+
    colScale+ #add this to modify color palette
    geom_point(aes(color = factor(Phenos)), alpha=1/2)+
    labs(list(y="Number of Phenotypes", x="Chromosome position", title=element_blank()))+
    guides(col = guide_legend(nrow = 8, title="Traits"))+
    scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
   expand_limits(y=0)
dev.off()

#Fig R9
#plot: just domestication
#need colorblind-friendly plot that matches domestication gene Venn Diagram
#myColors <- c("#CC79A7", "#009E73", "#E69F00")

jpeg("plots/MultiPlot/meta/FigR9_Sl_MAF20_Domestmeta.ManhattanPlot.jpg", width=8, height=4, units='in', res=600)
ggplot(plotdata.wide)+
  theme_bw()+
  #    colScale+ #add this for BW version only
  colScale+
  geom_point(aes(x=Index, y=(Effect.Domesticated), color = "Domesticated"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Effect.Wild), color = "Wild"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Effect.DmWoD), color = "(D-W)/D"), alpha=1/2)+
  labs(list(y="Effect size", x="Chromosome position", title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title="Traits"))+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
  expand_limits(y=0)
dev.off()

#keep a subset of rows only if they have Domestication phenotypes
#Phenos column =/= IndPlants
d<-d[!(d$A=="B" & d$E==0),]
plotdata.wide.domest <- plotdata.wide[!(plotdata.wide$Phenos=="IndPlants"),]
write.csv(plotdata.wide.domest, "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_domestONLY_wide.csv")
tapply(HEM.plotdata$Trait, HEM.plotdata$Trait, length)

