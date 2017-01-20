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
#AND Venn Diagrams for SNP overlaps between these lists
############################################################################
#Plotting the HEM results

#NEED TO CHECK CONTIGS FOR THIS

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_golong.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation

#Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Segment <- as.numeric(as.character(HEM.plotdata$Segment))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#sort dataframe rows in order of Chrom, then Pos
HEM.plotdata2 <- HEM.plotdata[with(HEM.plotdata, order(Chrom, Segment, Pos)), ]

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
#still need to fix plotting so that CONTIGS are sequential (HEM.plotdata$Segment)
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    #changed lastbase+tail to lastbase+max
    lastbase=+lastbase+max(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

write.csv(HEM.plotdata, "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_Indexed.csv")

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

write.csv(plotdata.wide, "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_wide.csv")

#make plots for each phenotype

#create a custom color scale for black and white version
#myColors <- c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
#names(myColors) <- levels(HEM.plotdata$Chrom)
#colScale <- scale_colour_manual(name = "Chrom",values = myColors)
plotdata.wide <- read.csv("data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_alltraits_wide.csv")
names(plotdata.wide)
#without the loop: all traits
jpeg("plots/MultiPlot/meta/Sl_LesionSize_MAF20_meta.ManhattanPlot.jpg", width=8, height=4, units='in', res=600)
  ggplot(plotdata.wide, aes(x=Index, y=abs(Trait.Num)))+
    theme_bw()+
#    colScale+ #add this for BW version only
    geom_point(aes(color = factor(Phenos)))+
    labs(list(y="Number of Phenotypes", x="Chromosome position", title="Number of Phenotypes with Sig Fx"))+
    guides(col = guide_legend(nrow = 8, title="Traits"))+
   expand_limits(y=0)
dev.off()

#without the loop: just domest
jpeg("plots/MultiPlot/meta/Sl_LesionSize_MAF20_Domestmeta.ManhattanPlot.jpg", width=8, height=4, units='in', res=600)
ggplot(plotdata.wide)+
  theme_bw()+
  #    colScale+ #add this for BW version only
  geom_point(aes(x=Index, y=(Effect.Domesticated), color = "Domesticated"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Effect.Wild), color = "Wild"), alpha=1/2)+
  geom_point(aes(x=Index, y=(Effect.DmWoD), color = "(D-W)/D"), alpha=1/2)+
  labs(list(y="Effect size", x="Chromosome position", title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title="Traits"))+
  expand_limits(y=0)
dev.off()

#keep a subset of rows only if they have Domestication phenotypes
#Phenos column =/= IndPlants
d<-d[!(d$A=="B" & d$E==0),]
plotdata.wide.domest <- plotdata.wide[!(plotdata.wide$Phenos=="IndPlants"),]
write.csv(plotdata.wide.domest, "data/GWAS_files/04_bigRRoutput/SNP_overlap/TopSNPs_domestONLY_wide.csv")
tapply(HEM.plotdata$Trait, HEM.plotdata$Trait, length)

install.packages("eulerr")
library(eulerr)
jpeg("plots/MultiPlot/meta/VennDia_domest.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(Do=154, Wi=223, Di=378, "Do&Wi" = 76, 
                "Di&Do" = 4, "Di&Wi" = 2, 
                "Do&Wi&Di" = 2))
plot(fit, fill_opacity=0.3)
dev.off()
#can also use venn() in gplots package
#by gene
jpeg("plots/MultiPlot/meta/VennDia_domest_gene.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c(Do=74, Wi=93, Di=129, "Do&Wi" = 56, 
                "Di&Do" = 40, "Di&Wi" = 45, 
                "Do&Wi&Di" = 30))
plot(fit, fill_opacity=0.3)
dev.off()

#summarize gene annotations
GeneDat <- read.csv("data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest_geneannot_forR.csv")
names(GeneDat)
GeneDat2 <- reshape(GeneDat, idvar = c("Chrom","Segment","Pos","Gene","TotalTraits","Annot1","Annot2","Class"), timevar = "Trait", direction = "wide")
write.csv(GeneDat2, "data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest_genesummarized.csv")

GeneDat3 <- GeneDat2[,c("Chrom","Segment","Pos","Gene","Annot1","Annot2","Class")]
GeneDat3 <- reshape(GeneDat3, idvar = c("Chrom","Segment","Gene","Annot1","Annot2","Class"), timevar = "Pos", direction = "wide")
write.csv(GeneDat3, "data/GWAS_files/04_bigRRoutput/domestication/TopSNPs_domest_justGenes.csv")
