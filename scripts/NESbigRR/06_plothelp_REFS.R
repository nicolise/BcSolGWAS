#Nicole E Soltis
#011917
#Cut and paste code for some manhattan plots etc.

#--------------------------------------------------------------------

#greyscale Manhattan plot for Botrytis (16 Chr)
#create a custom color scale
library(RColorBrewer)
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80")
names(myColors) <- levels(HEM.plotdata$Chrom)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#when plotting, add +colScale + as a line

#colorblind-friendly palette in R
#black, yellow orange, lt blue, green, yellow, blue, red orange, pink
#for subsets: avoid 2 oranges, yellow with yellow orange, green with blue
+scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

#reshape data to and from long // wide
Top50SNP.wide <- reshape(Top50SNP, 
                         timevar = "Plant",
                         idvar = c("Chrom","Segment","Pos","Index"),
                         direction = "wide")

Top50SNPs.long.DM <- reshape(HEM.plotdata, 
                             varying = c("Domesticated", "Wild", "DmWoD"),
                             v.names = "Effect",
                             timevar = "Trait",
                             times = c("Domesticated", "Wild", "DmWoD"),
                             direction = "long")

#original INDEXING without multiple segments per chromosome
#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

#want to figure out where to add +500 to draw breaks between chromosomes
#Redo the positions to make them sequential		-- accurate position indexing
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