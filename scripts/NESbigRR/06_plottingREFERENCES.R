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