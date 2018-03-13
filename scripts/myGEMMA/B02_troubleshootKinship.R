#Nicole E Soltis
#03/12/18

#troubleshoot kinship matrix
#www.r-graph-gallery.com for ideas
#------------------------------------------------------------------------
rm(list=ls())

setwd("~/Documents/GitRepos/BcSolGWAS/")
mykin <- read.table("data/GEMMA_files/04_GEMMAoutput/kmatrices/binMAF20NA10_kmatrix1.cXX.txt", header=FALSE)
mykin.num <- as.matrix(mykin)

library(RColorBrewer); library(gplots)
my_palette <- colorRampPalette(c("red", "orange","yellow", "green", "blue"))(n = 299)

# creates a 5 x 5 inch image

png("plots/resub/kmatrix1_v1_heatmap.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
heatmap.2(mykin.num, scale="row", 
        col=my_palette,
        main = "k-matrix Domestication",
        trace="none",    
        dendrogram="row")
dev.off()
