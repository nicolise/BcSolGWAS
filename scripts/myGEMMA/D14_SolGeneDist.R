#Nicole E Soltis
#05/07/18

#-----------------------------------------------------------------------------------------------
#calculate pairwise distances in R between tomato accessions based on SNP data
#https://popgen.nescent.org/2015-05-18-Dist-SNP.html
#non-evolutionary genetic distances
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
mydata <- read.csv("SlGenome/Sanger/solcap_tomato_opa_genotype_NES.csv")

library("poppr")
library("pegas")

#create a matrix with only genotype data
mygenos <- mydata[,c(4:51)]
ind <- as.character(mydata$Variety) # individual labels 
population <- rep(1,96) # dummy population labels-- can add this into myGenind or not
myGenind <- df2genind(mygenos, ploidy = 1, ind.names = ind, sep="")
myGenind2 <- genind2loci(myGenind)

#Euclidean distance in {adegenet}
distgenEUCL <- dist(myGenind, method = "euclidean", diag = FALSE, upper = FALSE, p = 1)
hist(distgenEUCL)
blah <- data.matrix(distgenEUCL)
setwd("~/Projects/BcSolGWAS")
jpeg("paper/plots/addGEMMA/SolDistance.jpg", width=10, height=10, units='in', res=600)
heatmap(blah)
dev.off()

#-----------------------------------------------------------------------------------------------------------
#here's another option : Sim 2012 
rm(list=ls())
setwd("~/Projects/BcSolGWAS/data")
mysim <- read.csv("SlGenome/Sim2012/TableS2_NES.csv")

library("poppr")
library("pegas")

#remove row 2 (key for my tomato dataset)
mygenotypes_sim <- mysim[c(1),]
mysim <- mysim[-c(1),]

#create a matrix with only genotype data
mysim$Position <- as.numeric(mysim$Position..Mbp.) * 1000
mysim$Locus <- paste(mysim$Chromosome, mysim$Position, sep=".")
mygenosim <- mysim[,c(433,6:431)]
Bcgenos.domest <- c("LA0410","LA2706","LA3008","LA3475","LA4345","LA4355","Gold.Nugget","Heinz.1706","M82","Moneymaker","San.Marzano")
Bcgenos.wild <- c("LA0480","LA1547","LA1589","LA1684","LA2093","LA2176")



#transpose
myloci <- mygenosim[,1]
mygenosim <- mygenosim[,2:427]
mygenosim_t <- as.data.frame(t(mygenosim))
names(mygenosim_t) <- myloci

ind <- as.character(names(mygenosim)) # individual labels 

myGenind <- df2genind(mygenosim_t, ploidy = 2, ind.names = ind, sep="")
myGenind2 <- genind2loci(myGenind)

#Euclidean distance in {adegenet}
distgenEUCL <- dist(myGenind, method = "euclidean", diag = FALSE, upper = FALSE, p = 1)
hist(distgenEUCL)
blah <- data.matrix(distgenEUCL)

library(gplots)
#color code names
mycolor.samples <- as.data.frame(ind)
mycolor.samples$colors <- ifelse(mycolor.samples$ind %in% Bcgenos.domest, "Blue", ifelse(mycolor.samples$ind %in% Bcgenos.wild, "Orange","Black"))
sample.types <- mycolor.samples$colors
setwd("~/Projects/BcSolGWAS")
jpeg("paper/plots/addGEMMA/SolDistance_Sim2012.jpg", width=30, height=30, units='in', res=600)
heatmap.2(blah, trace="none", colCol=sample.types)
dev.off()

#or just the dendrogram
setwd("~/Projects/BcSolGWAS")
jpeg("paper/plots/addGEMMA/SolDistance_dendro_Sim2012.jpg", width=40, height=10, units='in', res=600)
plot(hclust(distgenEUCL), cex=0.5)
dev.off()
