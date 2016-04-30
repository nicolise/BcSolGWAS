#Nicole E Soltis 
#from script by Jason A Corwin, Modified by Rachel Fordyce
#to run bigRR on Linux GPU for GWAS
#---------------------------------------------------------------

rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")

#########################
# This makes the bigRR_update run through the GPU
# You need to do this first to mask the native 'bigRR_update' in the bigRR package
# one alternative to family = gaussian(link = identity) is family = poisson(link = log)
bigRR_update <- function (obj, Z, family = poisson(link = log), tol.err = 1e-06, 
                          tol.conv = 1e-08) 
{
  w <- as.numeric(obj$u^2/(1 - obj$leverage))
  w[w < tol.err] <- tol.err
  bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE )
}
########################
#NOTE1 FROM RACHEL:  we need bigRR1.3-9 to get GPU option
# must download R-Forge version bigRR1.3-9tar.gz and manually install
# https://r-forge.r-project.org/R/?group_id=1301
# install.packages("bigRR", repos="http://R-Forge.R-project.org")

#NOTE2 FROM RACHEL: need package 'gputools' but CRAN version fails to install
# must first install Nvidia's CUDA toolkit -current version is 7.5
# installed from developer.nvidia.com/cuda-downloads

library(bigRR) #check if version is 1.3-9

#Get genotype data
SNPs <- read.csv("binSNP_bigRR_MAF20hp.csv", row.names = 1)
FullSNPs <- SNPs
SNPs <- FullSNPs
#add a column with position as chr.base
SNPs$Chr.Base <- do.call(paste, c(SNPs[c("X.CHROM","POS")], sep="."))

#Rachel's attempt------------------=====================================

rownames(SNPs) <- SNPs[,95] #set the new column of chrom.base as rownames - this could maybe be written as: rownames(SNPs) <- SNPs$Chr.Base?
SNPs <- SNPs[,4:94] #take out first three cols (X.CHROM, POS, REF) and new last col (Chr.Base). dim(SNPs) should now be [345485, 91], colnames(SNPs) are all Bc Isolates, rownames(SNPs) are all Chr.Base

#Rachel's attempt ----FIN-----============================================

#--Rachel commented out--SNPs <- SNPs[,-1]
#make SNPs numeric
# myconvert <- colnames(SNPs[1:91])
# for (i in myconvert)
#   {
#   SNPs[[i]] <- as.numeric(SNPs[[i]])
# }


#what does this do? 
#loads SNP dataframe as a matrix?
#makes SNP states numeric (also transposes SNP matrix)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
  SNPs[i,] <- as.numeric(SNPs[i,])
}


Phenos <- read.csv("Sl_Pheno_bigRR.csv", row.names = 1)
dat <- as.data.frame((Phenos[,2:13]))  #INSERT PHENOTYPE COLUMNS HERE
#e.g. LesionGreen as.data.frame(c(Phenos[,31:32],Phenos[,34:35]))

#should I remove reference (B05.10 I assume) phenotypes and genotypes from list?
# B05.10.Phenos <- dat[64,]
# dat <- dat[-64,]


outpt.BLUP <- colnames(SNPs)
outpt.HEM <- colnames(SNPs)
thresh.BLUP <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)
thresh.HEM <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)

#Calculate BLUP and HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { #i will be each isolate
  print(colnames(dat)[i])
  MyX <- matrix(1, dim(dat)[1], 1)
  
  Pheno.BLUP.result <- bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE)
  Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)
  
  
  outpt.BLUP <- cbind(outpt.BLUP, Pheno.BLUP.result$u)
  outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
  
  #Permute Thresholds for Phenos - this is what takes forever
  perm.u.BLUP <- vector()
  perm.u.HEM <- vector()
  for(p in 1:1000) {  
    if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
    temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE)
    try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE),silent = TRUE)
    temp.HEM <- bigRR_update(temp.BLUP, SNPs) #had to change this from Jason's script - was bigRR_update(Pheno.BLUP.result...
    perm.u.BLUP <- c(perm.u.BLUP, temp.BLUP$u) #had to change from ...c(perm.u.HEM...)
    perm.u.HEM <- c(perm.u.HEM, temp.HEM$u)
    
  }
  #write.csv(perm.u.HEM, paste("PermEffects_",colnames(dat)[i],".csv",sep=""))
  thresh.BLUP$"0.95Thresh"[i] <- quantile(perm.u.BLUP,0.95)
  thresh.BLUP$"0.975Thresh"[i] <- quantile(perm.u.BLUP,0.975)
  thresh.BLUP$"0.99Thresh"[i] <- quantile(perm.u.BLUP,0.99)
  thresh.BLUP$"0.999Thresh"[i] <- quantile(perm.u.BLUP,0.999)
  thresh.HEM$"0.95Thresh"[i] <- quantile(perm.u.HEM,0.95)
  thresh.HEM$"0.975Thresh"[i] <- quantile(perm.u.HEM,0.975)
  thresh.HEM$"0.99Thresh"[i] <- quantile(perm.u.HEM,0.99)
  thresh.HEM$"0.999Thresh"[i] <- quantile(perm.u.HEM,0.999)
  colnames(outpt.BLUP)[i+1] <- paste(colnames(dat)[i],"BLUP",sep=".")
  colnames(outpt.HEM)[i+1] <- paste(colnames(dat)[i],"HEM",sep=".")
}

#Give column names to the thresholds from the HEM list
for(j in 1:length(thresh.HEM)) {
  names(thresh.HEM[[j]]) <- colnames(dat)
}

#RF-give row names to thresh.HEM and thresh.BLUP so that threshhold values will line up correctly with phenotypes, and you can see which threshold value is displayed
thresh.BLUP$"0.95Thresh" <- c("0.95 Thresh", thresh.BLUP$"0.95Thresh")
thresh.BLUP$"0.975Thresh" <- c("0.975 Thresh", thresh.BLUP$"0.975Thresh")
thresh.BLUP$"0.99Thresh" <- c("0.99 Thresh", thresh.BLUP$"0.99Thresh")
thresh.BLUP$"0.999Thresh" <- c("0.999 Thresh", thresh.BLUP$"0.999Thresh")
thresh.HEM$"0.95Thresh" <- c("0.95 Thresh", thresh.HEM$"0.95Thresh")
thresh.HEM$"0.975Thresh" <- c("0.975 Thresh", thresh.HEM$"0.975Thresh")
thresh.HEM$"0.99Thresh" <- c("0.99 Thresh", thresh.HEM$"0.99Thresh")
thresh.HEM$"0.999Thresh" <- c("0.999 Thresh", thresh.HEM$"0.999Thresh")

#Write results to output
write.csv(rbind(thresh.BLUP$"0.95Thresh",thresh.BLUP$"0.975Thresh",thresh.BLUP$"0.99Thresh",thresh.BLUP$"0.999Thresh",outpt.BLUP),"SolanumLesionSizePoisson.BLUP.csv")
write.csv(rbind(thresh.HEM$"0.95Thresh",thresh.HEM$"0.975Thresh",thresh.HEM$"0.99Thresh",thresh.HEM$"0.999Thresh",outpt.HEM),"SolanumLesionSizePoisson.HEM.csv")

# #This part failed for me -- NES
# used to store upper percentile of SNPS
# #Write just the positive positions (RF- effect size greater than .99 thresh)
# sig.HEM <- data.frame()
# for(i in 1:dim(outpt.HEM)[1]) {
#   if(i %% 1000 == 0) {print(paste(i, "--", Sys.time()))}
#   if(any(abs(as.numeric(outpt.HEM[i,2:5]))-abs(as.numeric(thresh.HEM$"0.99Thresh"[2:5]))>0)) { 
#     # change output.HEM[1,2:5] to appropriate column dims
#     sig.HEM <- unname(as.matrix(rbind(sig.HEM,outpt.HEM[i,])))
#   }
# }
# colnames(sig.HEM) <- colnames(outpt.HEM)
# write.csv(rbind(thresh.HEM$"0.99Thresh",sig.HEM),"SolanumLesionSizePoisson.HEM.99Sig.csv")



# Phenos[13,] -> Col0.LSMeans
# Phenos[-13,] -> Phenos
# 
# MyX <- matrix(1, dim(Phenos)[1], 1)
# 
# Pheno.BLUP.result <- bigRR(y = Phenos[,1], X = MyX, Z = SNPs)
# 
# Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)




#RF the code below generates plots comparing .BLUP and .HEM methods - not necessary

split.screen(c(1, 2))
split.screen(c(2, 1), screen = 1)
screen(3); plot(abs(Pheno.BLUP.result$u), cex = .6, col = 'slateblue')
screen(4); plot(abs(Pheno.HEM.result$u), cex = .6, col = 'olivedrab')
screen(2); plot(abs(Pheno.BLUP.result$u), abs(Pheno.HEM.result$u), cex = .6, pch = 19, col = 'darkmagenta')

split.screen(c(1, 2))
split.screen(c(2, 1), screen = 1)
screen(3); plot(abs(Pheno.BLUP.result$leverage), cex = .6, col = 'slateblue')
screen(4); plot(abs(Pheno.HEM.result$leverage), cex = .6, col = 'olivedrab')

############################################################################


###Plotting the HEM results

#Load plotting package
library(ggplot2)
library(grid)

#Import data (reorganized from script ReformatBigRRouts.R)
HEM.plotdata <- read.csv("data/GWAS_files/04_bigRRoutput/Sl_LesionSize.HEM.PlotFormat.csv")

HEM.plotdata$Pos <- as.character(HEM.plotdata$Pos)#ensure that position data is not in scientific notation
HEM.plotdata <- HEM.plotdata[-c(1:2)]

#get threshhold values 
HEM.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/Sl_LesionSize.HEM.Thresh.csv")
HEM.thresh <- HEM.thresh[,-c(1:2)]
TH99 <- HEM.thresh[3,]
TH99_LA0410 <- as.numeric(TH99[2])
TH99_LA0480 <- as.numeric(TH99[3])
TH99_LA1547 <- as.numeric(TH99[4])
TH99_LA1589 <- as.numeric(TH99[5])
TH99_LA1684 <- as.numeric(TH99[6])
TH99_LA2093 <- as.numeric(TH99[7])
TH99_LA2176 <- as.numeric(TH99[8])
TH99_LA2706 <- as.numeric(TH99[9])
TH99_LA3008 <- as.numeric(TH99[10])
TH99_LA3475 <- as.numeric(TH99[11])
TH99_LA4345 <- as.numeric(TH99[12])
TH99_LA4355 <- as.numeric(TH99[13])

TH999 <- HEM.thresh[4,]
TH999_LA0410 <- as.numeric(TH999[2])
TH999_LA0480 <- as.numeric(TH999[3])
TH999_LA1547 <- as.numeric(TH999[4])
TH999_LA1589 <- as.numeric(TH999[5])
TH999_LA1684 <- as.numeric(TH999[6])
TH999_LA2093 <- as.numeric(TH999[7])
TH999_LA2176 <- as.numeric(TH999[8])
TH999_LA2706 <- as.numeric(TH999[9])
TH999_LA3008 <- as.numeric(TH999[10])
TH999_LA3475 <- as.numeric(TH999[11])
TH999_LA4345 <- as.numeric(TH999[12])
TH999_LA4355 <- as.numeric(TH999[13])

TH95 <- HEM.thresh[1,]
TH95_LA0410 <- as.numeric(TH95[2])
TH95_LA0480 <- as.numeric(TH95[3])
TH95_LA1547 <- as.numeric(TH95[4])
TH95_LA1589 <- as.numeric(TH95[5])
TH95_LA1684 <- as.numeric(TH95[6])
TH95_LA2093 <- as.numeric(TH95[7])
TH95_LA2176 <- as.numeric(TH95[8])
TH95_LA2706 <- as.numeric(TH95[9])
TH95_LA3008 <- as.numeric(TH95[10])
TH95_LA3475 <- as.numeric(TH95[11])
TH95_LA4345 <- as.numeric(TH95[12])
TH95_LA4355 <- as.numeric(TH95[13])


# #Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("Chromosome", "", HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))

#Make plotting variables
HEM.plotdata$Index = NA
ticks = NULL
lastbase = 0

##NES failing here
#Redo the positions to make them sequencial		-- accurate position indexing
##RFF code isn't working... replaced "Pos" with "pos"
for (i in unique(HEM.plotdata$Chrom)) {
  print(i)
  if (i==1) {
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
  }	else {
    lastbase=lastbase+tail(subset(HEM.plotdata,HEM.plotdata$Chrom==i-1)$Pos, 1)
    HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
  }
  ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#make plots for each phenotype
jpeg("plots/Sl_LesionSize_LA0410.ManhattanPlot.jpg")
qplot(Index,abs(LA0410), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA0410", colour=factor(Chrom)) +
 geom_hline(yintercept=TH99_LA0410) +
  geom_text(aes(0,TH99_LA0410, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black") +
geom_hline(yintercept=TH95_LA0410, colour = "blue") +
  geom_text(aes(0,TH95_LA0410, label = ".95 Threshold", vjust = 1.5, hjust = .05), col = "blue")
dev.off()

jpeg("plots/Sl_LesionSize_LA0410.ManhattanPlot_999.jpg")
qplot(Index,abs(LA0410), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA0410", colour=factor(Chrom)) +
  geom_hline(yintercept=TH999_LA0410) +
  geom_text(aes(0,TH999_LA0410, label = ".999 Threshold", vjust = 1.5, hjust = .05), col = "black")
dev.off()

#how many SNPs are above a certain threshhold?
sum(HEM.plotdata$LA0410 >= TH999_LA0410) #spits out number
#6491

#rest of plots
jpeg("plots/Sl_LesionSize_LA4355.ManhattanPlot_99.jpg")
qplot(Index,abs(LA4355), data=HEM.plotdata, ylab="SNP Effect Estimate" , 
      main = "LesionSize_LA4355", colour=factor(Chrom)) +
  geom_hline(yintercept=TH99_LA4355, colour = "black") +
  geom_text(aes(0,TH99_LA4355, label = ".99 Threshold", vjust = 1.5, hjust = .05), col = "black")
dev.off()

#not sure why chromosomes are overlapping
m <- ggplot(HEM.plotdata, aes(x=Index, fill=as.factor(Chrom)))
m + geom_bar() 

#indices ARE overlapping
range(subset(HEM.plotdata,HEM.plotdata$Chrom==1)$Index)
range(subset(HEM.plotdata,HEM.plotdata$Chrom==2)$Index)
range(subset(HEM.plotdata,HEM.plotdata$Chrom==12)$Index)
range(subset(HEM.plotdata,HEM.plotdata$Chrom==13)$Index)

#set my colors
mycols=rep(c("gray10","gray60"),max(HEM.plotdata$Chrom))

#Make plot (single)
# plot=qplot(Index,abs(Lesion.Size_Apple517.HEM),data=HEM.plotdata, ylab="SNP Effect Est" , colour=factor(Chrom))
# plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(HEM.plotdata$Chrom)))
# plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
# plot=plot+scale_colour_manual(values=mycols)
# plot=plot + geom_hline(yintercept=HEM.Perm.Threshold$Lesion.Size_Apple517.HEM, colour="red")


###############################
###Make a BiPlot for two phenotypes (Cam and Lesion)
#This code has several dependent objects from the single HEM plotting code above


for(k in c(3:7)) {
  #Make temp data for figure
  print(names(thresh.HEM[[1]])[k-2])
  
  CamLes.plot.data<-cbind("Cam",HEM.plotdata[,k],HEM.plotdata[,c(1,13)])
  colnames(CamLes.plot.data)[1:2] <- c("Pheno","Effect")
  CamLes.plot.data <- rbind(CamLes.plot.data,cbind(Pheno="Lesion",Effect=-HEM.plotdata[,k+5],HEM.plotdata[,c(1,13)]))
  CamLes.plot.data <- cbind(CamLes.plot.data,PhenoChrom=paste(CamLes.plot.data[,1],CamLes.plot.data[,3]))
  
  #Scale the effects and threshold
  temp.scaled <- scale(c(thresh.HEM[[3]][k-2],CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Cam"]))
  Cam.Scaled.Thresh <- abs(temp.scaled[1])
  CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Cam"] <- abs(temp.scaled[-1])
  
  temp.scaled <- scale(c(thresh.HEM[[3]][k+3],CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Lesion"]))
  Les.Scaled.Thresh <- -abs(temp.scaled[1])
  CamLes.plot.data$Effect[CamLes.plot.data$Pheno == "Lesion"] <- -abs(temp.scaled[-1])
  
  rm(temp.scaled)
  
  
  #Plot Data
  
  mycols2=c("royalblue","royalblue4","royalblue","royalblue4","royalblue","olivedrab3","olivedrab","olivedrab3","olivedrab","olivedrab3")
  
  
  if(grepl("Ctrl",colnames(HEM.plotdata)[k])){
    plot2=qplot(pos,Effect,data=CamLes.plot.data[CamLes.plot.data$Pheno=="Cam",], ylab="|Scaled SNP Effect Estimate|" , colour=factor(PhenoChrom))
    
    ##All the thresholds
    #thresholds <- c(thresh.HEM[[1]][k-2],thresh.HEM[[2]][k-2],thresh.HEM[[3]][k-2],thresh.HEM[[4]][k-2])
    #thresholds <- c(thresh.HEM[[3]][k-2])
    
    plot2=plot2 + geom_hline(yintercept=Cam.Scaled.Thresh, colour="red")
    
  } else {
    plot2=qplot(pos,Effect,data=CamLes.plot.data, ylab="|Scaled SNP Effect Estimate|" , colour=factor(PhenoChrom))
    ##All the thresholds
    #thresholds <- c(thresh.HEM[[1]][c(k-2,k+3)],thresh.HEM[[2]][c(k-2,k+3)],thresh.HEM[[3]][c(k-2,k+3)],thresh.HEM[[4]][c(k-2,k+3)])
    #thresholds <- c(thresh.HEM[[3]][c(k-2,k+3)])
    #thresholds[grepl("Lesion",names(thresholds))] <- thresholds[grepl("Lesion",names(thresholds))] * -1
    plot2=plot2 + geom_hline(yintercept=c(Cam.Scaled.Thresh,Les.Scaled.Thresh), colour="red")
    
    
  }
  #plot2=qplot(pos,Effect,data=CamLes.plot.data, ylab="SNP Effect Estimate" , colour=factor(PhenoChrom))
  plot2=plot2+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(HEM.plotdata$Chrom)))
  #plot2=plot2+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
  plot2=plot2+scale_colour_manual(values=mycols2) + theme(legend.position="none", axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"))
  #plot2=plot2 + geom_hline(yintercept=c(as.numeric(thresh.HEM[k-2]),-as.numeric(thresh.HEM[k+3])), colour="red")
  
  #write to file.  May want something in a vector graphic format
  jpeg(filename = paste(strsplit(names(thresh.HEM[[1]][k-2]),"_")[[1]][2],".jpeg",sep=""), width = 1500, height = 800)
  print(plot2)
  dev.off()
}



###############################
###Make a HEM plot of a region (Cam and Lesion)
#This code has several dependent objects from the single HEM plotting code above

#PRLM3 (Chr4 Pos 9560135:9565454) in this example.  Just change the Chrom and positions to find another gene.
Reg.min <- 9560135 - 50000
Reg.max <- Reg.min + 1000000


Reg.min <- 13750000 - 50000
Reg.max <- Reg.min + 100000


Region2Plot <- HEM.plotdata[HEM.plotdata$Chrom == 4 & HEM.plotdata$Pos %in% c(Reg.min:Reg.max),]
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

plot.reg <- plot.reg + geom_segment(aes(x = gene.start, y = -0.001, xend = gene.stop, yend = -.001), colour = "red")

#####
###For some reason, I had to use the following code at the terminal to get this to print right.
#xwd > PolyMorphicRegOnChr4.xwd
#convert PolyMorphicRegOnChr4.xwd PolyMorphicRegOnChr4.jpg

jpeg(filename = "PolyMorphicRegionOnChr4.jpg", width = 1500, height = 800)
print(plot2)
dev.off()