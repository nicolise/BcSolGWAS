
#########################
# This makes the bigRR_update run through the GPU
	# You need to do this first to mask the native 'bigRR_update' in the bigRR package

bigRR_update <- function (obj, Z, family = gaussian(link = identity), tol.err = 1e-06, 
    tol.conv = 1e-08) 
{
    w <- as.numeric(obj$u^2/(1 - obj$leverage))
    w[w < tol.err] <- tol.err
    bigRR(y = obj$y, X = obj$X, Z = Z, family = family, weight = w, 
        tol.err = tol.err, tol.conv = tol.conv, GPU = TRUE )
}
########################

library(bigRR)

#Get data
SNPs <- read.csv("../3_Emma/JulinsAcc3_MAF20.csv", row.names = 1)
SNPs <- as.matrix(t(SNPs))
for(i in 1:dim(SNPs)[1]) {
	SNPs[i,] <- as.numeric(SNPs[i,])
}

Phenos <- read.csv("../3_Emma/Phenos.LSMeans.csv", row.names = 1)
dat1 <- read.csv("../2_LSMeans2/Cam.LSMeans.csv",row.names=1)
dat1 <- dat1[-73,]
dat2 <- read.csv("../2_LSMeans2/Les.LSmeans.csv",row.names=1)
dat <- cbind(dat1,(dat2)[,21:25])

Col0.Phenos <- dat[13,] 
dat <- dat[-13,]

outpt.BLUP <- colnames(SNPs)
outpt.HEM <- colnames(SNPs)
thresh.BLUP <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)
thresh.HEM <- list("0.95Thresh" = NA, "0.975Thresh" = NA, "0.99Thresh" = NA, "0.999Thresh" = NA)

#Calculate BLUP and HEMs for all phenotypes
for(i in 1:dim(dat)[2]) { 
	print(colnames(dat)[i])
	MyX <- matrix(1, dim(dat)[1], 1)

	Pheno.BLUP.result <- bigRR(y = dat[,i], X = MyX, Z = SNPs, GPU = TRUE)
	Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)
	
	
	outpt.BLUP <- cbind(outpt.BLUP, Pheno.BLUP.result$u)
	outpt.HEM <- cbind(outpt.HEM, Pheno.HEM.result$u)
	
	#Permute Thresholds for Phenos
	perm.u.BLUP <- vector()
	perm.u.HEM <- vector()
	for(p in 1:1000) {
		if(p %% 10 == 0) {print(paste("Thresh sample:", p, "--", Sys.time()))}
		temp.Pheno <- sample(dat[,i], length(dat[,i]), replace = FALSE)
		try(temp.BLUP  <- bigRR(y = temp.Pheno, X = MyX, Z = SNPs, GPU = TRUE),silent = TRUE)
		temp.HEM <- bigRR_update(Pheno.BLUP.result, SNPs)
		perm.u.BLUP <- c(perm.u.HEM, temp.BLUP$u)
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

#Write results to output
write.csv(rbind(thresh.BLUP$"0.95Thresh",thresh.BLUP$"0.975Thresh",thresh.BLUP$"0.99Thresh",thresh.BLUP$"0.999Thresh",outpt.BLUP),"CamLes.BLUP.csv")
write.csv(rbind(thresh.HEM$"0.95Thresh",thresh.HEM$"0.975Thresh",thresh.HEM$"0.99Thresh",thresh.HEM$"0.999Thresh",outpt.HEM),"CamLes.HEM.csv")

#Write just the positive positions
sig.HEM <- data.frame()
for(i in 1:dim(outpt.HEM)[1]) {
	if(i %% 1000 == 0) {print(paste(i, "--", Sys.time()))}
	if(any(abs(as.numeric(outpt.HEM[i,2:11]))-abs(thresh.HEM$"0.99Thresh")>0)) {
		sig.HEM <- unname(as.matrix(rbind(sig.HEM,outpt.HEM[i,])))
	}
}
colnames(sig.HEM) <- colnames(outpt.HEM)
write.csv(rbind(thresh.HEM$"0.99Thresh",sig.HEM),"CamLes.HEM.99Sig.csv")



Phenos[13,] -> Col0.LSMeans
Phenos[-13,] -> Phenos

MyX <- matrix(1, dim(Phenos)[1], 1)

Pheno.BLUP.result <- bigRR(y = Phenos[,1], X = MyX, Z = SNPs)

Pheno.HEM.result <- bigRR_update(Pheno.BLUP.result, SNPs)






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

#Import data
HEM.plotdata <- read.csv("CamLes.HEM.PlottingFormat.csv")

#Extract thresholds
HEM.Perm.Threshold <- as.vector(HEM.plotdata[1,])
HEM.plotdata <- HEM.plotdata[-1,]

#Reformat Chromosomes and Positions
HEM.plotdata$Chrom <- gsub("IV", 4, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("V", 5, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("III", 3, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("II", 2, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- gsub("I", 1, HEM.plotdata$Chrom)
HEM.plotdata$Chrom <- as.numeric(as.character(HEM.plotdata$Chrom))
HEM.plotdata$Pos <- as.numeric(as.character(HEM.plotdata$Pos))


#Give column names to the thresholds from the HEM list
for(j in 1:length(thresh.HEM)) {
	names(thresh.HEM[[j]]) <- colnames(dat)
}


#Make plotting variables
HEM.plotdata$Index=NA
ticks=NULL
lastbase=0
		
#Redo the positions to make them sequencial		
for (i in unique(HEM.plotdata$Chrom)) {
	print(i)
	if (i=="Chromosome1") {
		HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos
	}	else {
		previous = lastbase=lastbase+tail(subset(HEM.plotdata,Chrom==i-1)$Pos, 1)
		lastbase=lastbase+tail(subset(HEM.plotdata,Chrom==previous)$Pos, 1)
		HEM.plotdata[HEM.plotdata$Chrom==i, ]$pos=HEM.plotdata[HEM.plotdata$Chrom==i, ]$Pos+lastbase
	}
	ticks=c(ticks, HEM.plotdata[HEM.plotdata$Chrom==i, ]$index[floor(length(HEM.plotdata[HEM.plotdata$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(HEM.plotdata$Index),max(HEM.plotdata$Index))

#set my colors
mycols=rep(c("gray10","gray60"),max(HEM.plotdata$Chrom))

#Make plot (single)
plot=qplot(pos,abs(Camalexin.ng.LesPer._Apple517.HEM),data=HEM.plotdata, ylab="SNP Effect Estimate" , colour=factor(Chrom))
plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(HEM.plotdata$Chrom)))
plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
plot=plot+scale_colour_manual(values=mycols)
plot=plot + geom_hline(yintercept=HEM.Perm.Threshold$Camalexin.ng.LesPer._Apple517.HEM, colour="red")



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








