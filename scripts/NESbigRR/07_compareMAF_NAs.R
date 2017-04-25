#Nicole E Soltis
#041717
#compare bigRR SNP fx estimates from trueMAF20 with cutoffs of 10% NA sites (imputed) and 20% NA sites (imputed)
#------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcSolGWAS/")
library(plyr); library(ggplot2); library(eulerr)

HEM.plotdata.20NA <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_domest_trueMAF20_20NA.HEM.PlotFormat.csv") 
HEM.thresh.20NA <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_20NA/SlBc_domest_trueMAF20_20NA.HEM.Thresh.csv")
HEM.plotdata.20NA <- HEM.plotdata.20NA[,-c(1)]
HEM.thresh.20NA <- HEM.thresh.20NA[,-c(1)]

HEM.plotdata.10NA <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.PlotFormat.csv") 
HEM.thresh.10NA <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.Thresh.csv")
HEM.plotdata.10NA <- HEM.plotdata.10NA[,-c(1)]
HEM.thresh.10NA <- HEM.thresh.10NA[,-c(1)]

HEM.plotdata.10NA$Chrom.Seg.Pos <- paste(HEM.plotdata.10NA$Chrom, HEM.plotdata.10NA$Segment, HEM.plotdata.10NA$Pos, sep=".")
HEM.plotdata.20NA$Chrom.Seg.Pos <- paste(HEM.plotdata.20NA$Chrom, HEM.plotdata.20NA$Segment, HEM.plotdata.20NA$Pos, sep=".")

#---------------------------------------------------------
#Part 1:
#venn diagrams of genes > threshold for each
#take the top 50 over the threshold for each phenotype
TH95pos.10 <- HEM.thresh.10NA[1,]
for (i in 2:ncol(TH95pos.10)){
  assign(paste("TH95pos10_", names(TH95pos.10[i]), sep=""),as.numeric(TH95pos.10[i]))
}
TH95neg.10 <- HEM.thresh.10NA[5,]
for (i in 2:ncol(TH95neg.10)){
  assign(paste("TH95neg10_", names(TH95neg.10[i]), sep=""),as.numeric(TH95neg.10[i]))
}
TH99pos.10 <- HEM.thresh.10NA[3,]
for (i in 2:ncol(TH99pos.10)){
  assign(paste("TH99pos10_", names(TH99pos.10[i]), sep=""),as.numeric(TH99pos.10[i]))
}
TH99neg.10 <- HEM.thresh.10NA[7,]
for (i in 2:ncol(TH99neg.10)){
  assign(paste("TH99neg10_", names(TH99neg.10[i]), sep=""),as.numeric(TH99neg.10[i]))
}
TH999pos.10 <- HEM.thresh.10NA[4,]
for (i in 2:ncol(TH999pos.10)){
  assign(paste("TH999pos10_", names(TH999pos.10[i]), sep=""),as.numeric(TH999pos.10[i]))
}
TH999neg.10 <- HEM.thresh.10NA[8,]
for (i in 2:ncol(TH999neg.10)){
  assign(paste("TH999neg10_", names(TH999neg.10[i]), sep=""),as.numeric(TH999neg.10[i]))
}
#for 20
TH95pos.20 <- HEM.thresh.20NA[1,]
for (i in 2:ncol(TH95pos.20)){
  assign(paste("TH95pos20_", names(TH95pos.20[i]), sep=""),as.numeric(TH95pos.20[i]))
}
TH95neg.20 <- HEM.thresh.20NA[5,]
for (i in 2:ncol(TH95neg.20)){
  assign(paste("TH95neg20_", names(TH95neg.20[i]), sep=""),as.numeric(TH95neg.20[i]))
}
TH99pos.20 <- HEM.thresh.20NA[3,]
for (i in 2:ncol(TH99pos.20)){
  assign(paste("TH99pos20_", names(TH99pos.20[i]), sep=""),as.numeric(TH99pos.20[i]))
}
TH99neg.20 <- HEM.thresh.20NA[7,]
for (i in 2:ncol(TH99neg.20)){
  assign(paste("TH99neg20_", names(TH99neg.20[i]), sep=""),as.numeric(TH99neg.20[i]))
}
TH999pos.20 <- HEM.thresh.20NA[4,]
for (i in 2:ncol(TH999pos.20)){
  assign(paste("TH999pos20_", names(TH999pos.20[i]), sep=""),as.numeric(TH999pos.20[i]))
}
TH999neg.20 <- HEM.thresh.20NA[8,]
for (i in 2:ncol(TH999neg.20)){
  assign(paste("TH999neg20_", names(TH999neg.20[i]), sep=""),as.numeric(TH999neg.20[i]))
}

#All groups (4:6)
#keep only: SNPs over 99% Threshold
#now very few over 99.9% Thr
for (i in 4:6){
  assign(paste("HEMpos10.", names(HEM.plotdata.10NA[i]), sep=""), subset(HEM.plotdata.10NA, HEM.plotdata.10NA[i] > get(paste("TH99pos10_", names(HEM.plotdata.10NA[i]), sep="")), select=c(Chrom,Segment,Pos,Chrom.Seg.Pos,Index,i)))
  assign(paste("HEMneg10.", names(HEM.plotdata.10NA[i]), sep=""), subset(HEM.plotdata.10NA, HEM.plotdata.10NA[i] < get(paste("TH99neg10_", names(HEM.plotdata.10NA[i]), sep="")), select=c(Chrom,Segment,Pos,Chrom.Seg.Pos,Index,i)))
}

#20
for (i in 4:6){
  assign(paste("HEMpos20.", names(HEM.plotdata.20NA[i]), sep=""), subset(HEM.plotdata.20NA, HEM.plotdata.20NA[i] > get(paste("TH99pos20_", names(HEM.plotdata.20NA[i]), sep="")), select=c(Chrom,Segment,Pos,Chrom.Seg.Pos,Index,i)))
  assign(paste("HEMneg20.", names(HEM.plotdata.20NA[i]), sep=""), subset(HEM.plotdata.20NA, HEM.plotdata.20NA[i] < get(paste("TH99neg20_", names(HEM.plotdata.20NA[i]), sep="")), select=c(Chrom,Segment,Pos,Chrom.Seg.Pos,Index,i)))
}

#combine pos and neg by group
HEM.DmWoD.10 <- rbind(HEMpos10.DmWoD, HEMneg10.DmWoD)
HEM.Domesticated.10 <- rbind(HEMpos10.Domesticated, HEMneg10.Domesticated)
HEM.Wild.10 <- rbind(HEMpos10.Wild, HEMneg10.Wild)
HEM.DmWoD.20 <- rbind(HEMpos20.DmWoD, HEMneg20.DmWoD)
HEM.Domesticated.20 <- rbind(HEMpos20.Domesticated, HEMneg20.Domesticated)
HEM.Wild.20 <- rbind(HEMpos20.Wild, HEMneg20.Wild)

#then combine
HEM.Domesticated.10 <- rename(HEM.Domesticated.10, c("Domesticated" ="Effect"))
HEM.Domesticated.10$Trait <- "Domesticated"
HEM.Wild.10 <- rename(HEM.Wild.10, c("Wild" ="Effect"))
HEM.Wild.10$Trait <- "Wild"
HEM.DmWoD.10 <- rename(HEM.DmWoD.10, c("DmWoD" ="Effect"))
HEM.DmWoD.10$Trait <- "DmWoD"
HEM.Domesticated.20 <- rename(HEM.Domesticated.20, c("Domesticated" ="Effect"))
HEM.Domesticated.20$Trait <- "Domesticated"
HEM.Wild.20 <- rename(HEM.Wild.20, c("Wild" ="Effect"))
HEM.Wild.20$Trait <- "Wild"
HEM.DmWoD.20 <- rename(HEM.DmWoD.20, c("DmWoD" ="Effect"))
HEM.DmWoD.20$Trait <- "DmWoD"

HEM.topSNPs.10 <- rbind(HEM.Domesticated.10, HEM.Wild.10, HEM.DmWoD.10)
HEM.topSNPs.20 <- rbind(HEM.Domesticated.20, HEM.Wild.20, HEM.DmWoD.20)

#check for overlaps within sets
HEM.Domesticated.allNA <- merge(HEM.Domesticated.10, HEM.Domesticated.20, by="Chrom.Seg.Pos")
HEM.DmWoD.allNA <- merge(HEM.DmWoD.10, HEM.DmWoD.20, by="Chrom.Seg.Pos")
HEM.Wild.allNA <- merge(HEM.Wild.10, HEM.Wild.20, by="Chrom.Seg.Pos")

#venn diagram sample script:
#Domesticated
jpeg("paper/plots/ActualPaper/Supp/VennDia_99thr_Domest.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c("NA10"=(25144), "NA20"=(60526), "NA10&NA20"=(25109)))
plot(fit, fill_opacity=0)
dev.off()

#Wild
jpeg("paper/plots/ActualPaper/Supp/VennDia_99thr_Wild.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c("NA10"=(13132), "NA20"=(38157), "NA10&NA20"=(13117)))
plot(fit, fill_opacity=0)
dev.off()

#DmWoD
jpeg("paper/plots/ActualPaper/Supp/VennDia_99thr_DmWoD.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c("NA10"=(199), "NA20"=(274), "NA10&NA20"=(127)))
plot(fit, fill_opacity=0)
dev.off()

#top 200 snps/ pheno
top200.Wild.10 <- head(arrange(HEM.Wild.10, desc(abs(Effect))), n=200)
top200.Domesticated.10 <- head(arrange(HEM.Domesticated.10, desc(abs(Effect))), n=200)
top200.Wild.20 <- head(arrange(HEM.Wild.20, desc(abs(Effect))), n=200)
top200.Domesticated.20 <- head(arrange(HEM.Domesticated.20, desc(abs(Effect))), n=200)
top200.DmWoD.20 <- head(arrange(HEM.DmWoD.20, desc(abs(Effect))), n=200)

HEM.Domesticated.top200 <- merge(top200.Domesticated.10, top200.Domesticated.20, by="Chrom.Seg.Pos")
HEM.DmWoD.top200 <- merge(HEM.DmWoD.10, top200.DmWoD.20, by="Chrom.Seg.Pos")
HEM.Wild.top200 <- merge(top200.Wild.10, top200.Wild.20, by="Chrom.Seg.Pos")

#venn diagram sample script:
#Domesticated
jpeg("paper/plots/ActualPaper/Supp/VennDia_top200_Domest.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c("NA10"=(200), "NA20"=(200), "NA10&NA20"=(175)))
plot(fit, fill_opacity=0)
dev.off()

#Wild
jpeg("paper/plots/ActualPaper/Supp/VennDia_top200_Wild.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c("NA10"=(200), "NA20"=(200), "NA10&NA20"=(191)))
plot(fit, fill_opacity=0)
dev.off()

#DmWoD
jpeg("paper/plots/ActualPaper/Supp/VennDia_top200_DmWoD.jpg", width=4, height=4, units='in', res=600)
fit <- eulerr(c("NA10"=(199), "NA20"=(200), "NA10&NA20"=(102)))
plot(fit, fill_opacity=0)
dev.off()
#dev.off()
#----------------------------------------------------------------
#part 2: scatterplots to compare effect sizes genomewide
#merge HEM.plotdata.20NA with HEM.plotdata.10NA based on Index
HEM.10NA <- HEM.plotdata.10NA[,c("Chrom", "Segment", "Pos", "Domesticated", "Wild", "DmWoD", "Index")]
HEM.20NA <- HEM.plotdata.20NA[,c("Chrom", "Segment", "Pos", "Domesticated", "Wild", "DmWoD", "Index")]

HEM.10NA <- rename(HEM.10NA, c("Chrom"="Chrom.10", "Segment"="Segment.10", "Pos"="Pos.10", "Domesticated"="Domesticated.10", "Wild"="Wild.10", "DmWoD"="DmWoD.10"))
HEM.20NA <- rename(HEM.20NA, c("Chrom"="Chrom.20", "Segment"="Segment.20", "Pos"="Pos.20", "Domesticated"="Domesticated.20", "Wild"="Wild.20", "DmWoD"="DmWoD.20"))

#merge on Chrom.Seg.Pos
HEM.10NA$Chrom.Seg.Pos <- paste(HEM.10NA$Chrom.10, HEM.10NA$Segment.10, HEM.10NA$Pos.10, sep=".")
HEM.20NA$Chrom.Seg.Pos <- paste(HEM.20NA$Chrom.20, HEM.20NA$Segment.20, HEM.20NA$Pos.20, sep=".")
HEM.allNA <- merge(HEM.10NA, HEM.20NA, by="Chrom.Seg.Pos")

jpeg("paper/plots/ActualPaper/Supp/CompareNAthr_Domest.jpg", width=8, height=5, units='in', res=600)
ggplot(HEM.allNA, aes(x=Domesticated.10, y=Domesticated.20))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom.10)), alpha=0.2)+
  labs(list(y="Effect size estimate (20% NA)", title="Lesion Size on Domesticated Tomato", x="Effect size estimate (10% NA)"))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))
dev.off()
  
jpeg("paper/plots/ActualPaper/Supp/CompareNAthr_Wild.jpg", width=8, height=5, units='in', res=600)
ggplot(HEM.allNA, aes(x=Wild.10, y=Wild.20))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom.10)), alpha=0.2)+
  labs(list(y="Effect size estimate (20% NA)", title="Lesion Size on Wild Tomato", x="Effect size estimate (10% NA)"))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))
dev.off()

jpeg("paper/plots/ActualPaper/Supp/CompareNAthr_DmWoD.jpg", width=8, height=5, units='in', res=600)
ggplot(HEM.allNA, aes(x=DmWoD.10, y=DmWoD.20))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom.10)), alpha=0.2)+
  labs(list(y="Effect size estimate (20% NA)", title="Domestication Sensitivity of Lesion Size on Tomato", x="Effect size estimate (10% NA)"))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))
dev.off()

#basic version
plot(HEM.allNA$DmWoD.10, HEM.allNA$DmWoD.20)
mylm<-lm(HEM.allNA$DmWoD.20 ~ HEM.allNA$DmWoD.10)
abline(mylm, lty=2)
a <- predict(mylm, interval="confidence")
#lines(HEM.allNA$DmWoD.10, a[,2], lty=2)
#lines(HEM.allNA$DmWoD.10, a[,3], lty=2)
