#Nicole E Soltis
#030818

#B05_QQ_GEMMAbigRR.R
#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

myGEMMA <- read.csv("data/GEMMA_files/04_analysis/GEMMA_allDWS.csv")
mybigRR <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.PlotFormat.csv")
mybigRR <- mybigRR[,-c(1)]
bigRR.thresh <- read.csv("data/GWAS_files/04_bigRRoutput/trueMAF_10NA/SlBc_domest_trueMAF20_10NA.HEM.Thresh.csv")
bigRR.thresh <- bigRR.thresh[,-c(1)]
myGEMMA <- myGEMMA[,-c(1)]

head(myGEMMA)

#for filtering bigRR SNPs:
TH99pos <- bigRR.thresh[3,]
for (i in 2:ncol(TH99pos)){
  assign(paste("TH99pos_", names(TH99pos[i]), sep=""),as.numeric(TH99pos[i]))
}
TH99neg <- bigRR.thresh[7,]
for (i in 2:ncol(TH99neg)){
  assign(paste("TH99neg_", names(TH99neg[i]), sep=""),as.numeric(TH99neg[i]))
}

mybigRR$TotTraits <- ifelse(abs(mybigRR$Domesticated) > TH99pos$Domesticated & abs(mybigRR$Wild) > TH99pos$Wild & abs(mybigRR$DmWoD) > TH99pos$DmWoD, "ALL", 
                                 ifelse(abs(mybigRR$Domesticated) > TH99pos$Domesticated & abs(mybigRR$Wild) > TH99pos$Wild, "DW",
                                        ifelse(abs(mybigRR$Wild) > TH99pos$Wild & abs(mybigRR$DmWoD) > TH99pos$DmWoD, "WS",
                                               ifelse(abs(mybigRR$Domesticated) > TH99pos$Domesticated & abs(mybigRR$DmWoD) > TH99pos$DmWoD, "DS",
                                                      ifelse(abs(mybigRR$Domesticated) > TH99pos$Domesticated, "D",
                                                             ifelse(abs(mybigRR$Wild) > TH99pos$Wild, "W", 
                                                                    ifelse(abs(mybigRR$DmWoD) > TH99pos$DmWoD, "S", "none")))))))

table(mybigRR$TotTraits)

#for comparisons: T4 and B05.10 have different alignments & different chromosome numbers
#CAN rank all SNPs by absolute value of effect size and plot them side by side.
#but mybigRR has an extra 30,000 SNPs- how to compare? 

#1. randomly remove SNPs from mybigRR until equal numbers
#2. rank SNPs positionally by index = posrank.G posrank.B
#3. pair SNPs  on posrank
#4. plot pscore (GEMMA) vs. effect estimate (bigRR) for each phenotype

#1. randomly remove SNPs from mybigRR until equal numbers
#going to only remove them from the "none" significant traits
#split df by TotTraits=none, else
#randomly remove 35577 rows from TotTraits=none
#combine back together
mybigRR.ss <- mybigRR[mybigRR$TotTraits!="none",]
mybigRR.na <- mybigRR[mybigRR$TotTraits=="none",]
mybigRR.na <- mybigRR.na[-sample(1:nrow(mybigRR.na), 35577),]
mybigRR.mt <- rbind(mybigRR.ss, mybigRR.na)

#2. rank SNPs positionally by index = posrank.G posrank.B
mybigRR.mt <- mybigRR.mt[order(mybigRR.mt$Index),]
mybigRR.mt$posrank.B <- 1:nrow(mybigRR.mt)
myGEMMA.mt <- myGEMMA[order(myGEMMA$Index),]
myGEMMA.mt$posrank.G <- 1:nrow(myGEMMA.mt)
#3. pair SNPs  on posrank
myQQcomp <- cbind(mybigRR.mt, myGEMMA.mt)
#double check that they're ordered appropriately
plot(myQQcomp$posrank.B, myQQcomp$posrank.G)
#4. plot pscore (GEMMA) vs. effect estimate (bigRR) for each phenotype
plot(abs(myQQcomp$Domesticated), -log(myQQcomp$pscore.D))
plot(abs(myQQcomp$Wild), -log(myQQcomp$pscore.W))
plot(abs(myQQcomp$DmWoD), -log(myQQcomp$pscore.S))

plot(myQQcomp$Domesticated, myQQcomp$beta.D)
plot(myQQcomp$Wild, myQQcomp$beta.D)
plot(myQQcomp$DmWoD, myQQcomp$beta.S)

plot(abs(myQQcomp$Domesticated), abs(myQQcomp$beta.D))
plot(abs(myQQcomp$Wild), abs(myQQcomp$beta.W))
plot(abs(myQQcomp$DmWoD), abs(myQQcomp$beta.S))

#try it if we sort on effect size and then plot on location
mybigRR.mt <- mybigRR.mt[order(abs(mybigRR.mt$Domesticated)),]
mybigRR.mt$Drank.B <- 1:nrow(mybigRR.mt)
mybigRR.mt <- mybigRR.mt[order(abs(mybigRR.mt$Wild)),]
mybigRR.mt$Wrank.B <- 1:nrow(mybigRR.mt)
mybigRR.mt <- mybigRR.mt[order(abs(mybigRR.mt$DmWoD)),]
mybigRR.mt$Srank.B <- 1:nrow(mybigRR.mt)

#these two are positively correlated, and kinda sorta linearly so.
plot(-log(myGEMMA$pscore.D), abs(myGEMMA$beta.D))
names(myGEMMA)[9] <- "Index.G"
myGEMMA.mt <- myGEMMA[order(-log(myGEMMA$pscore.D)),]
myGEMMA.mt$Drank.G <- 1:nrow(myGEMMA.mt)
myGEMMA.mt <- myGEMMA.mt[order(-log(myGEMMA.mt$pscore.W)),]
myGEMMA.mt$Wrank.G <- 1:nrow(myGEMMA.mt)
myGEMMA.mt <- myGEMMA.mt[order(-log(myGEMMA.mt$pscore.S)),]
myGEMMA.mt$Srank.G <- 1:nrow(myGEMMA.mt)
#3. pair SNPs  on posrank
myQQcomp.S <- cbind(mybigRR.mt, myGEMMA.mt)
myGEMMA.mt <- myGEMMA.mt[order(myGEMMA.mt$Drank.G),]
mybigRR.mt <- mybigRR.mt[order(mybigRR.mt$Drank.B),]
myQQcomp.D <- cbind(mybigRR.mt, myGEMMA.mt)
myGEMMA.mt <- myGEMMA.mt[order(myGEMMA.mt$Wrank.G),]
mybigRR.mt <- mybigRR.mt[order(mybigRR.mt$Wrank.B),]
myQQcomp.W <- cbind(mybigRR.mt, myGEMMA.mt)
#double check that they're ordered appropriately
plot(myQQcomp.S$Srank.B, myQQcomp.S$Srank.G)
plot(myQQcomp.W$Wrank.B, myQQcomp.W$Wrank.G)
plot(myQQcomp.D$Drank.B, myQQcomp.D$Drank.G)
#4. plot location (GEMMA) vs. location (bigRR) for each phenotype
plot(myQQcomp.D$Index, myQQcomp.D$Index.G)
plot(myQQcomp.W$Index, myQQcomp.W$Index.G)
plot(myQQcomp.S$Index, myQQcomp.S$Index.G)

#nicer plots
jpeg(paste("paper/plots/addGEMMA/SlBc_MAF20_10NA_GEMMA_Domest.jpg", sep=""), width=8, height=5, units='in', res=600)
print(ggplot(myGEMMA, aes(x=Index, y=beta))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(chr),alpha=0.001))+
        labs(list(y="SNP Effect Estimate", title="Domesticated"))+
        theme(legend.position="none")+
        expand_limits(y=0))
dev.off()
