#Nicole E Soltis
#030818

#B05_QQ_GEMMAbigRR.R
#---------------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#original file, no accounting for pop str
#only doing this with kmat1 because the results are identical
myGEMMA <- read.csv("data/GEMMA_files/D_08_results/GEMMA_allDWS_kmat1_99thr.csv")
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
table(myGEMMA$TotTraits)

#for comparisons: T4 and B05.10 have different alignments & different chromosome numbers
#CAN rank all SNPs by absolute value of effect size and plot them side by side.
#but mybigRR has an extra 30,000 SNPs- how to compare? 
#qqplot can deal with differing numbers of observations per trait. Do not need to align data frames.
#If trying to align, could go off of Index. Should approximately align genomes in order

#and for the qq plot:
qqplot(myGEMMA$pscore.D, mybigRR$Domesticated)
qqplot(myGEMMA$pscore.W, mybigRR$Wild)
qqplot(myGEMMA$pscore.S, mybigRR$DmWoD)

qqplot(myGEMMA$beta.D, mybigRR$Domesticated)
qqplot(myGEMMA$beta.W, mybigRR$Wild)
qqplot(myGEMMA$beta.S, mybigRR$DmWoD)

#other options:
qplot

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
