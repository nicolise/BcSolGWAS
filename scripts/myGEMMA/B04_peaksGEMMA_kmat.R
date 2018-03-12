#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
#first round: just one file at a time. Then convert to loops
#1 is Domest, 2 is Wild, 3 is Sensitivity
#original files, no accounting for pop str
#myGEMMA.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_1.assoc.txt", header=TRUE)
#myGEMMA.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_2.assoc.txt", header=TRUE)
#myGEMMA.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_3.assoc.txt", header=TRUE)

myGEMMA.2.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_1.assoc.txt", header=TRUE)
myGEMMA.2.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_2.assoc.txt", header=TRUE)
myGEMMA.2.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_3.assoc.txt", header=TRUE)

myGEMMA.1.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_1.assoc.txt", header=TRUE)
myGEMMA.1.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_2.assoc.txt", header=TRUE)
myGEMMA.1.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_3.assoc.txt", header=TRUE)

#run through this once for myGEMMA.2, once for myGEMMA.1
names(myGEMMA.2.D)
myGEMMA <- myGEMMA.2.D[,c(1,3,9,13)]
names(myGEMMA)[3] <- "beta.D"
names(myGEMMA)[4] <- "pscore.D"
myGEMMA <- cbind(myGEMMA,myGEMMA.2.W[,c(9,13)])
names(myGEMMA)[5] <- "beta.W"
names(myGEMMA)[6] <- "pscore.W"
myGEMMA <- cbind(myGEMMA,myGEMMA.2.S[,c(9,13)])
names(myGEMMA)[7] <- "beta.S"
names(myGEMMA)[8] <- "pscore.S"
#Make plotting variables
myGEMMA$Index = NA
lastbase = 0
#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(myGEMMA$chr)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    myGEMMA[myGEMMA$chr==i, ]$Index=myGEMMA[myGEMMA$chr==i, ]$ps
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(myGEMMA,myGEMMA$chr==i-1)$ps, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    myGEMMA[myGEMMA$chr==i, ]$Index=myGEMMA[myGEMMA$chr==i, ]$ps+lastbase
  }
}

myGEMMA.fulldat <- myGEMMA
write.csv(myGEMMA.fulldat, "data/GEMMA_files/04_analysis/GEMMA_allDWS.csv")
#select just top SNPs for comparison to bigRR T4
#conditionally replace nonsig values with zero
hist(myGEMMA$beta.D)
myGEMMA$beta.D[myGEMMA$pscore.D > 0.01] <- 0
myGEMMA$beta.W[myGEMMA$pscore.W > 0.01] <- 0
myGEMMA$beta.S[myGEMMA$pscore.S > 0.01] <- 0
#remove rows if all 3 = 0 
myGEMMA_2 <- myGEMMA[!(myGEMMA$beta.D==0 & myGEMMA$beta.W==0 & myGEMMA$beta.S==0),]
#now add counting variable
myGEMMA_2$TotTraits <- ifelse(myGEMMA_2$beta.D != 0 & myGEMMA_2$beta.W != 0 & myGEMMA_2$beta.S != 0, "ALL",
                              ifelse(myGEMMA_2$beta.D != 0 & myGEMMA_2$beta.W != 0, "DW",
                                     ifelse(myGEMMA_2$beta.W != 0 & myGEMMA_2$beta.S != 0, "WS",
                                            ifelse(myGEMMA_2$beta.D != 0 & myGEMMA_2$beta.S != 0, "DS",
                                                   ifelse(myGEMMA_2$beta.D != 0, "D",
                                                          ifelse(myGEMMA_2$beta.W != 0, "W", "S"))))))

table(myGEMMA_2$TotTraits)



jpeg(paste("paper/plots/addGEMMA/SlBc_MAF20_10NA_GEMMA_kmat2_Domest.jpg", sep=""), width=8, height=5, units='in', res=600)
  print(ggplot(myGEMMA, aes(x=Index, y=beta))+
          theme_bw()+
          #colScale+
          geom_point(aes(color = factor(chr),alpha=0.001))+
          labs(list(y="SNP Effect Estimate", title="Domesticated"))+
          guides(col = guide_legend(nrow = 8, title="Chromosome"))+
          geom_hline(yintercept=myp05, colour = "black", lty=2) +
          geom_hline(yintercept=-myp05, colour = "black", lty=2) +
          geom_text(aes(0,myp05, label = "p < 0.05", vjust = 1.2, hjust = -0.1), col = "white")+
          #0.0412 for Sens (3), 0.0673 for Dom (1)
          geom_hline(yintercept=myp01, lty=2) +
          geom_hline(yintercept=-myp01, lty=2) +
          geom_text(aes(0,myp01, label = "p < 0.01", vjust = 1, hjust = -0.1), col = "white")+
          theme(legend.position="none")+
  #same for all 3 phenos
          scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
          expand_limits(y=0))
  dev.off()


#get chromosome midpoints
  my.chroms <- as.data.frame(myGEMMA[!duplicated(myGEMMA$chr, fromLast=FALSE), "Index"]) #Lower Bounds
  names(my.chroms)[1] <- "Chr.Start"
  my.chroms$Chr.End <- myGEMMA[!duplicated(myGEMMA$chr, fromLast=TRUE), "Index"] # Upper Bounds
  my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
  
  
#now read in all 3, make combination file for meta-analysis
