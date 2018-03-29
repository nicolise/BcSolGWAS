#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
#first round: just one file at a time. Then convert to loops
#1 is Domest, 2 is Wild, 3 is Sensitivity

#on laptop:
myGEMMA <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_2.assoc.txt", header=TRUE)
#on linux:
myGEMMA <- read.table("GEMMA_lmm_v2/binMAF20NA10_PLINK_kmat1_pheno1.assoc.txt", header=TRUE)

library(ggplot2); 

#virtually identical, just checking.
plot(myGEMMA$p_wald, myGEMMA$p_score)
plot(myGEMMA$p_wald, myGEMMA$p_lrt)
#wald test, likelihood ratio test, or score test

#let's try a manhattan plot. Choosing score test for now.

#create a custom color scale
myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60", "grey20", "grey60")
names(myColors) <- levels(myGEMMA$chr)
colScale <- scale_colour_manual(name = "Chrom",values = myColors)

#sort dataframe rows in order of Chrom, then Pos
str(myGEMMA)
myGEMMA$ps <- as.numeric(myGEMMA$ps)
myGEMMA$chr <- as.numeric(myGEMMA$chr)
myGEMMA <- myGEMMA[with(myGEMMA, order(chr, ps)), ]

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

#comparing positions to B05.10 outputs in 06_transPatterns_B05.R
#pos should be about 0 to 14,000
#ps goes 121 to 408,6062 ... is one of the chromosomes too long somehow?
#Index should be about 1000 to 41,997,650
#Index goes 4224 to 40,421,038
hist(myGEMMA$ps)
hist(myGEMMA$Index)
#positions look fine...

#thresholds: p < 0.05 , p < 0.01, p < 0.001
myp001 <- mean(abs(myGEMMA[myGEMMA$p_score > 0.00099 & myGEMMA$p_score < 0.00101,]$beta), na.rm=TRUE) 
myp01 <- mean(abs(myGEMMA[myGEMMA$p_score > 0.0099 & myGEMMA$p_score < 0.0101,]$beta), na.rm=TRUE) # for 3, beta = 0.0412
#for 1, beta = 0.0673

jpeg(paste("paper/plots/addGEMMA/SlBc_MAF20_10NA_GEMMA_Wild_kmat2_beta.jpg", sep=""), width=8, height=5, units='in', res=600)
  print(ggplot(myGEMMA, aes(x=Index, y=beta))+
#print(ggplot(myGEMMA, aes(x=Index, y=(-log(p_score))))+
          theme_bw()+
          #colScale+
          geom_point(aes(color = factor(chr),alpha=0.001))+
          labs(list(y="-log(p value)", title="Wild"))+
          guides(col = guide_legend(nrow = 8, title="Chromosome"))+
          geom_hline(yintercept=myp001, colour = "black", lty=2) +
          geom_hline(yintercept=-myp001, colour = "black", lty=2) +
          geom_text(aes(0,myp001, label = "p < 0.001", vjust = 1.2, hjust = -0.1), col = "black")+
          geom_hline(yintercept=myp01, lty=2) +
          geom_hline(yintercept=-myp01, lty=2) +
          geom_text(aes(0,myp01, label = "p < 0.01", vjust = 1, hjust = -0.1), col = "black")+
=======
jpeg(paste("paper/plots/addGEMMA/SlBc_MAF20_10NA_GEMMA_Domest_kmat_p.jpg", sep=""), width=8, height=5, units='in', res=600)
  #print(ggplot(myGEMMA, aes(x=Index, y=beta))+
print(ggplot(myGEMMA, aes(x=Index, y=(-log(p_score))))+
          theme_bw()+
          colScale+
          geom_point(aes(color = factor(chr),alpha=0.001))+
          labs(list(y="-log(p value)", title="Domesticated"))+
          guides(col = guide_legend(nrow = 8, title="Chromosome"))+
          #geom_hline(yintercept=myp001, colour = "black", lty=2) +
          #geom_hline(yintercept=-myp001, colour = "black", lty=2) +
          #geom_text(aes(0,myp001, label = "p < 0.001", vjust = 1.2, hjust = -0.1), col = "black")+
          #geom_hline(yintercept=myp01, lty=2) +
          #geom_hline(yintercept=-myp01, lty=2) +
          #geom_text(aes(0,myp01, label = "p < 0.01", vjust = 1, hjust = -0.1), col = "white")+
>>>>>>> 4fb0550a39cf5dbefd32b7db002779179c7257ff
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
