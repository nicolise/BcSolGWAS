#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")

#LA2093 is pheno 6 (see Phenos_match)
Phenos_match <- read.csv("data/GEMMA_files/D_02_randGEMMA/binMAF20NA10_fam.csv")

myGEMMA <- read.table("data/GEMMA_files/D_04_randphenos/binMAF20NA10_fullrand_kmat1_pheno6_LA2093.assoc.txt", header=TRUE)

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
#myGEMMA$ps <- as.numeric(myGEMMA$ps)
#myGEMMA$chr <- as.numeric(myGEMMA$chr)
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

#need to add thresholds later? once I do permutation analysis
#thresholds: p < 0.05 , p < 0.01, p < 0.001

#get chromosome midpoints
my.chroms <- as.data.frame(myGEMMA[!duplicated(myGEMMA$chr, fromLast=FALSE), "Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- myGEMMA[!duplicated(myGEMMA$chr, fromLast=TRUE), "Index"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2

jpeg(paste("paper/plots/addGEMMA/SlBc_MAF20_10NA_GEMMArand_LA2093_kmat1_pscore.jpg", sep=""), width=8, height=5, units='in', res=600)
  #print(ggplot(myGEMMA, aes(x=Index, y=beta))+
print(ggplot(myGEMMA, aes(x=Index, y=(-log10(p_score))))+
          theme_bw()+
          colScale+
          geom_point(aes(color = factor(chr),alpha=0.001))+
          labs(list(y=expression('-log'[10]*'p'), title="LA2093"))+
          guides(col = guide_legend(nrow = 8, title="Chromosome"))+
        # geom_hline(yintercept=-log(0.01), colour = "black", lty=2)+
        # geom_hline(yintercept=-log(0.001), colour = "black", lty=2)+
        # geom_text(aes(0,-log(0.001), label="p = 0.001", vjust = 1, hjust = -0.1), col= "black")+
          theme(legend.position="none")+
          theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #same for all 3 phenos
          scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
          expand_limits(y=0)
    # geom_vline(xintercept=11314745, lty=2)+
    # geom_vline(xintercept=36318265, lty=2)+
    # geom_vline(xintercept=29494631, lty=2)+
    # geom_vline(xintercept=15533048, lty=2)+
    # geom_vline(xintercept=33521184, lty=2)+
    # geom_vline(xintercept=31206188, lty=2)+
    # geom_vline(xintercept=16109750, lty=2)+
    # geom_vline(xintercept=23134240, lty=2)+
    # geom_vline(xintercept=40158898, lty=2)+
    # geom_vline(xintercept=37183274, lty=2)+
    # geom_vline(xintercept=10959311, lty=2)
    )
  dev.off()
  
#------------------------------------------------------------------------------ 
## start analysis here once I have permutations
  #now plot for 4b
  setwd("~/Projects/BcSolGWAS/data/GEMMA_files/D_04_randphenos/")
        
  #read in each of 12 phenotypes
  # my.files <- list.files(pattern = c("assoc"))
  
  #check these: 
  names(Phenos_match)
  my.names <- c("1_LA0410", "10_LA3475", "11_LA4345", "12_LA4355", "13_Domest", "14_Wild", "15_Sens", "2_LA0480", "3_LA1547", "4_LA1589","5_LA1684","6_LA2093","7_LA2176","8_LA2706","9_LA3008")
  #rename all files
  # for(i in 1:length(my.files)) {
  #   my.file <- read.csv(my.files[i])
  #   file.rename(from=file.path(my.files[i]), to=file.path(paste(file.path(my.files[i]),my.names[i],".txt",sep="")))
  # }
  
  #now read in files and combine
  my.files <- list.files(pattern = c("assoc"))
  for (i in 1:length(my.files)){
    my.file <- read.table(my.files[i], header=TRUE)
    ifelse(i == 1, full.file <- my.file, full.file <- cbind(full.file, my.file[,c("beta","se","p_score")]))
           ifelse(i == 1, names(full.file)[9] <- paste(my.names[i], "_beta", sep=""), names(full.file)[(ncol(full.file)-2)] <- paste(my.names[i], "_beta", sep=""))
           ifelse(i == 1, names(full.file)[10] <- paste(my.names[i], "_se", sep=""), names(full.file)[(ncol(full.file)-1)] <- paste(my.names[i], "_se", sep=""))
           ifelse(i == 1, names(full.file)[13] <- paste(my.names[i], "_pscore", sep=""), names(full.file)[(ncol(full.file))] <- paste(my.names[i], "_pscore", sep=""))
  }
full.file <- full.file[,-c(2,11,12)]
pheno.bin <- full.file

 #now, PosPhenos = 1 if p < 0.01, 0 if p > 0.01
  ##unclear which phenotypes are significant without permutation analysis -- need to wait on this step until determining permutation threshold.
setwd("~/Projects/BcSolGWAS/data/GEMMA_files")
  for (i in c(1:12)){
    mybeta = 5 + (3*i)
    mypscore = 7 + (3*i)
  pheno.bin[,paste(colnames(pheno.bin[mybeta]),"bin", sep="")] <- ifelse(pheno.bin[,mypscore] < 0.01, 1, 0)
  }

names(pheno.bin)
pheno.bin$SUMM <- rowSums(pheno.bin[,c(44:55)], na.rm=T)

write.csv(pheno.bin, "D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_1kpermut_kmat1.csv")
table(pheno.bin$SUMM)
  
#high overlap SNP list for annotation
HOSNP <- pheno.bin[pheno.bin$SUMM > 6,]
write.csv(HOSNP, "D_08_results/12Plants_HiOverlapSNPs_trueMAF20_10NA_GEMMA_1kpermut_kmat1.csv")

#and top 1000 SNPs > Thr per genotype
#then wide format
top1ksnp <- pheno.bin
#*_betabin columns tells SNP p < 0.01 Thr or not
table(top1ksnp$`11_LA4355_betabin`)
#SUMM column tells which SNPs have no sig phenos
#and number of sig phenos per SNP
#sig phenos only p < 0.01:
top1ksnp <- top1ksnp[top1ksnp$SUMM>0,]
#now only keep rows with top 1000 per phenotype
#remove column 20:22 for failed genotype
top1ksnp <- top1ksnp[,c(1:19,23:43,56,57)]

for (y in c(10,13,16,19,22,25,28,31,34,37,40)){
  #sort data frame with small to high p values for each phenotype
  top1ksnp <- top1ksnp[order(top1ksnp[,y]),]
  #then reassign beta to be snp rank 
  top1ksnp[,(y-2)] <- 1:71754
}


#only beta <= 1000 are significant
top1ksnp$top1k <- ifelse(top1ksnp$`1_LA1547_beta` < 1001 | top1ksnp$`3_LA1684_beta` < 1001 | top1ksnp$`4_LA2093_beta` < 1001 | top1ksnp$`5_LA2176_beta` < 1001 | top1ksnp$`6_LA2706_beta` < 1001 | top1ksnp$`7_LA3008_beta` < 1001 | top1ksnp$`8_LA3475_beta` < 1001 | top1ksnp$`9_LA0410_beta` < 1001 | top1ksnp$`10_LA4345_beta` < 1001 | top1ksnp$`11_LA4355_beta` < 1001 | top1ksnp$`12_LA0480_beta` < 1001, "topset", "omit")
table(top1ksnp$top1k)

top1ksnp <- top1ksnp[top1ksnp$top1k == "topset",]
write.csv(top1ksnp, "D_08_results/12Plants_top1kSNPs_MAF20_10NA_GEMMA_kmat1.csv")

#------------------------------------------------------------------  
pheno.bin <- read.csv("D_08_results/12Plants_allSNPs_MAF20NA10_GEMMA_kmat1.csv")

#make plots 
  setwd("~/Projects/BcSolGWAS")
  jpeg("paper/plots/addGEMMA/FigS1b_sigphenos_ManhattanPlot.jpg", width=8, height=5, units='in', res=600)
  #SUMMtemp <- subset(SUMM.plot[SUMM.plot$Chrom==c(5,6),])
  ggplot(pheno.bin, aes(x=Index, y=SUMM))+
    colScale+ #remove for rainbow plot
    theme_bw()+
    #    scale_x_continuous(breaks = ticks)+
    geom_point(aes(color = factor(chr)))+
    labs(list(y="Plant Genotypes per Significant SNP", x="Chromosome position"))+
    #nrow=8
    theme(legend.position="none")+
    theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
    guides(col = guide_legend(nrow = 8, title=element_blank()))+
    scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
    scale_y_continuous(breaks= c(0,2,4,6,8,10,12), labels=c("0","2","4","6","8","10","12"))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
            geom_vline(xintercept=11314745, lty=2)+
            geom_vline(xintercept=36318265, lty=2)+
            geom_vline(xintercept=29494631, lty=2)+
            geom_vline(xintercept=15533048, lty=2)+
            geom_vline(xintercept=33521184, lty=2)+
            geom_vline(xintercept=31206188, lty=2)+
            geom_vline(xintercept=16109750, lty=2)+
            geom_vline(xintercept=23134240, lty=2)+
            geom_vline(xintercept=40158898, lty=2)+
            geom_vline(xintercept=37183274, lty=2)+
            geom_vline(xintercept=10959311, lty=2)
  dev.off()
  
#find top SNPs to mark: pheno.bin$SUMM and LA2093 pscore
  pheno.bin.top <- pheno.bin[pheno.bin$`6_LA2093_pscore` < 0.001, ]
  #top 100 beta out of these
  library(plyr)
  pheno.bin.top <- head(arrange(pheno.bin.top, desc(`6_LA2093_beta`)), n=100)
  pheno.bin.top <- pheno.bin.top[pheno.bin.top$SUMM>6,]
pheno.bin.top$Index

#add figure 5a / b for GEMMA (S3): SNP overlap across phenotypes
table(pheno.bin$SUMM)

jpeg("paper/plots/addGEMMA/S3a_topSNPOverlap_12Plants_GEMMA.jpg", width=8, height=5, units='in', res=600)
ggplot(pheno.bin, aes(pheno.bin$SUMM)) + 
  geom_bar()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(name= "Number of SNPs", limits = c(0,62000))+
  scale_x_continuous(name= "Plant Genotypes per Candidate SNP", breaks=c(1,2,3,4,5,6,7,8,9,10,11,12),labels=c(1,2,3,4,5,6,7,8,9,10,11,12), limits = c(0, 12))+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
dev.off()


#small inset plot
jpeg("paper/plots/addGEMMA/S3a_topSNPOverlap_12Plants_GEMMA_inset.jpg", width=4, height=3, units='in', res=600)
ggplot(pheno.bin, aes(pheno.bin$SUMM)) + 
  geom_bar()+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(name= "", limits = c(0,250))+
  scale_x_continuous(breaks=c(6,7,8,9,10,11,12),labels=c(6,7,8,9,10,11,12), limits = c(5, 13), name="")+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
dev.off()

setwd("~/Projects/BcSolGWAS/paper/plots/addGEMMA")
#can add this later if needed. calculation needs fixing.
myprobs <- read.csv("GEMMA_probabilities.csv")

for (i in 44:55){
  myg <- sum(pheno.bin[,i])
  print (i)
  print(myg)}