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
#based on Manhattan plots, kmat1 and kmat2 are identical
myGEMMA.2.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_1.assoc.txt", header=TRUE)
myGEMMA.2.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_2.assoc.txt", header=TRUE)
myGEMMA.2.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat2_3.assoc.txt", header=TRUE)

myGEMMA.1.D <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_1.assoc.txt", header=TRUE)
myGEMMA.1.W <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_2.assoc.txt", header=TRUE)
myGEMMA.1.S <- read.table("data/GEMMA_files/03_GEMMAouts/binMAF20NA10_PLINK_kmat1_3.assoc.txt", header=TRUE)

#run through this once for myGEMMA.2, once for myGEMMA.1
names(myGEMMA.1.D)
myGEMMA <- myGEMMA.1.D[,c(1,3,9,13)]
names(myGEMMA)[3] <- "beta.D"
names(myGEMMA)[4] <- "pscore.D"
myGEMMA <- cbind(myGEMMA,myGEMMA.1.W[,c(9,13)])
names(myGEMMA)[5] <- "beta.W"
names(myGEMMA)[6] <- "pscore.W"
myGEMMA <- cbind(myGEMMA,myGEMMA.1.S[,c(9,13)])
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

#add a tottraits variable
myGEMMA$TotTraits <- ifelse(myGEMMA$pscore.D < 0.01 & myGEMMA$pscore.W < 0.01 & myGEMMA$pscore.S <0.01, "ALL",
                              ifelse(myGEMMA$pscore.D < 0.01 & myGEMMA$pscore.W < 0.01, "DW",
                                     ifelse(myGEMMA$pscore.W < 0.01 & myGEMMA$pscore.S <0.01, "WS",
                                            ifelse(myGEMMA$pscore.D < 0.01 & myGEMMA$pscore.S <0.01, "DS",
                                                   ifelse(myGEMMA$pscore.W < 0.01, "W", 
                                                          ifelse(myGEMMA$pscore.S < 0.01,"S", 
                                                                 ifelse(myGEMMA$pscore.D < 0.01, "D", "none")))))))

table(myGEMMA$TotTraits)

myGEMMA.fulldat <- myGEMMA
write.csv(myGEMMA.fulldat, "data/GEMMA_files/04_analysis/GEMMA_allDWS_kmat1.csv")
myGEMMA.fulldat <- read.csv("data/GEMMA_files/04_analysis/GEMMA_allDWS_kmat1.csv")
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

write.csv(myGEMMA_2, "data/GEMMA_files/04_analysis/GEMMA_peaksDWS_kmat1.csv")
myGEMMA_2 <- read.csv("data/GEMMA_files/04_analysis/GEMMA_peaksDWS_kmat1.csv")


#-------------------------------------------------------------------------------------
#manhattan plot
#meta analysis plot

#keep only top 1000 SNPs/ phenotype

myGEMMA_3 <- myGEMMA_2[order(myGEMMA_2$pscore.D),]
myGEMMA_3$Drank <- 1:nrow(myGEMMA_3)
myGEMMA_3 <- myGEMMA_3[order(myGEMMA_3$pscore.W),]
myGEMMA_3$Wrank <- 1:nrow(myGEMMA_3)
myGEMMA_3 <- myGEMMA_3[order(myGEMMA_3$pscore.S),]
myGEMMA_3$Srank <- 1:nrow(myGEMMA_3)

myGEMMA_3$tops <- ifelse(myGEMMA_3$Drank < 1001 | myGEMMA_3$Wrank < 1001 | myGEMMA_3$Srank < 1001, "keep", "omit")

myGEMMA_3 <- myGEMMA_3[myGEMMA_3$tops == "keep",]

myGEMMA_3$plot.D <- ifelse(myGEMMA_3$Drank < 1001, myGEMMA_3$pscore.D, NA)
myGEMMA_3$plot.W <- ifelse(myGEMMA_3$Wrank < 1001, myGEMMA_3$pscore.W, NA)
myGEMMA_3$plot.S <- ifelse(myGEMMA_3$Srank < 1001, myGEMMA_3$pscore.S, NA)

myColors <- c("#050505", "#1C86EE", "#EE7600")
#names(myColors) <- levels(HEM.plotdata$Phenos)
colScale <- scale_colour_manual(name="Phenotype", values=myColors)

jpeg("paper/plots/addGEMMA/S4A_DWSmanhattan.jpg", width=7.5, height=5, units='in', res=600)
ggplot(myGEMMA_3, aes(x=Index))+
  theme_bw()+
  colScale+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title="Phenotype"))+
  geom_point(aes(x=Index, y=-log(plot.D), color = "Domesticated"), alpha=1/4)+
  geom_point(aes(x=Index, y=-log(plot.W), color = "Wild"), alpha=1/4)+
  geom_point(aes(x=Index, y=-log(plot.S), color = "Sensitivity"), alpha=1/4)+
  geom_hline(yintercept=-log(0.01), colour = "black", lty=2)+
  geom_hline(yintercept=-log(0.001), colour = "black", lty=2)+
  geom_text(aes(0,-log(0.001), label="p = 0.001", vjust = 1, hjust = -0.1), col= "black")+
  geom_text(aes(0,-log(0.01), label="p = 0.01", vjust = 1, hjust = -0.1), col= "black")+
  theme(legend.position="none")+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
#theme(legend.position="none")
dev.off()