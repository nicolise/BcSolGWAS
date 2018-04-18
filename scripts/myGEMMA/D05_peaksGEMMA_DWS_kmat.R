#Nicole E Soltis
#030618

#read in GEMMA outputs with goal of comparison to bigRR for Bc x Solanum GWAS

#--------------------------------------------------------------
rm(list=ls())
setwd("~/Documents/GitRepos/BcSolGWAS/")
setwd("~/Projects/BcSolGWAS/")
#first round: just one file at a time. Then convert to loops

#get thresholds here 
mythrs <- read.csv("data/GEMMA_files/D_07_randOUTS/GEMMA_1krand_thresholds.csv")
mythrs

#check phenotype order here
Phenos_match <- read.csv("data/GEMMA_files/D_02_randGEMMA/binMAF20NA10_fam.csv")
names(Phenos_match)
#also phenos have been renamed in script D05B_manhattanGEMMA_fig4_fig5a.R

#13 is Domest, 14 is Wild, 15 is Sensitivity
#original files, no accounting for pop str
myGEMMA.D <- read.table("data/GEMMA_files/D_04_randphenos/binMAF20NA10_fullrand_kmat1_pheno13_Domest.assoc.txt", header=TRUE)
myGEMMA.W <- read.table("data/GEMMA_files/D_04_randphenos/binMAF20NA10_fullrand_kmat1_pheno14_Wild.assoc.txt", header=TRUE)
myGEMMA.S <- read.table("data/GEMMA_files/D_04_randphenos/binMAF20NA10_fullrand_kmat1_pheno15_Sens.assoc.txt", header=TRUE)
#based on Manhattan plots, kmat1 and kmat2 are identical

names(myGEMMA.D)
myGEMMA <- myGEMMA.D[,c(1,3,9,13)]
names(myGEMMA)[3] <- "beta.D"
names(myGEMMA)[4] <- "pscore.D"
myGEMMA <- cbind(myGEMMA,myGEMMA.W[,c(9,13)])
names(myGEMMA)[5] <- "beta.W"
names(myGEMMA)[6] <- "pscore.W"
myGEMMA <- cbind(myGEMMA,myGEMMA.S[,c(9,13)])
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

mythrs99 <- mythrs[mythrs$SNPnum==2500,]
mythrs999 <- mythrs[mythrs$SNPnum==250,]
curthrs <- mythrs99
thrWi <- curthrs[curthrs$pheno=="Wild",3]
thrDo <- curthrs[curthrs$pheno=="Domesticated",3]
thrSe <- curthrs[curthrs$pheno=="DmWoD",3]

#add a tottraits variable
##need to wait until permutation is done for this step. come back to this
myGEMMA$TotTraits <- ifelse(myGEMMA$pscore.D < thrDo & myGEMMA$pscore.W < thrWi & myGEMMA$pscore.S < thrSe, "ALL",
                              ifelse(myGEMMA$pscore.D < thrDo & myGEMMA$pscore.W < thrWi, "DW",
                                     ifelse(myGEMMA$pscore.W < thrWi & myGEMMA$pscore.S < thrSe, "WS",
                                            ifelse(myGEMMA$pscore.D < thrDo & myGEMMA$pscore.S < thrSe, "DS",
                                                   ifelse(myGEMMA$pscore.W < thrWi, "W", 
                                                          ifelse(myGEMMA$pscore.S < thrSe,"S", 
                                                                 ifelse(myGEMMA$pscore.D < thrDo, "D", "none")))))))

table(myGEMMA$TotTraits)

myGEMMA.fulldat <- myGEMMA
write.csv(myGEMMA.fulldat, "data/GEMMA_files/D_08_results/GEMMA_allDWS_kmat1_99thr.csv")
myGEMMA.fulldat <- read.csv("data/GEMMA_files/D_08_results/GEMMA_allDWS_kmat1_99thr.csv")

myGEMMA <- myGEMMA.fulldat
#select just top SNPs for comparison to bigRR T4
#conditionally replace nonsig values with zero
hist(myGEMMA$beta.D)
myGEMMA$beta.D[myGEMMA$pscore.D > thrDo] <- 0
myGEMMA$beta.W[myGEMMA$pscore.W > thrWi] <- 0
myGEMMA$beta.S[myGEMMA$pscore.S > thrSe] <- 0
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

write.csv(myGEMMA_2, "data/GEMMA_files/D_08_results/GEMMA_peaksDWS_kmat1_99thr.csv")
myGEMMA_2 <- read.csv("data/GEMMA_files/D_08_results/GEMMA_peaksDWS_kmat1_99thr.csv")


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

#ee is orange, 05 is black, 1c is blue
myColors <- c("#1C86EE","#050505", "#EE7600")
names(myColors) <- levels(myGEMMA_3$Phenos)
colScale <- scale_colour_manual(name="Phenotype", values=myColors)

#double check color: low to high group on y axis is Sens, Wild, Domest 
#should be S = black, W = orange, D = blue

#get chromosome midpoints
my.chroms <- as.data.frame(myGEMMA[!duplicated(myGEMMA$chr, fromLast=FALSE), "Index"]) #Lower Bounds
names(my.chroms)[1] <- "Chr.Start"
my.chroms$Chr.End <- myGEMMA[!duplicated(myGEMMA$chr, fromLast=TRUE), "Index"] # Upper Bounds
my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2

jpeg("paper/plots/addGEMMA/S4A_DWSmanhattan.jpg", width=7.5, height=5, units='in', res=600)
ggplot(myGEMMA_3, aes(x=Index))+
  theme_bw()+
  colScale+
  labs(list(y=expression('-log'[10]*'p'), x="Chromosome position", title=element_blank()))+
  guides(col = guide_legend(nrow = 8, title="Phenotype"))+
  geom_point(aes(x=Index, y=-log10(plot.D), color = "Domesticated"), alpha=1/4)+
  geom_point(aes(x=Index, y=-log10(plot.W), color = "Wild"), alpha=1/4)+
  geom_point(aes(x=Index, y=-log10(plot.S), color = "Sensitivity"), alpha=1/4)+
  theme(legend.position="none")+
  theme(text = element_text(size=14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name="Chromosome", breaks = c(2045143, 5763240, 9045566, 11884449, 14590093, 17417481, 20093765, 22716437, 25291433, 27764370, 30138572, 32480630, 34788869, 36988057, 39090468, 40253384), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
#theme(legend.position="none")
dev.off()