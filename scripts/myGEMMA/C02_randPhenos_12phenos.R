#Nicole E Soltis 
#03/28/18
#randomize phenotypes once each to try out permutations for GEMMA

#------------------------------------------------------------------------------
rm(list=ls())

#advice on permutations here:
#https://github.com/genetics-statistics/GEMMA/issues/93

#pipeline note:
#1. run B01_TABtoPEDnMAP.R
#2. copy plink to GEMMA_files
#3. in command prompt: cd to GEMMA_files
#4. RUN ./plink --noweb --file B_01_PLINK/dpbinMAF20NA10 --maf 0.2 --make-bed --out binMAF20NA10 } did this ONCE for all BcSolGWAS. NEXT STEP is customized by DWS/ 12phenos
#5. run this script (B01_12p_FAMaddPhenos.R)
#6. cd to GEMMA_files
#7. copy edited .fam, .bim, .bam to randtest/
#8. calculate k-matrix with: bash randphenos_GEMMA_kmatrix.sh
#9. run GEMMA: bash randphenos_GEMMA_kmatrix_run.sh

#10: Pheno order is "V2"           "LA0410"       "LA0480"       "LA1547"       "LA1589"
#"LA1684"       "LA2093"       "LA2176"       "LA2706"       "LA3008"      
#"LA3475"       "LA4345"       "LA4355"       "Domesticated" "Wild"   "DmWoD"

setwd("~/Projects/BcSolGWAS/data")
setwd("~/Documents/GitRepos/BcSolGWAS/data/")
Phenos <- read.csv("GWAS_files/02_csvPrep/phenos/NewModel0711/BcSl_lsmeans_forbigRR.csv")
D.Phenos <- read.csv("GWAS_files/02_csvPrep/phenos/Domestication/BcSl_lsmeans_domest_forbigRR.csv")
#myFAM is the PLINK output of converting *.ped and *.map (01_TABtoPEDnMAP.R) to *.bed and *.bim and *.fam
myFAM <- read.table("GEMMA_files/B_01_PLINK/binMAF20NA10.fam")
#GEMMA only needs column 1 (individual ID), column 6 (phenotype)
#n (num) with -n 1 means column 6 is pheno, -n 2 means 7… etc.

#first merge Phenos, D.Phenos
setdiff(Phenos$Igeno, D.Phenos$Igeno) #cool
myPhenos <- cbind(Phenos,D.Phenos)
myPhenos <- myPhenos[,-c(14)]
Phenos <- myPhenos
#col2 = V2 = Isolate
names(Phenos)[1] <- "V2"

myFAM_match <- myFAM
Phenos_match <- Phenos

setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)

#add an empty variables 1.01.12 and 1.02.19 and Triple7 to Phenos
for (y in 2:length(Phenos_match)){
  Phenos_match[96,y] <- mean(Phenos_match[,y], na.rm=TRUE)
  Phenos_match[97,y] <- mean(Phenos_match[,y], na.rm=TRUE)
  Phenos_match[98,y] <- mean(Phenos_match[,y], na.rm=TRUE)
}
levels(Phenos_match$V2) <- c(levels(Phenos_match$V2), "1.01.12", "1.02.19","Triple7", "MEAPGG")
Phenos_match[96,1] <- "1.01.12"
Phenos_match[97,1] <- "1.02.19"
Phenos_match[98,1] <- "Triple7"

#change Phenos_match$V2 "MEAP6G" to "MEAPGG"
Phenos_match[80, 1] <- "MEAPGG"

#remove 1.02.05, 94.4
Phenos_match <- Phenos_match[c(1:10,12:94,96:98),]

#reorder phenos
Phenos_match$V2
Phenos_match <- Phenos_match[c(1:5,94,6:15,95,16:44,93,45:91,96,92),]

#double check now and match orders
myFAM_match$delete <- c(1:97)
myFAM_match <- myFAM_match[ order(myFAM_match$V2), ]
setdiff(myFAM_match$V2,Phenos_match$V2)
setdiff(Phenos_match$V2,myFAM_match$V2)
intersect(myFAM_match$V2,Phenos_match$V2)
myFAM_match <- myFAM_match[1:96,]

#now randomize each phenotype! Wheeee
Phenos_rand <- transform(Phenos_match, LA1547 = sample(LA1547), LA1589 = sample(LA1589), LA1684 = sample(LA1684), LA2093 = sample(LA2093), LA2176 = sample(LA2176), LA2706 = sample(LA2706), LA3008 = sample(LA3008), LA3475 = sample(LA3475), LA0410 = sample(LA410), LA4345 = sample(LA4345), LA4355 = sample(LA4355), LA0480 = sample(LA480), Domesticated = sample(Domesticated), Wild = sample(Wild), DmWoD = sample(DmWoD))
Phenos_rand <- Phenos_rand[,-c(10,13)]
Phenos_rand <- Phenos_rand[,c(1,15,16,2:14)]

#now add Phenos_match onto myFAM_match
myFAM_match2 <- cbind(myFAM_match, Phenos_rand)
myFAM_match2 <- myFAM_match2[order(myFAM_match$delete),]
myFAM_match2 <- myFAM_match2[,c(1:5,11:ncol(myFAM_match2))]

write.table(myFAM_match2, "GEMMA_files/B_01_PLINK/ind12plants/binMAF20NA10_randphenos.fam", row.names=FALSE, col.names=TRUE)
