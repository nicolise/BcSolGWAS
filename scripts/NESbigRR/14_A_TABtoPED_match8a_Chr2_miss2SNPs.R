#Nicole E Soltis
#091317
#plot of SNPs along gene of interest
#and now: haplotype plots
#convert .tab SNP data to binary .csv

#------------------------------------------------------
rm(list=ls())
#on linux desktop
setwd("~/Documents/GitRepos/BcSolGWAS/")
#on windows laptop
setwd("~/Projects/BcSolGWAS/")
#convert .tab SNP file to .csv
#sticking to MAF20, 10% missingness

#convert each file to binary
SNPsMAF20 <- read.csv("data/GEMMA_files/02_csvPrep/snps_maf20.csv")
mySNPs <- SNPsMAF20

mySNPs[] <- lapply(mySNPs, as.character)
mySNPs[mySNPs=="."]<-NA #this is a true NA
allSNPs<- mySNPs
mySNPs_Chr2 <- mySNPs[grep("Chromosome2.2", mySNPs$X.CHROM), ]
mySNPs_Chr2 <- mySNPs_Chr2[order(mySNPs_Chr2$POS),]

#remove low MAFs!
names(mySNPs_Chr2)

mySNPs_Chr2$Freq.A <- rowSums(mySNPs_Chr2[,c(4:100)] =="A", na.rm=T)
mySNPs_Chr2$Freq.C <- rowSums(mySNPs_Chr2[,c(4:100)] =="C", na.rm=T)
mySNPs_Chr2$Freq.T <- rowSums(mySNPs_Chr2[,c(4:100)] =="T", na.rm=T)
mySNPs_Chr2$Freq.G <- rowSums(mySNPs_Chr2[,c(4:100)] =="G", na.rm=T)
mySNPs_Chr2$All <- (mySNPs_Chr2$Freq.A + mySNPs_Chr2$Freq.C + mySNPs_Chr2$Freq.G + mySNPs_Chr2$Freq.T)
#97 isolates total
#remove all loci (rows) with >10% missingness (mySNPs_Chr16$All <87.3)
mySNPs_Chr2 <- mySNPs_Chr2[mySNPs_Chr2$All>87,] 

#now remove all loci (rows) with MAF <20 (max > 77.6)
mySNPs_Chr2 <- transform(mySNPs_Chr2, max = pmax(Freq.A, Freq.C, Freq.G, Freq.T))
mySNPs_Chr2 <- mySNPs_Chr2[mySNPs_Chr2$max < 78,] #this is fine, keep all.
mySNPs_Chr2$missing <- apply(mySNPs_Chr2[,c(4:100)], 1, function(x) sum(is.na(x)))

#MAF is too low if: max + missing > 78
#removes a few more
mySNPs_Chr2 <- mySNPs_Chr2[(mySNPs_Chr2$max + mySNPs_Chr2$missing) < 78,]

#NOW filter to just region of interest on Chromosome 2
#823306 to 828345
mySNPs_Chr2 <- mySNPs_Chr2[order(mySNPs_Chr2$POS),]
mySNPs_Chr2 <- mySNPs_Chr2[which(mySNPs_Chr2$POS > 823305),]
mySNPs_Chr2 <- mySNPs_Chr2[which(mySNPs_Chr2$POS < 828346),]

#originally: 544 SNPs here. But need to re-filter to match figure 8a
#filter to only keep SNPs from Fig8a: 12_singleGeneManhattan_figR8.R
mySNPlist_8a <- read.csv("data/genome/chr2_analysis/SNPlistFig8a.csv")
mySNPs_Chr2.x <- mySNPs_Chr2[mySNPs_Chr2$POS %in% mySNPlist_8a$x, ]


#and remove any duplicated POS here
mySNPs_Chr2 <- mySNPs_Chr2[!duplicated(mySNPs_Chr2$POS), ]
mySNPs_Chr2 <- mySNPs_Chr2[order(mySNPs_Chr2$POS),]
#excellent, now this matches the SNPs file

#save a list of SNPs to match with 14_B_SNP.FILE.R where I'm writing SNP.FILE
mySNPlist <- as.data.frame(mySNPs_Chr2$POS)

write.csv(mySNPlist, "data/genome/chr2_analysis/SNPlistfromPED.csv")

#and now for making PED format for PLINK!
  #do not need positional info: just SNP states for PED
#turn df sideways (individuals as rows, SNPs as columns)
#split each genotype into 2 identical columns (PED assumes diploid)
#keeping REF variable for now, to make binary phenotype for GENOTYPE.FILE
mySNPs_Chr2_2 <- mySNPs_Chr2[,-c(1:2, 101:107)]
#write.csv(mySNPs_Chr16_2, "GEMMA_files/02_csvPrep/hp_charMAF20_10NA_forPED.csv")
#df2[,c(1,3,2,4)]

#get SNP names for matching figure 8a to figure 8b
mydf <- as.data.frame(mySNPs_Chr2$POS)
names(mydf)[1] <- "Pos"
mydf$SNPnum <- c(1:32)
write.csv(mydf, "data/genome/chr2_analysis/PLINK/MatchDrawLines.csv")

#duplicate all the rows, to fake a diploid genome
#this is slow (but only takes ~4 seconds if 600 SNPs)
Sys.time()
mySNPs_Chr2_3 <- as.data.frame(NULL)
for(i in 1:nrow(mySNPs_Chr2_2)){
  mySNPs_Chr2_3<-rbind(mySNPs_Chr2_3, mySNPs_Chr2_2[i,])
  mySNPs_Chr2_3<-rbind(mySNPs_Chr2_3, mySNPs_Chr2_2[i,])
}
Sys.time()

#transpose and format for PED
mySNPs_Chr2_4 <- as.data.frame(t(mySNPs_Chr2_3))
#add binary phenotype = 1 (6)
mySNPs_Chr2_4 <- cbind("Pheno" = 1, mySNPs_Chr2_4)
#add individual sex = 1 (5)
mySNPs_Chr2_4 <- cbind("sex" = 1, mySNPs_Chr2_4)
#add Mother = 0 (4)
mySNPs_Chr2_4 <- cbind("Mother" = 0, mySNPs_Chr2_4)
#add Father = 0 (3)
mySNPs_Chr2_4 <- cbind("Father" = 0, mySNPs_Chr2_4)
#turn row names into column 2
mySNPs_Chr2_4 <- cbind(rownames(mySNPs_Chr2_4), mySNPs_Chr2_4)
colnames(mySNPs_Chr2_4)[1] <- 'Isolate'
#add the fam column (1)
mySNPs_Chr2_4 <- cbind("FAM" = "FAM1", mySNPs_Chr2_4)
#remove ref for myPED
myPED <- mySNPs_Chr2_4[-c(1),] 

#add a real phenotype to PED
#check same phenos as 04_runbigRR_indplants.R
myPhenos <- read.csv("data/GWAS_files/03_bigRRinput/NewModel0711/Sl_Pheno_bigRR_trueMAF20_10NA.csv")

SNPnames <- read.csv("data/GWAS_files/02_csvPrep/Key_SNPnames.csv")
SNPnames <- SNPnames[c(2,5)]

names(SNPnames)[1] <- "Igeno"
names(SNPnames)[2] <- "Isolate"

myPhenos <- merge(myPhenos,SNPnames, by="Igeno")
addPhenos <- myPhenos[,c("Isolate","LA1547")]

#match files on myFAM$V2 (original )
#need to keep rows of myPED with no pheno match!
myPED.2 <- merge(myPED,addPhenos, by="Isolate", all=TRUE)
#fill in missing phenos with average
mean(myPED.2$LA1547, na.rm=TRUE)
#move phenotype column and remove dummy
myPED.2 <- myPED.2[,c(1:5,71,7:70)]
myPED.2$LA1547[is.na(myPED.2$LA1547)] <- mean(myPED.2$LA1547, na.rm=TRUE)
#replace all NA with 0
#first add 0 as a valid level
for (i in 7:70){
  levels(myPED.2[,i]) <- c(levels(myPED.2[,i]),0)
}
myPED.2[is.na(myPED.2)] <- 0

myPED.bin <- myPED.2
myPED.bin$LA1547 <- ifelse( myPED.bin$LA1547 < mean(myPED.2$LA1547), 0, 1)
myPED.null <- myPED.2
myPED.null$LA1547 <- 1
#myPED.null <- myPED.null[,c(1:5,71,6:70)]

#finally, make the GENOTYPE.FILE for SNPplotter
#transposed and formatted for modified PED
myGENOTYPE.FILE <- mySNPs_Chr2_4
rownames(myGENOTYPE.FILE)
# family ID, individual ID, father ID, mother ID, sex, and affection status followed by marker loci coded as binary factors, as shown in the example below. This file should not have column headers.
#make marker loci binary

#by column
for (i in 7:ncol(myGENOTYPE.FILE)){
  #first have to add levels to the factor!!
  levels(myGENOTYPE.FILE[,i]) <- c(levels(myGENOTYPE.FILE[,i]),1,0)
  #by row
  for (j in 2:nrow(myGENOTYPE.FILE)){
    myGENOTYPE.FILE[j,i] <- ifelse(myGENOTYPE.FILE[j,i]==myGENOTYPE.FILE[1,i], 0,1)
  }
}
#now remove column headers and REF
myGENOTYPE.FILE <- myGENOTYPE.FILE[-c(1),]

#make a MAP file for plink (need it to make the bed (binary ped) file from ped)
myMAP <- mySNPs_Chr2[,c("X.CHROM","POS")]
#remove "Chromosome" from X.CHROM
myMAP <- as.data.frame(lapply(myMAP, function(x) {
                 gsub("Chromosome", "", x)
              }))

#1. Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
#2. Variant identifier
#3. Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
#4. Base-pair coordinate
myMAP$SNPid <- paste("snp",c(1:nrow(myMAP)), sep="")
myMAP$cm <- 0
myMAP <- myMAP[,c(1,3,4,2)]
myMAP$X.CHROM <- 2 #can't use 2.2 for HaploView

#----------------------------------------------------------------------
#save all of the files
#MAP
write.table(myMAP, file="data/genome/chr2_analysis/PLINK/myCHR2_A.map", quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

#PED
write.table(myPED.2, file="data/genome/chr2_analysis/PLINK/myCHR2_A.ped", quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
write.table(myPED.null, file="data/genome/chr2_analysis/PLINK/myCHR2_A.nullPheno.ped", quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
write.table(myPED.bin, file="data/genome/chr2_analysis/PLINK/myCHR2_A.binPheno.ped", quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

#GENOTYPE.FILE
write.table(myGENOTYPE.FILE, file="data/genome/chr2_analysis/GENOTYPE.FILE_8aMatch.txt",quote=FALSE, sep=' ',col.names=FALSE,row.names=FALSE)
