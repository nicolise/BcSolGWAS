#Nicole E Soltis
#10/02/18

#----------------------------------------------------------------------------
#calculate site frequency spectrum/ allele frequency spectrum from SNP information about Bc isolates

#create DNAbin object from VCF
#https://knausb.github.io/vcfR_documentation/dnabin.html

#another AFS/ SFS approach
#https://botany.natur.cuni.cz/hodnocenidat/Lesson_05_tutorial.pdf 

#an R package to analyze SFS 
#https://github.com/andrewparkermorgan/sfsr 

#for this, need sequences as object of class "DNAbin" or "spectrum"
#https://rdrr.io/cran/pegas/man/site.spectrum.html

#another homebrew script for this...
#http://www.molecularecologist.com/wp-content/uploads/2012/03/Allelefrequency_calculations2.txt 


#get vcf
#setwd("~/Projects/BcGenome/data/97_isolates_vcf/")
setwd("~/Documents/GitRepos/BcGenome/data/97_isolates_vcf")
library(vcfR)
myvcf <- read.vcfR("big_set_v97iso_SNPs.vcf")

#trying a filtering step in vcftools
#vcftools --gzvcf raw.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out raw.g5mac3

#vcftools --vcf big_set_v97iso_SNPs.vcf --max-missing 1.0 --mac 1 --recode --recode-INFO-all --out NES_vcf/SNPs.nomiss

#Max-missing is the minimum call-rate allowed so obviously any value greater than 0 wil remove variants and max-missing 0.006 can not give the same results than max-missing 0. ==> to only keep variants with ZERO missing genotypes, max-missing 1.0
#this keeps 1131744 out of a possible 1487721 Sites

#assume that missing calls are a problem for SFS. Try removing any SNPs with missing calls
#later, could impute missing data?
myvcf.nomiss <- read.vcfR("NES_vcf/SNPs.nomiss.recode.vcf")

#example of DNAbin conversion:
#my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, ref.seq = dna[, gff[record, 4]:gff[record, 5]], start.pos = gff[record, 4], verbose = FALSE)

#too big a vector on laptop!
myDNAbin_test <- vcfR2DNAbin(myvcf.nomiss)

#now try sfs
library(pegas)
my_SFS <- site.spectrum(myDNAbin_test, folded = TRUE)
#arbitrarily trying to choose an outgroup (must be an integer)
my_SFS <- site.spectrum(myDNAbin_test, folded = FALSE, outgroup=1)
plot(my_SFS)
#this failed: all = 0

#try just using minor allele frequency!
#plot number of sites vs. minor allele frequency
#element = 1 is major allele, element = 2 is second allele
mymaf2 <- maf(myvcf, element = 2)
gt <- extract.gt(myvcf)
mymafdf_all2 <- as.data.frame(mymaf2)

mymaf1 <- maf(myvcf, element = 1)
mymafdf_all1 <- as.data.frame(mymaf1)

write.csv(mymaf1,"NES_vcf/iso97vcf_mymaf1.csv")
write.csv(mymaf2,"NES_vcf/iso97vcf_mymaf2.csv")

#plot it
library(ggplot2)
hist(mymafdf$Frequency)

#double check from MAF file for consistency
mymaf_frq <- read.csv("MAFs/freq_analysis.csv")
library(stringr)
all1sp <- as.data.frame(str_split_fixed(mymaf_frq$Allele1.F, ":", 2))
names(all1sp) <- c("Allele1.st", "Allele1.frq")
all2sp <- as.data.frame(str_split_fixed(mymaf_frq$Allele2.F, ":", 2))
names(all2sp) <- c("Allele2.st", "Allele2.frq")
all3sp <- as.data.frame(str_split_fixed(mymaf_frq$Allele3.F, ":", 2))
names(all3sp) <- c("Allele3.st", "Allele3.frq")
mymaf_frq_rn <- cbind(mymaf_frq[,c(1:4)], all1sp, all2sp, all3sp)
mymaf_frq_rn$Allele1.frq <- as.numeric(paste(mymaf_frq_rn$Allele1.frq))
mymaf_frq_rn$Allele2.frq <- as.numeric(paste(mymaf_frq_rn$Allele2.frq))
mymaf_frq_rn$Allele3.frq <- as.numeric(paste(mymaf_frq_rn$Allele3.frq))

hist(mymaf_frq_rn$Allele2.frq)
write.csv(mymaf_frq_rn, "NES_vcf/freq_analysis_3alleles.csv")

#--------------------------------------------------------------------
#draw plots on laptop
