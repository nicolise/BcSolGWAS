#Nicole E Soltis
#10/02/18

#----------------------------------------------------------------------------
#calculate site frequency spectrum/ allele frequency spectrum from SNP information about Bc isolates

#an R package to analyze SFS 
#https://github.com/andrewparkermorgan/sfsr 

#https://rdrr.io/cran/pegas/man/site.spectrum.html
#for this, need sequences as object of class "DNAbin" or "spectrum"

#create DNAbin object from VCF
#https://knausb.github.io/vcfR_documentation/dnabin.html

#get vcf
setwd("~/Projects/BcGenome/data/97_isolates_vcf/")
library(vcfR)
myvcf <- read.vcfR("big_set_v97iso_SNPs.vcf")
#my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, ref.seq = dna[, gff[record, 4]:gff[record, 5]], start.pos = gff[record, 4], verbose = FALSE)

#too big a vector!
myDNAbin_test <- vcfR2DNAbin(myvcf)
