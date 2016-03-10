#Nicole E Soltis
#03/04/16
#converting vcf file to csv for BigRR input

#http://www.bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf
#http://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html
#https://www.biostars.org/p/51076/

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
remove.packages("RSQLite")
install.packages("RSQLite")
biocLite("GenomicFeatures")
#can add option ,lib.loc="/home/nesoltis/R/x86_64-pc-linux-gnu-library/3.0/"
biocLite("VariantAnnotation")
browseVignettes("VariantAnnotation")
library(VariantAnnotation)

setwd("C:/cygwin64/home/nesoltis/data/97_isolates_vcf")
myfile <- system.file("vcf", "big_set_v97iso_SNPs.vcf", package="VariantAnnotation")
myfile <- 
vcf <- readVcf("big_set_v97iso_SNPs.vcf", "mygenomex")
 #508559360 bytes = 0.5 Gb = 509 Mb
memory.limit(size = NA) #16295 Mb = 16.3 Gb
#set memory to at least 16804 Mb
#limit for 64-bit windows is 4Gb if 32-bit R, 8Tb if 64-bit R
memory.limit(size=64000) #32000 was not enough, 64000 is too much for my computer

# vcf <- readVcf("big_set_v97iso_SNPs.vcf", "mygenomex")
# Error: scanVcf: 'Realloc' could not re-allocate memory (508559360 bytes)
# path: C:\cygwin64\home\nesoltis\data\97_isolates_vcf\big_set_v97iso_SNPs.vcf
# In addition: Warning messages:
#   1: In doTryCatch(return(expr), name, parentenv, handler) :
#   Reached total allocation of 16295Mb: see help(memory.size)
# 2: In doTryCatch(return(expr), name, parentenv, handler) :
#   Reached total allocation of 16295Mb: see help(memory.size)
