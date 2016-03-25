#Nicole E Soltis
#03/04/16
#converting vcf file to csv for BigRR input
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
remove.packages("RSQLite")
install.packages("RSQLite")
biocLite("GenomicFeatures",lib.loc="/home/nesoltis/R/x86_64-pc-linux-gnu-library/3.0/")
biocLite("VariantAnnotation")
browseVignettes("VariantAnnotation")
library(VariantAnnotation)

#http://www.bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf
#http://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html

view contents of tar
$tar -tf filename.tar.gz

this didn't work
$vim filename.tar.gz

extract the file
$gzip -d filename.tar.gz

$tar -xvf filename.tar

compress the file / directory
$tar -cvzf file.tar.gz inputfile1 inputfile2

read the top few lines
$nano filename.vcf| head

kill a process
Ctrl + Alt + Esc 
kill


