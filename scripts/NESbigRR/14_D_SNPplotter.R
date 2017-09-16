#Nicole E Soltis
#081117
#plot of SNPs along gene of interest

#12_singleGeneManhattan.R
#-------------------------------------------------------------------------
#START HERE, ACTUALLY
#now to draw that darn LD plot
#this requires 3 file formats, none of which I have
library("snp.plotter")
setwd("~/Documents/GitRepos/BcSolGWAS/data/genome/chr16_analysis/")

snp.plotter(config.file = "Chr16_config.txt")
#----------------------------------------------------------------------------------
#some stuff I tried before (archive)

#this might make DNAbin files, from vcf files
#let's try just a 4kb region on Chromosome 16
#goes from POS 341 to POS 664510
#try POS 344785 to POS 347542
library("vcfR")
vcf <- read.vcfR("data/genome/chr16_analysis/chr16_analysis_seg.recode.vcf", verbose = FALSE)
my_dnabin1 <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE, verbose = FALSE)
my_dnabin1
#plot it by site
ape::image.DNAbin(my_dnabin1[,ape::seg.sites(my_dnabin1)])
#save it as a fasta
write.dna( my_dnabin1, file = 'data/genome/chr16_analysis/chr16_analysis_seg.fasta', format = 'fasta' )
#this might make haplotype files, from DNAbin files
#pegas fails on ubuntu -- run on windows PC
library("pegas")
my_haplo <- haplotype(my_dnabin1)
plot(sort(my_haplo))
attr(my_haplo, "index")[[2]]
print(my_haplo)
my_haplo_dist <- dist.haplotype.loci(my_haplo)
?haplotype()

subset(my_haplo, minfreq=20)
as.data.frame(my_haplo)
