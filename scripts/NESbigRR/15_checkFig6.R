#Nicole E Soltis
#01/07/19
#checking figure 6a SNP block limits -- esp block 2, 3
#----------------------------------------------------------

#limits from 12_singleGeneManhattan...

#exon 1 825306	826178
geom_rect(mapping=aes(ymin=-0.001, ymax=0.001, xmin=825306, xmax=826178), alpha=0.01, fill="darkturquoise")+
  #exon 2 826235	826345
  geom_rect(mapping=aes(ymin=-0.001, ymax=0.001, xmin=826235, xmax=826345), alpha=0.01, fill="darkturquoise")+
  #block 1
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=823323, xmax=823506), alpha=0.01, fill="red")+
  #block 2
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=823848, xmax=824176), alpha=0.01, fill="red")+
  #block 3
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=824908, xmax=827148), alpha=0.01, fill="red")+
  #block 4
  geom_rect(mapping=aes(ymin=0.010, ymax=0.012, xmin=827641, xmax=828305), alpha=0.01, fill="red")+
  
  
  #SNPs from 14_A_TABtoPED_match8a_chr2...
  setwd("~/Projects/BcSolGWAS")
#original MAP file is faulty...
myMAP <- read.table("data/Bcgenome/chr2_analysis/PLINK/myCHR2_A.map")

#weird thing to check: 
#21  2 snp21  0 825228
#22  2 snp22  0  82555
#23  2 snp23  0 826119
#24  2 snp24  0 826250
#25  2 snp25  0   8264
#26  2 snp26  0  82642
#27  2 snp27  0 826529
#are these SNP position numbers losing zeros and still ordered correctly, or does the range incorrectly include SNPs from a different location on the chromosome?
#it incorrectly includes SNPs from earlier on the chromosome, due to the error of sorting on POS as a character vector rather than numeric. This is corrected in 14_A_TABtoPED_...EDIT.R

#redrew figure 6b with the files generated in 14_A_TABtoPED_...EDIT.R
#the method in haploview is Linkage format
#the input files are 
#Data file: data/BcGenome/chr2_analysis/PLINK/myCHR2_A.nullPheno_EDIT_handEd.ped
    #handEd is FAM before Isolate, and quotes added
#Locus information file: data/BcGenome/chr2_analysis/PLINK/myCHR2_A_EDIT.map

#figure blocks:

#block 1: snp 1:4
#from edited MAP, snp 1 = 823323, snp 4 = 823506
#on plot, xmin=823323, xmax=823506

#block 2: snp 5:11
#from MAP, snp 5 = 823848, snp 11 = 824176
#on plot, xmin=823848, xmax=824176

#block 3: snp 13:26
#from MAP, snp 13 = 824908, snp 26 = 827148
#on plot, xmin=824908 (snp 15), xmax=827148 (snp 31)

#block 4: snp 28:32
#from MAP, snp 28 = 827641, snp 32 = 828305
# on plot, xmin=827641 (snp 33), xmax=828305 (snp 37)

#--------------------------------------
#HAPLOVIEW help
#from BcAtRNAGWAS, used linkage format
#data file: C:\Program Files (x86)\HaploView\NESfiles\binMAF20NA10_chr1_boa.ped
#locus information file: same, but .info

#here,
#BcSolGWAS/data/BcGenome/chr2_analysis/
#LSplot2.2V2.svg is the final image 6b // with labels is LDplot2.2V1.lg.png
#PLINK
