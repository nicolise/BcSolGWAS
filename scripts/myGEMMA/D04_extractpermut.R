#Nicole E Soltis
#04/05/18

#-------------------------------------------------------------------------
#list all files in /media/ outputs
#loop through all 1000 files
# for each of 12 phenotypes, create new dataframe that keeps only selected SNPs:
  #df 1: ordered on SNP p-value top 1% = 2720 SNPs ish
  #df 2: ordered on SNP p-value small -> large, take 1st, 2nd, 25th, 250th, 2500th SNP
#save df 1 to media
#save df 2 to ~/Projects
#loop through all df2 to build a new df across all 1000 permutations using rbind
  #columns: chromosome, position, p-value on phenotype 1...12, which SNP (1 / 2 / 25 / 250 / 2500)
  #
