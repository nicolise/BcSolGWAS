#Nicole E Soltis
#R loop to run to generate randomized phenotypes, then feed into GEMMA
#-------------------------------------------------------------------------

#now randomize each phenotype! Wheee
myFAM_match <- read.csv("GEMMA_files/D_02_randGEMMA/binMAF20NA10_fam.csv")

#test on a small run 
for (i in 1:10){
  Phenos_rand <- transform(Phenos_match, LA0410 = sample(LA410), LA0480 = sample(LA480), LA1547 = sample(LA1547), LA1589 = sample(LA1589), LA1684 = sample(LA1684), LA2093 = sample(LA2093), LA2176 = sample(LA2176), LA2706 = sample(LA2706), LA3008 = sample(LA3008), LA3475 = sample(LA3475), LA4345 = sample(LA4345), LA4355 = sample(LA4355),  Domesticated = sample(Domesticated), Wild = sample(Wild), DmWoD = sample(DmWoD))
  Phenos_rand <- Phenos_rand[,-c(10,13)]
  #only keep randomized phenotypes, not originals
  Phenos_rand <- Phenos_rand[,c(1,15,16,2:14)]
}


#now add Phenos_match onto myFAM_match
myFAM_match2 <- cbind(myFAM_match, Phenos_rand)
myFAM_match2 <- myFAM_match2[order(myFAM_match$delete),]
myFAM_match2 <- myFAM_match2[,c(1:5,11:ncol(myFAM_match2))]

write.table(myFAM_match2, "GEMMA_files/D_02_randGEMMA/binMAF20NA10_randphenos.fam", row.names=FALSE, col.names=TRUE)