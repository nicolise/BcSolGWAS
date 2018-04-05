#!/bin/sh

#script is run in GEMMA_files/
#files are structured as GEMMA_files/D_05_bigrand/rand1k_1 ... through 3 (but later, 1000)
#and still need a secondary loop through phenotypes in each file.

for i in {1..1000}
do
  mybpath='D_05_bigrand/rand1k_'"$i"'/binMAF20NA10_rand' 
  myoutpath='rand1k_'"$i"
  mkdir -p ./output/$myoutpath;

    for j in {1..15}
    do
        echo "Looping ... phenotype $j"
        ~/gemma/bin/gemma -bfile $mybpath -k D_03_kmat/binMAF20NA10_PLINK_randtest_kmatrix1.cXX.txt -n $j -miss 0.1 -maf 0.2 -lm 4 -o $myoutpath"/pheno${j}"
    done
done
