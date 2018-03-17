#!/bin/sh
for i in {1..3}
do
  echo "Looping ... phenotype $i"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -k 04_GEMMAoutput/kmatrices/binMAF20NA10_PLINK_kmatrix1_pheno"${i}".cXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_PLINK_kmat1"_pheno${i}"
  echo "k matrix 2"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -k 04_GEMMAoutput/kmatrices/binMAF20NA10_PLINK_kmatrix2_pheno"${i}".sXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_PLINK_kmat2"_pheno${i}"
done
