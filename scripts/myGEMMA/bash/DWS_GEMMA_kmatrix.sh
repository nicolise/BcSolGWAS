#!/bin/sh
for i in {1..3}
do
  echo "Looping ... phenotype $i"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -gk 1 -n $i -o binMAF20NA10_PLINK_kmatrix1"_pheno${i}"
  echo "k matrix 2"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -gk 2 -n $i -o binMAF20NA10_PLINK_kmatrix2"_pheno${i}"
done
