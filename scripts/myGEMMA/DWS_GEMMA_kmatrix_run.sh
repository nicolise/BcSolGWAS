#!/bin/sh
for i in {1..3}
do
  echo "Looping ... phenotype $i"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -k output/binMAF20NA10_kmatrix1.cXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_PLINK_kmat1"_${i}"
  echo "k matrix 2"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -k output/binMAF20NA10_kmatrix2.sXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_PLINK_kmat2"_${i}"
done
