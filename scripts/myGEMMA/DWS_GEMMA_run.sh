#!/bin/sh
for i in {1..3}
do
  echo "Looping ... phenotype $i"
  ~/gemma/bin/gemma -bfile B_01_PLINK/binMAF20NA10 -n $i -miss 0.1 -maf 0.2 -lm 4 -o binMAF20NA10_PLINK"_${i}"
done
