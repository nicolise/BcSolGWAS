#!/bin/sh
  echo "only running for default phenotype: DWS told us the k-matrices for all phenos are same!"
  echo "k matrix 1"
  ~/gemma/bin/gemma -bfile B_01_PLINK/ind12plants/binMAF20NA10 -gk 1 -o binMAF20NA10_PLINK_ind12plant_kmatrix1
  echo "k matrix 2"
  ~/gemma/bin/gemma -bfile B_01_PLINK/ind12plants/binMAF20NA10 -gk 2 -o binMAF20NA10_PLINK_ind12plant_kmatrix2
done
