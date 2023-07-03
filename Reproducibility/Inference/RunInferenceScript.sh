#!/bin/bash

# Example: Run the inference scripts for all the HepG2 and K562 scRNA-seq datasets
# using the data with 500 HVGs + eCLIP RBPs provided in the inputs folder

#SETS=('TF_RNA500' 'TF_RNA1000' 'RBP_RNA500' 'RBP_RNA1000')

REGSETS=('RBP_RNA500')
DATASETS=('HepG2_Smartseq2' 'HepG2_9CL_SCAN' 'K562_CELseq' 'K562_STORMseq' 'K562_Smartseq3' 'K562_9CL_SCAN' 'K562_UMI200_SCAN')

for i in "${REGSETS[@]}";
   do
   for j in "${DATASETS[@]}";
      do
      python InferencePipeline.py --folder ./ --dataID ${j}_$i --norm ${j}NormalizedData_$i.csv --raw ${j}RawData_$i.csv --pseudo ${j}_PseudoTime.csv --select 100
      mkdir ./outputs/${j}_$i/DeePSEM/
      mkdir ./outputs/${j}_$i/ARACNe/
   done
done

# HepG2 DNBelab has more than 1k cells, we run TENET using the 20% of the cells
DATASETS=('HepG2_DNBelab')

for i in "${REGSETS[@]}";
   do
   for j in "${DATASETS[@]}";
      do
      python InferencePipeline.py --folder ./ --dataID ${j}_$i --norm ${j}NormalizedData_$i.csv --raw ${j}RawData_$i.csv --pseudo ${j}_PseudoTime.csv --select 20
      mkdir ./outputs/${j}_$i/DeePSEM/
      mkdir ./outputs/${j}_$i/ARACNe/
   done
done
