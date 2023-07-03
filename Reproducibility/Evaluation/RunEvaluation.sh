#!/bin/bash

# Cell types

CELLTYPE=('HepG2_Smartseq2' 'HepG2_DNBelab' 'HepG2_9CL_SCAN')

REGTYPE=('RBP_RNA500' 'RBP_RNA1000')

### Define your python3 path
my_python_path='/mnt/large/jfiorentino/INTERACTomics/Beeline/BEELINE/bin/python'

GT='../Ground_Truth_Networks/'

GTFLAGS=('eCLIP_HepG2_single_FC_1_P_3_canonical-network.csv' 'shRNA_seq_indirect_HepG2_strict-network.csv')


for i in "${CELLTYPE[@]}";
  do
  for j in "${REGTYPE[@]}";
    do
    for l in "${GTFLAGS[@]}";
	  do
	  echo $l
	  $my_python_path EvaluationPipeline.py --dataID "${i}_${j}" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	  done
   done
done


REGTYPE=('TF_RNA500' 'TF_RNA1000')
GTFLAGS=('ChIP_seq_strict_HepG2-network.csv')
for i in "${CELLTYPE[@]}";
  do
  for j in "${REGTYPE[@]}";
    do
    for l in "${GTFLAGS[@]}";
	  do
	  echo $l
	  $my_python_path EvaluationPipeline.py --dataID "${i}_${j}" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	  done
   done
done



CELLTYPE=('K562_CELseq' 'K562_STORMseq' 'K562_Smartseq3' 'K562_9CL_SCAN' 'K562_UMI200_SCAN')
REGTYPE=('RBP_RNA500' 'RBP_RNA1000')

GTFLAGS=('eCLIP_K562_single_FC_1_P_3_canonical-network.csv' 'shRNA_seq_indirect_K562_strict-network.csv')

for i in "${CELLTYPE[@]}";
  do
  for j in "${REGTYPE[@]}";
    do
    for l in "${GTFLAGS[@]}";
	  do
	  echo $l
	  $my_python_path EvaluationPipeline.py --dataID "${i}_${j}" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	  done
   done
done


REGTYPE=('TF_RNA500' 'TF_RNA1000')
GTFLAGS=('ChIP_seq_strict_K562-network.csv')
for i in "${CELLTYPE[@]}";
  do
  for j in "${REGTYPE[@]}";
    do
    for l in "${GTFLAGS[@]}";
	  do
	  echo $l
	  $my_python_path EvaluationPipeline.py --dataID "${i}_${j}" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	  done
   done
done
