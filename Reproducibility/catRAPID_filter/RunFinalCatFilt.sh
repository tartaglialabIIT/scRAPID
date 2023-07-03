#!/bin/bash

# Cell types
REGTYPE=('RBP_RNA500' 'RBP_RNA1000')

## Define your python3 path
my_python_path='/mnt/large/jfiorentino/INTERACTomics/Beeline/BEELINE/bin/python'

GT='../Ground_Truth_Networks/'

for(( k=0; k<31; k+=5 ));
  do
  echo $k
  $my_python_path catRAPID_filter.py --threshold $k

  GTFLAGS=('eCLIP_HepG2_single_FC_1_P_3_canonical-network.csv')
  ##GTFLAGS=('shRNA_seq_indirect_HepG2_strict-network.csv')


  CELLTYPE=('HepG2_Smartseq2' 'HepG2_DNBelab' 'HepG2_9CL_SCAN')

  for i in "${CELLTYPE[@]}";
    do
    for j in "${REGTYPE[@]}";
      do
      for l in "${GTFLAGS[@]}";
	    do
	    echo $l
	    $my_python_path ../EvaluationPipeline.py --dataID "${i}_${j}_cfilt" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	    #cd /mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/MAIN/${l/-network.csv/""}${i}_${j}_cfilt
	    #mv EPr.csv ${i}_${j}_EPr_${k}.csv
	    #mv EPrTgt.csv ${i}_${j}_EPrTgt_${k}.csv
	    #mv ${i}_${j}_EPr_${k}.csv /mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/catRAPID_filter_results/
	    #mv ${i}_${j}_EPrTgt_${k}.csv /mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/catRAPID_filter_results/
	    #cd /mnt/large/jfiorentino/INTERACTomics/Beeline/
	    done
     done
  done

  CELLTYPE=('K562_CELseq' 'K562_STORMseq' 'K562_Smartseq3' 'K562_9CL_SCAN' 'K562_UMI200_SCAN')

  GTFLAGS=('eCLIP_K562_single_FC_1_P_3_canonical-network.csv')

  for i in "${CELLTYPE[@]}";
    do
    for j in "${REGTYPE[@]}";
      do
      for l in "${GTFLAGS[@]}";
	    do
	    echo $l
	    $my_python_path ../EvaluationPipeline.py --dataID "${i}_${j}_cfilt" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	    #cd /mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/MAIN/${l/-network.csv/""}${i}_${j}_cfilt
	    #mv EPr.csv ${i}_${j}_EPr_${k}.csv
	    #mv EPrTgt.csv ${i}_${j}_EPrTgt_${k}.csv
	    #mv ${i}_${j}_EPr_${k}.csv /mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/catRAPID_filter_results/
	    #mv ${i}_${j}_EPrTgt_${k}.csv /mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/catRAPID_filter_results/
	    #cd /mnt/large/jfiorentino/INTERACTomics/Beeline/
	    done
     done
  done
done

$my_python_path postProcess_catRAPID_filter.py
