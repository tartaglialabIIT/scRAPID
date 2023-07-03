#!/bin/bash

# Cell types
REGTYPE=('RBP_mRNA400')

## Define your python3 path
my_python_path='/mnt/large/jfiorentino/INTERACTomics/Beeline/BEELINE/bin/python'
my_input_folder='/mnt/large/jfiorentino/INTERACTomics/Beeline/inputs/mRNA_lncRNA/'
my_output_folder='${my_main_folder}mRNA_lncRNA/'
my_main_folder="/mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/"

GT='../Ground_Truth_Networks/'

mkdir ${my_main_folder}mRNA_downsampling_results/

for(( k=0; k<100; k+=1 ));
  do
  echo $k


  cp -r ${my_input_folder}HepG2_Smartseq2_RBP_mRNA400/ ${my_input_folder}HepG2_Smartseq2_RBP_mRNA400_$((k+1))/
  cp -r ${my_input_folder}HepG2_DNBelab_RBP_mRNA400/ ${my_input_folder}HepG2_DNBelab_RBP_mRNA400_$((k+1))/
  cp -r ${my_input_folder}K562_CELseq_RBP_mRNA400/ ${my_input_folder}K562_CELseq_RBP_mRNA400_$((k+1))/
  cp -r ${my_input_folder}K562_STORMseq_RBP_mRNA400/ ${my_input_folder}K562_STORMseq_RBP_mRNA400_$((k+1))/
  cp -r ${my_input_folder}K562_Smartseq3_RBP_mRNA400/ ${my_input_folder}K562_Smartseq3_RBP_mRNA400_$((k+1))/


  cp -r ${my_output_folder}HepG2_Smartseq2_RBP_mRNA400/ ${my_output_folder}HepG2_Smartseq2_RBP_mRNA400_$((k+1))/
  cp -r ${my_output_folder}HepG2_DNBelab_RBP_mRNA400/ ${my_output_folder}HepG2_DNBelab_RBP_mRNA400_$((k+1))/
  cp -r ${my_output_folder}K562_CELseq_RBP_mRNA400/ ${my_output_folder}K562_CELseq_RBP_mRNA400_$((k+1))/
  cp -r ${my_output_folder}K562_STORMseq_RBP_mRNA400/ ${my_output_folder}K562_STORMseq_RBP_mRNA400_$((k+1))/
  cp -r ${my_output_folder}K562_Smartseq3_RBP_mRNA400/ ${my_output_folder}K562_Smartseq3_RBP_mRNA400_$((k+1))/

  cp -r ${my_output_folder}HepG2_Smartseq2_RBP_mRNA400_cfilt/ ${my_output_folder}HepG2_Smartseq2_RBP_mRNA400_cfilt_$((k+1))/
  cp -r ${my_output_folder}HepG2_DNBelab_RBP_mRNA400_cfilt/ ${my_output_folder}HepG2_DNBelab_RBP_mRNA400_cfilt_$((k+1))/
  cp -r ${my_output_folder}K562_CELseq_RBP_mRNA400_cfilt/ ${my_output_folder}K562_CELseq_RBP_mRNA400_cfilt_$((k+1))/
  cp -r ${my_output_folder}K562_STORMseq_RBP_mRNA400_cfilt/ ${my_output_folder}K562_STORMseq_RBP_mRNA400_cfilt_$((k+1))/
  cp -r ${my_output_folder}K562_Smartseq3_RBP_mRNA400_cfilt/ ${my_output_folder}K562_Smartseq3_RBP_mRNA400_cfilt_$((k+1))/

  GTFLAGS=('eCLIP_HepG2_single_FC_1_P_3_canonical-network.csv')


  CELLTYPE=('HepG2_Smartseq2' 'HepG2_DNBelab')

  for i in "${CELLTYPE[@]}";
    do
    for j in "${REGTYPE[@]}";
      do
      for l in "${GTFLAGS[@]}";
	    do
	    echo $l
	    
	    $my_python_path EvaluationPipeline_downsampling.py --folder mRNA_lncRNA --dataID "${i}_${j}_$((k+1))" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	    $my_python_path EvaluationPipeline_downsampling.py --folder mRNA_lncRNA --dataID "${i}_${j}_cfilt_$((k+1))" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT	    
	    cd ${my_output_folder}${l/-network.csv/""}${i}_${j}_$((k+1))
	    
	    mv EPrTgt.csv ${i}_${j}_EPrTgt_$((k+1)).csv
	    mv ${i}_${j}_EPrTgt_$((k+1)).csv ${my_main_folder}mRNA_downsampling_results/
	    
	    cd ${my_output_folder}${l/-network.csv/""}${i}_${j}_cfilt_$((k+1))
	   
	    mv EPrTgt.csv ${i}_${j}_EPrTgt_cfilt_$((k+1)).csv
	    mv ${i}_${j}_EPrTgt_cfilt_$((k+1)).csv ${my_main_folder}mRNA_downsampling_results/
	    
	    cd ${my_input_folder}${i}_${j}_$((k+1))/
	    
	    mv ${l} ${l}_${i}_${j}_$((k+1)).txt 
	    mv ${l/-network.csv/""}_stat_tgt.txt ${l/-network.csv/""}_stat_${i}_${j}_$((k+1)).txt 
	    mv ${l/-network.csv/""}_stat_${i}_${j}_$((k+1)).txt ${my_main_folder}mRNA_downsampling_results/
	    
	    cd /mnt/large/jfiorentino/INTERACTomics/Beeline/
	    
	    rm -r ${my_output_folder}${i}_${j}_$((k+1))/
	    rm -r ${my_output_folder}${i}_${j}_cfilt_$((k+1))/
	    rm -r ${my_output_folder}${l/-network.csv/""}${i}_${j}_$((k+1))/
	    rm -r ${my_output_folder}${l/-network.csv/""}${i}_${j}_cfilt_$((k+1))/
	    rm -r ${my_input_folder}${i}_${j}_cfilt_$((k+1))/
	    done
     done
  done

  

  CELLTYPE=('K562_CELseq' 'K562_STORMseq' 'K562_Smartseq3')

  GTFLAGS=('eCLIP_K562_single_FC_1_P_3_canonical-network.csv')

  for i in "${CELLTYPE[@]}";
    do
    for j in "${REGTYPE[@]}";
      do
      for l in "${GTFLAGS[@]}";
	    do
	    echo $l
	    $my_python_path EvaluationPipeline_downsampling.py --folder mRNA_lncRNA --dataID "${i}_${j}_$((k+1))" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT
	    $my_python_path EvaluationPipeline_downsampling.py --folder mRNA_lncRNA --dataID "${i}_${j}_cfilt_$((k+1))" --norm ${i}NormalizedData_${j}.csv --pseudo ${i}PseudoTime.csv --ref $l --gt_folder $GT	    
	    
	    cd ${my_output_folder}${l/-network.csv/""}${i}_${j}_$((k+1))
	    
	    mv EPrTgt.csv ${i}_${j}_EPrTgt_$((k+1)).csv
	    mv ${i}_${j}_EPrTgt_$((k+1)).csv ${my_main_folder}mRNA_downsampling_results/
	    
	    cd ${my_output_folder}${l/-network.csv/""}${i}_${j}_cfilt_$((k+1))
	    
	    mv EPrTgt.csv ${i}_${j}_EPrTgt_cfilt_$((k+1)).csv
	    mv ${i}_${j}_EPrTgt_cfilt_$((k+1)).csv ${my_main_folder}mRNA_downsampling_results/
	    
	    cd ${my_input_folder}${i}_${j}_$((k+1))/
	    
	    mv ${l} ${l}_${i}_${j}_$((k+1)).txt 
	    mv ${l/-network.csv/""}_stat_tgt.txt ${l/-network.csv/""}_stat_${i}_${j}_$((k+1)).txt 
	    mv ${l/-network.csv/""}_stat_${i}_${j}_$((k+1)).txt ${my_main_folder}mRNA_downsampling_results/
	    
	    cd /mnt/large/jfiorentino/INTERACTomics/Beeline/
	    
	    rm -r ${my_output_folder}${i}_${j}_$((k+1))/
	    rm -r ${my_output_folder}${i}_${j}_cfilt_$((k+1))/
	    rm -r ${my_output_folder}${l/-network.csv/""}${i}_${j}_$((k+1))/
	    rm -r ${my_output_folder}${l/-network.csv/""}${i}_${j}_cfilt_$((k+1))/
	    rm -r ${my_input_folder}${i}_${j}_cfilt_$((k+1))/

	    done
     done
  done
done

$my_python_path PostProcessDownsampling.py -d mRNA_downsampling_results -f mRNA_lncRNA

$my_python_path MakeBoxplot_DownSample.py -d mRNA_downsampling_results -i ${my_input_folder} -r ${my_output_folder}
