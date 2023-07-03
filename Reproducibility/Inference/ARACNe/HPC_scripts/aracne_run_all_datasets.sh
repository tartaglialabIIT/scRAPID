#!/bin/bash

cd MAIN/

DATASETS=('HepG2_Smartseq2' 'HepG2_DNBelab' 'HepG2_9CL_SCAN' 'K562_CELseq' 'K562_STORMseq' 'K562_Smartseq3' 'K562_UMI200_SCAN')

GENES1=('RBP_RNA500')# 'RBP_RNA1000')
#GENES2=('TF_RNA500' 'TF_RNA1000')

for i in "${DATASETS[@]}";
	do
	for j in "${GENES1[@]}";
	    do
	    bash ./aracne-script.sh -n ARACNE-RUN-${i}_${j} -b /work/jfiorentino/ARACNe/MAIN/network_${i}_${j}/ -i /work/jfiorentino/ARACNe/MAIN/Processed_scdata/${i}_${j}processed_${i}_${j}.rds -r /work/jfiorentino/ARACNe/human_hugo/RBPs_human.txt
	    
	    done
	#for k in "${GENES2[@]}";
	    #do
	    #bash ./aracne-script.sh -n ARACNE-RUN-${i}_${k} -b /work/jfiorentino/ARACNe/MAIN/network_${i}_${k}/ -i /work/jfiorentino/ARACNe/MAIN/Processed_scdata/${i}_${j}processed_${i}_${k}.rds -r /work/jfiorentino/ARACNe/human_hugo/TFs_ALL.txt
	    
	    #done
	
done
