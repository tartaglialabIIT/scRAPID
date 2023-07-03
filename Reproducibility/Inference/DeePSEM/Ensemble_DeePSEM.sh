#!/bin/bash


# Example script for running DeePSEM on all the scRNA-seq datasets and different gene sets
# We run DeePSEM 10 times (Ensemble DeePSEM) and average the edge weights as recommended by the authors 
# We assume that the working directory contains the subfolders input and output 
END=10

#GENES=('TF_RNA500' 'TF_RNA1000' 'RBP_RNA500' 'RBP_RNA1000')
GENES=('RBP_RNA500')
DATASETS=('HepG2_Smartseq2' 'HepG2_DNBelab' 'HepG2_9CL_SCAN' 'K562_CELseq' 'K562_STORMseq' 'K562_Smartseq3' 'K562_9CL_SCAN' 'K562_UMI200_SCAN')

for i in "${GENES[@]}";
	do
	for k in "${DATASETS[@]}";
	do
	python transpose_input.py --folder ./input/${k}Input/ --filename ${k}NormalizedData_$i.csv	
	
	mkdir -p /data01/jonathan/DeePSEM/output/DeePSEM_${k}_$i	

	for ((j=1;j<=END;j++));
		do 
		echo $j
        /usr/bin/time -v -o ./DeePSEM/output/DeePSEM_${k}_$i/time_$j.txt python main.py --task celltype_GRN --data_file ./input/${k}Input/${k}NormalizedData_$i.csv --setting test --alpha 0.1 --beta 0.01 --n_epoch 150 --save_name ./output/DeePSEM_${k}_$i --out_file_name /res_$j.tsv;		
		
		done
	python DeePSEM_results.py --out_folder ./output/DeePSEM_${k}_$i
		
	done
done
