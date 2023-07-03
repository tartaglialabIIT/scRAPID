# Rome, 28/04/2023. Jonathan Fiorentino

# This python script runs a pipeline for GRN inference from single-cell transcriptomic data.
# It is a wrapper for BEELINE pipeline and it assumes that BEELINE and the associated packages are installed.
 
# It takes as inputs the matrix of normalized counts, the matrix of raw counts 
# and the pseudotime values with the cells ordered as in the count matrices

# Then it runs 4 GRN inference methods through the BEELINE pipeline: PIDC, SINCERITIES and GRNBOOST2 from the BEELINE pipeline, plus TENET. 
# This is achieved using the script BLRunner.py with the configuration file config.yaml.

# The configuration file for running BEELINE is created automatically within the script

# The inputs are the matrix of the normalized counts (genes x cells) and the pseudotime values in csv files

# Additionally, we run TENET, which takes in input the matrix of raw counts (note that it is cells x genes),
# the pseudotime values and a cell selection file. By default all the cells are selected, but a desired 
# percentage can be defined by the user.

# The output is a ranked list of edges for each GRN inference method (a folder is created for each method)

import argparse
import yaml
import numpy as np
import pandas as pd
import os
import subprocess
import time

def get_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description='Run Gene Regulatory Network inference pipeline from single-cell transcriptomic data.')

    parser.add_argument('--dataID', default='MyData',help='Name of the dataset')
        
    parser.add_argument('--norm', default='NormalizedCounts.csv',help='Path to normalized counts file (genes x cells)')
        
    parser.add_argument('--raw', default='RawCounts.csv',help='Path to raw counts file (genes x cells)')
    
    parser.add_argument('--pseudo', default='PseudotimeValues.csv',help='Path to pseudotime values file (cells ordered as in the count files)')
        
    parser.add_argument('--select', default=100,type=int,help='Percentage of cells to be selected (default = all the cells) for TENET')
    
    parser.add_argument('--folder', default='ExampleFolder',help='Folder containing the input data')

    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts

def CreateConfigFile(ID, norm_data_f,pseudotime_f,folder):
	data={'input_settings': {'input_dir': 'inputs','dataset_dir': folder,'datasets': [{'name': ID,'exprData': norm_data_f,'cellData': pseudotime_f,'trueEdges': 'refNetwork.csv'}],'algorithms': [{'name': 'PIDC', 'params': {'should_run': [True]}},{'name': 'GRNBOOST2', 'params': {'should_run': [True]}},{'name': 'PPCOR', 'params': {'should_run': [False], 'pVal': [0.01]}},{'name': 'SINCERITIES', 'params': {'should_run': [True], 'nBins': [10]}}]},'output_settings': {'output_dir': 'outputs', 'output_prefix': ID}}

	with open('./config-files/'+ID+'_config.yaml', 'w') as outfile:
		yaml.dump(data, outfile, default_flow_style=False)

def main():
	opts = parse_arguments()
	datasetID=opts.dataID
	norm_data_file = opts.norm
	raw_data_file = opts.raw
	pseudotime_file=opts.pseudo
	sel_cells=opts.select
	my_folder=opts.folder
	
	print('--------Creating configuration file for BEELINE------------')
	# Create the yaml configuration file for BEELINE
	CreateConfigFile(datasetID,norm_data_file,pseudotime_file,my_folder)
	conf_file='./config-files/'+datasetID+'_config.yaml'
	print('--------Configuration file created------------')
	
	# Run the BEELINE script
	print('--------Running BEELINE---------')
	BEELINE_cmd = 'python BLRunner.py --config '+conf_file
	(out, err) = subprocess.Popen(BEELINE_cmd, stdout=subprocess.PIPE, shell=True).communicate()        

	# Prepare data for TENET
	# Read the raw data and transpose them (TENET needs a cells x genes matrix)
	raw_data=pd.read_csv('./inputs/'+my_folder+'/'+datasetID+'/'+raw_data_file,index_col=0)
	raw_data_t=raw_data.T
	raw_data_t.to_csv('./TENET/'+datasetID+'raw_data_t.csv')
	
	# Save the pseudotime values as txt
	pseudot=pd.read_csv('./inputs/'+my_folder+'/'+datasetID+'/'+pseudotime_file,index_col=0)
	np.savetxt('./TENET/'+datasetID+'pseudotime.txt',np.c_[pseudot.mean(axis=1,skipna=True)],fmt='%f')

	# Create the cell selection file and save it
	if sel_cells==100:
		my_selection=np.ones(np.array(raw_data_t).shape[0])
	else:
		# Randomly select a sel_cells percentage of the cells and generate the matrix of 0 and 1
		my_selection=np.zeros(np.array(raw_data_t).shape[0])
		indices=np.random.choice(np.arange(0,len(my_selection),1), size=int((sel_cells/100.0)*raw_data_t.shape[0]), replace=False)
		my_selection[indices]=1
	
	# Save the file for cell selection
	np.savetxt('./TENET/'+datasetID+'_select_cells.txt',np.c_[my_selection],delimiter='\t',fmt='%d')
	
	# Run TENET
	# File names for TENET
	rawdata_TENET=datasetID+'raw_data_t.csv'
	pseudo_TENET=datasetID+'pseudotime.txt'
	cells_TENET=datasetID+'_select_cells.txt'
	os.chdir('./TENET/')
	TENET='./outputs/'+my_folder+'/'+datasetID+ '/TENET/'
	if os.path.isdir(TENET)==False:
		os.mkdir(TENET)
	print('--------Running TENET inference---------')
	
	TENET_cmd = 'time -v -o ./outputs/'+my_folder+'/'+ datasetID+ '/TENET/time.txt ' + './TENET '+rawdata_TENET+' 10 '+ pseudo_TENET +' '+cells_TENET +' 1'
	(out, err) = subprocess.Popen(TENET_cmd, stdout=subprocess.PIPE, shell=True).communicate()        
	mv_cmd1='mv TE_result_matrix.txt ./outputs/'+my_folder+'/'+ datasetID+ '/TENET/'
	(out, err) = subprocess.Popen(mv_cmd1, stdout=subprocess.PIPE, shell=True).communicate()  
	
	TENET_A='./outputs/'+my_folder+'/'+datasetID+'/'+'TENET_A/'
	TENET_B='./outputs/'+my_folder+'/'+datasetID+'/'+'TENET_B/'
		
	if os.path.isdir(TENET_A)==False:
		os.mkdir(TENET_A)
		
	if os.path.isdir(TENET_B)==False:
		os.mkdir(TENET_B)
	FDR_threshold='0.01'
	TENET_GRN_cmd = 'python makeGRN.py '+FDR_threshold+' '+TENET
	(out, err) = subprocess.Popen(TENET_GRN_cmd, stdout=subprocess.PIPE, shell=True).communicate()
		
	mv_cmd1='mv '+TENET+'TE_result_matrix.fdr'+FDR_threshold+'.sif '+ TENET_A
	(out, err) = subprocess.Popen(mv_cmd1, stdout=subprocess.PIPE, shell=True).communicate() 
		
	print('-----------Trimming indirect edges-----------------')
	TENET_trim_cmd='python trim_indirect.py TE_result_matrix.fdr'+FDR_threshold+'.sif -0.1 '+TENET_A
	(out,err) = subprocess.Popen(TENET_trim_cmd, stdout=subprocess.PIPE, shell=True).communicate()     
		
	FDR_threshold='0.5'
	TENET_GRN_cmd = 'python makeGRN.py '+FDR_threshold+' '+TENET
	(out, err) = subprocess.Popen(TENET_GRN_cmd, stdout=subprocess.PIPE, shell=True).communicate()
		
	mv_cmd1='mv '+TENET+'TE_result_matrix.fdr'+FDR_threshold+'.sif ' +TENET_B
	(out, err) = subprocess.Popen(mv_cmd1, stdout=subprocess.PIPE, shell=True).communicate() 
        
	print('-----------Trimming indirect edges-----------------')
	TENET_trim_cmd='python trim_indirect.py TE_result_matrix.fdr'+FDR_threshold+'.sif -0.1 '+TENET_B
	(out,err) = subprocess.Popen(TENET_trim_cmd, stdout=subprocess.PIPE, shell=True).communicate()

if __name__ == '__main__':
	main()
