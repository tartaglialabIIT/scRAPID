# Rome, 28/04/2023. Jonathan Fiorentino

# This python script runs the evaluation pipeline on the GRN inferred from single-cell transcriptomic data.
# It is a wrapper for the functions provided in BEELINE, but I also added some additional evaluation metrics (such as the topkranking)


import argparse
import yaml
import numpy as np
import pandas as pd
import os
import subprocess
import time

def get_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description='Run Gene Regulatory Network evaluation pipeline.')
    
    parser.add_argument('--dataID', default='MyData',help='Name of the dataset')
        
    parser.add_argument('--norm', default='NormalizedCounts.csv',help='Path to normalized counts file (genes x cells)')
            
    parser.add_argument('--pseudo', default='PseudotimeValues.csv',help='Path to pseudotime values file (cells ordered as in the count files)')
    
    parser.add_argument('--ref', default='refNetwork.csv',help='file name for the ground truth network')
    
    parser.add_argument('--gt_folder', default='./GT_folder/',help='Path to the ground truth network file')
    
    parser.add_argument('--folder', default='ExampleFolder',help='Folder containing the input data')
    
    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts

def CreateConfigFile(ID, norm_data_f,pseudotime_f, ref_net,folder):
	data={'input_settings': {'input_dir': 'inputs','dataset_dir': folder,'datasets': [{'name': ID,'exprData': norm_data_f,'cellData': pseudotime_f,'trueEdges': ref_net}],'algorithms': [{'name': 'PIDC', 'params': {'should_run': [True]}},{'name': 'GRNBOOST2', 'params': {'should_run': [True]}},{'name': 'SINCERITIES', 'params': {'should_run': [True], 'nBins': [10]}},{'name': 'TENET', 'params': {'should_run': [True]}},{'name': 'TENET_A', 'params': {'should_run': [True]}},{'name': 'TENET_B', 'params': {'should_run': [True]}},{'name': 'DeePSEM_CT', 'params': {'should_run': [True]}},{'name': 'ARACNe_FDR', 'params': {'should_run': [True]}}]},'output_settings': {'output_dir': 'outputs', 'output_prefix': ref_net.replace('-network.csv','')+ID}}

	with open('../Inference/config-files/'+ref_net.replace('-network.csv','')+ID+'_config.yaml', 'w') as outfile:
		yaml.dump(data, outfile, default_flow_style=False)

def FilterGroundTruth(GT_folder,refnet_file,norm_data_file, dataset_ID):
	netDF = pd.read_csv(GT_folder+refnet_file)
	if 'Gene1' not in netDF.columns:
		netDF.rename(columns={'reg': 'Gene1', 'target': 'Gene2'}, inplace=True)
	expr_df = pd.read_csv('../Inference/inputs/'+dataset_ID+'/'+norm_data_file, header=0, index_col=0)
	netDF = netDF[(netDF.Gene1.isin(expr_df.index)) & (netDF.Gene2.isin(expr_df.index))]
	# Remove self-loops.
	netDF = netDF[netDF.Gene1 != netDF.Gene2]
	netDF.drop_duplicates(keep = 'first', inplace=True)
	netDF.to_csv('../Inference/inputs/'+dataset_ID+'/'+refnet_file, index=False)
	allNodes = set(netDF.Gene1.unique()).union(set(netDF.Gene2.unique()))
	nTFs = expr_df[expr_df.index.isin(netDF.Gene1.unique())].shape[0]
	nGenes = expr_df[expr_df.index.isin(allNodes)].shape[0]
	refnet_string=refnet_file.replace('-network.csv','')
	np.savetxt('../Inference/inputs/'+dataset_ID+'/'+refnet_string+'_stat.txt',np.c_[[nTFs,nGenes,netDF.shape[0],netDF.shape[0]/((nTFs*nGenes)-nTFs)]],delimiter='\t')
	
	# Save the statistics considering only reg-reg interactions
	netDF1=netDF.loc[netDF.Gene2.isin(set(netDF.Gene1))]
	allNodes1 = set(netDF1.Gene1.unique()).union(set(netDF1.Gene2.unique()))
	nTFs1 = expr_df[expr_df.index.isin(netDF1.Gene1.unique())].shape[0]
	nGenes1 = expr_df[expr_df.index.isin(allNodes1)].shape[0]
	np.savetxt('../Inference/inputs/'+dataset_ID+'/'+refnet_string+'_stat_reg.txt',np.c_[[nTFs1,nGenes1,netDF1.shape[0],netDF1.shape[0]/((nTFs1*nGenes1)-nTFs1)]],delimiter='\t')

	# Save the statistics considering only reg-RNA interactions (no interactions involving mRNAs coding for regulators (TFs or RBPs))
	RBPs_1=np.loadtxt('./RBPs_HepG2.txt',dtype=str)
	RBPs_2=np.loadtxt('./RBPs_K562.txt',dtype=str)
	RBPs=list(set(list(RBPs_1)+list(RBPs_2)))
	netDF2=netDF.loc[~netDF.Gene2.isin(set(RBPs))]
	allNodes2 = set(netDF2.Gene1.unique()).union(set(netDF2.Gene2.unique()))
	nTFs2 = expr_df[expr_df.index.isin(netDF2.Gene1.unique())].shape[0]
	nGenes2 = expr_df[expr_df.index.isin(set(netDF2.Gene2.unique())-set(RBPs))].shape[0]
	np.savetxt('../Inference/inputs/'+dataset_ID+'/'+refnet_string+'_stat_tgt.txt',np.c_[[nTFs2,nGenes2,netDF2.shape[0],netDF2.shape[0]/(nTFs2*nGenes2)]],delimiter='\t')
	
def main():
	opts = parse_arguments()
	datasetID=opts.dataID
	norm_data_file = opts.norm
	pseudotime_file=opts.pseudo
	refnet_file=opts.ref
	GT_folder=opts.gt_folder
	my_folder=opts.folder
	
	print('Intersection of ground truth with scRNA-seq data')
	# The results after catRAPID filter have a "cfilt" string at the end of the dataset ID
	if 'cfilt' not in datasetID:
		FilterGroundTruth(GT_folder,refnet_file,norm_data_file,datasetID)
	else:
		head, sep, tail = datasetID.partition('cfilt')
		head=head[:-1]
		tail=tail[2:]
		ID1=head+tail
		ID2=datasetID
		INPUT_cmd = 'cp -r ../Inference/inputs/'+ID1+' ../Inference/inputs/'+ID2
		(out, err) = subprocess.Popen(INPUT_cmd, stdout=subprocess.PIPE, shell=True).communicate()
	print('--------Creating configuration file for evaluation------------')
	# Create the yaml configuration file for BEELINE (now it includes TENET, TENET_A, TENET_B, DeePSEM and ARACNe)
	CreateConfigFile(datasetID,norm_data_file,pseudotime_file,refnet_file,my_folder)
	conf_file='../Inference/config-files/'+refnet_file.replace('-network.csv','')+datasetID+'_config.yaml'
	
	print('--------Configuration file created------------')
	
	print('--------Running BEELINE evaluation---------')
	# Here we compute the EPR, the percentage of true positives vs ranking and hubs, other options are available
	# e.g. --auc --motifs --paths --borda --experimental_topology --hubs --eprtgt --topk_ranking_tgt --hubstgt --cointer

	BEELINE_cmd = 'python BLEvaluator_adapted.py --config '+conf_file +' --epr --topk_ranking --hubs --cointer'
	(out, err) = subprocess.Popen(BEELINE_cmd, stdout=subprocess.PIPE, shell=True).communicate() 	

if __name__ == '__main__':
	main()
