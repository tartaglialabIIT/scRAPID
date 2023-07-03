# Rome, 16/11/2021. Jonathan Fiorentino


import argparse
import yaml
import numpy as np
import pandas as pd
import os
import subprocess
import time
import random

def get_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description='Run Gene Regulatory Network inference pipeline.')
    
    parser.add_argument('--folder', default='MyData',help='Name of the dataset')
    
    parser.add_argument('--dataID', default='MyData',help='Name of the dataset')
        
    parser.add_argument('--norm', default='NormalizedCounts.csv',help='Path to normalized counts file (genes x cells)')
            
    parser.add_argument('--pseudo', default='PseudotimeValues.csv',help='Path to pseudotime values file (cells ordered as in the count files)')
    
    parser.add_argument('--ref', default='refNetwork.csv',help='Path to the ground truth network file')
    
    parser.add_argument('--gt_folder', default='./GT_folder/',help='Path to the ground truth network file')
    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts

def CreateConfigFile(ID, norm_data_f,pseudotime_f, ref_net, myfolder):
	#if 'cfilt' in ID:
	#	head, sep, tail = ID.partition('cfilt')
	#	head=head[:-1]
	#	tail=tail[2:]
	#	ID1=head+tail
	#	ID2=ID
	#else:
	#	ID1=ID
	#	ID2=ID
	#print(ID1,ID2)
	data={'input_settings': {'input_dir': 'inputs','dataset_dir': myfolder,'datasets': [{'name': ID,'exprData': norm_data_f,'cellData': pseudotime_f,'trueEdges': ref_net}],'algorithms': [{'name': 'PIDC', 'params': {'should_run': [True]}},{'name': 'GRNBOOST2', 'params': {'should_run': [True]}},{'name': 'SINCERITIES', 'params': {'should_run': [True], 'nBins': [10]}},{'name': 'TENET', 'params': {'should_run': [True]}},{'name': 'TENET_A', 'params': {'should_run': [True]}},{'name': 'TENET_B', 'params': {'should_run': [True]}},{'name': 'DeePSEM_CT', 'params': {'should_run': [True]}},{'name': 'DeePSEM_NCT', 'params': {'should_run': [True]}},{'name': 'ARACNe_MI', 'params': {'should_run': [True]}},{'name': 'ARACNe_FDR', 'params': {'should_run': [True]}}]},'output_settings': {'output_dir': 'outputs', 'output_prefix': ref_net.replace('-network.csv','')+ID}}

	with open('../Inference/config-files/'+ref_net.replace('-network.csv','')+ID+'_config.yaml', 'w') as outfile:
		yaml.dump(data, outfile, default_flow_style=False)

def FilterGroundTruth(GT_folder,refnet_file,norm_data_file, dataset_ID, myfolder):
	origIDlist=dataset_ID.split('_')[:-1]
	origID='_'.join(origIDlist)
	netDF = pd.read_csv(GT_folder+refnet_file)
	if 'Gene1' not in netDF.columns:
		netDF.rename(columns={'reg': 'Gene1', 'target': 'Gene2'}, inplace=True)
	expr_df = pd.read_csv('../Inference/inputs/'+myfolder+'/'+dataset_ID+'/'+norm_data_file, header=0, index_col=0)
	netDF = netDF[(netDF.Gene1.isin(expr_df.index)) & (netDF.Gene2.isin(expr_df.index))]
	# Remove self-loops.
	netDF = netDF[netDF.Gene1 != netDF.Gene2]
	# Remove duplicates (there are some repeated lines in the ground-truth networks!!!). 
	netDF.drop_duplicates(keep = 'first', inplace=True)
	
	# Read the corresponding statistic for the lncRNA dataset
	dataset_ID2= origID.replace('mRNA','lnc')
	net_lnc = pd.read_csv('../Inference/inputs/'+myfolder+'/'+dataset_ID2+'/'+refnet_file)
	RBPs_1=np.loadtxt('../Evaluation/RBPs_HepG2.txt',dtype=str)
	RBPs_2=np.loadtxt('../Evaluation/RBPs_K562.txt',dtype=str)
	RBPs=list(set(list(RBPs_1)+list(RBPs_2)))
	net_lnc=net_lnc.loc[~net_lnc.Gene2.isin(set(RBPs))]
	netDF=netDF.loc[~netDF.Gene2.isin(set(RBPs))]
	rbp_mRNA=list(netDF.groupby('Gene1').count().sort_values('Gene2',ascending=False).index)
	
	nedges_lncRNA=len(net_lnc)

	if len(set(net_lnc.Gene1))<len(set(netDF.Gene1)):
	    ntgt_lncRNA=len(set(net_lnc.Gene2))
	    nRBP_lncRNA=len(set(net_lnc.Gene1))
	    sampled_RBPs=random.sample(list(set(netDF.Gene1)),nRBP_lncRNA)
	else: 
	    ntgt_lncRNA=int((len(set(net_lnc.Gene2))*len(set(net_lnc.Gene1)))/len(set(netDF.Gene1)))
	    nRBP_lncRNA=len(set(netDF.Gene1))
	    sampled_RBPs=list(set(netDF.Gene1))


	tmpNetDf=netDF.copy()
	tmpNetDf['Edges']=tmpNetDf.Gene1+'|'+tmpNetDf.Gene2
	lsampled_edges=[]
	lsampled_RBPs=[]
	lsampled_tgts=[]
	nsampled_RBPs=0
	nsampled_tgts=0
	    
	flag=0
	    
	while flag==0:
	    # Sample edges recursively
	    # start sampling one edge per rbp
	    for rbp in sampled_RBPs:
	        if len(tmpNetDf[tmpNetDf.Gene1==rbp])>0:
	            if nsampled_tgts<ntgt_lncRNA:
	                lsampled_tgt=random.sample(list(tmpNetDf[tmpNetDf.Gene1==rbp]['Gene2']),1)
	                lsampled_edges.append(rbp+'|'+lsampled_tgt[0])
	                lsampled_RBPs.append(rbp)
	                lsampled_tgts.append(lsampled_tgt[0])
	                lsampled_edges=list(set(lsampled_edges))
	                lsampled_RBPs=list(set(lsampled_RBPs))
	                lsampled_tgts=list(set(lsampled_tgts))
	            else:
	                if len(set(lsampled_tgts).intersection(set(list(tmpNetDf[tmpNetDf.Gene1==rbp]['Gene2']))))>0:
	                    possible_tgts=list(set(lsampled_tgts).intersection(set(list(tmpNetDf[tmpNetDf.Gene1==rbp]['Gene2']))))
	                    lsampled_tgt=random.sample(possible_tgts,1)
	                    lsampled_edges.append(rbp+'|'+lsampled_tgt[0])
	                    lsampled_RBPs.append(rbp)
	                    lsampled_tgts.append(lsampled_tgt[0])
	                    lsampled_edges=list(set(lsampled_edges))
	                    lsampled_RBPs=list(set(lsampled_RBPs))
	                    lsampled_tgts=list(set(lsampled_tgts))                    
	                
	            nsampled_RBPs=len(lsampled_RBPs)
	            nsampled_tgts=len(lsampled_tgts)
	            
	            tmpNetDf=tmpNetDf[tmpNetDf.Edges!=rbp+'|'+lsampled_tgt[0]]
	            
	            if (nsampled_RBPs==nRBP_lncRNA) & (nsampled_tgts==ntgt_lncRNA) & (len(lsampled_edges)==nedges_lncRNA):
	                flag=1
	                break
	            if len(lsampled_edges)>nedges_lncRNA:
	                flag=1
	                break
	                
	netDF_mRNA = pd.DataFrame(data={'Gene1':[i.split('|')[0] for i in lsampled_edges],'Gene2':[i.split('|')[1] for i in lsampled_edges]})    
	
	netDF_mRNA.to_csv('../Inference/inputs/'+myfolder+'/'+dataset_ID+'/'+refnet_file, index=False)
	allNodes = set(netDF_mRNA.Gene1.unique()).union(set(netDF_mRNA.Gene2.unique()))
	nTFs = len(set([i.split('|')[0] for i in lsampled_edges]))
	nGenes = len(set([i.split('|')[1] for i in lsampled_edges]))
	refnet_string=refnet_file.replace('-network.csv','')
	np.savetxt('../Inference/inputs/'+myfolder+'/'+dataset_ID+'/'+refnet_string+'_stat.txt',np.c_[[nTFs,nGenes,netDF_mRNA.shape[0],netDF_mRNA.shape[0]/((nTFs*nGenes)-nTFs)]],delimiter='\t')
	

	# Save the statistics considering only reg-RNA interactions (no interactions involving mRNAs coding for regulators (TFs or RBPs))
	
	netDF2=netDF.loc[~netDF.Gene2.isin(set(RBPs))]
	allNodes2 = set(netDF2.Gene1.unique()).union(set(netDF2.Gene2.unique()))
	nTFs2 = len(set([i.split('|')[0] for i in lsampled_edges]))
	nGenes2 = len(set([i.split('|')[1] for i in lsampled_edges]))
	np.savetxt('../Inference/inputs/'+myfolder+'/'+dataset_ID+'/'+refnet_string+'_stat_tgt.txt',np.c_[[nTFs2,nGenes2,netDF_mRNA.shape[0],netDF_mRNA.shape[0]/(nTFs2*nGenes2)]],delimiter='\t')
	
def main():
	opts = parse_arguments()
	datasetID=opts.dataID
	norm_data_file = opts.norm
	pseudotime_file=opts.pseudo
	refnet_file=opts.ref
	GT_folder=opts.gt_folder
	
	print('Intersection of ground truth with scRNA-seq data')
	if 'cfilt' not in datasetID:
		FilterGroundTruth(GT_folder,refnet_file,norm_data_file,datasetID,opts.folder)
	else:
		head, sep, tail = datasetID.partition('cfilt')
		#head=head[:-1]
		#tail=tail[2:]
		ID1=head[:-1]+tail
		ID2=datasetID
		INPUT_cmd = 'cp -r ../Inference/inputs/'+opts.folder+'/'+ID1+' ../Inference/inputs/'+opts.folder+'/'+ID2
		(out, err) = subprocess.Popen(INPUT_cmd, stdout=subprocess.PIPE, shell=True).communicate()
	print('--------Creating configuration file for evaluation------------')
	# Create the yaml configuration file for BEELINE
	CreateConfigFile(datasetID,norm_data_file,pseudotime_file,refnet_file,opts.folder)
	conf_file='../Inference/config-files/'+refnet_file.replace('-network.csv','')+datasetID+'_config.yaml'
	print('--------Configuration file created------------')
	
	# Run the BEELINE script (interface with bash here!)
	print('--------Running BEELINE evaluation---------')

	# Remove eCLIP interactions that are not in the downsampled ground truth from the rankings
	FILTER_cmd = '/mnt/large/jfiorentino/INTERACTomics/Beeline/BEELINE/bin/python Filter_Rankings.py --dataset ' + datasetID + ' --folder '+opts.folder
	(out, err) = subprocess.Popen(FILTER_cmd, stdout=subprocess.PIPE, shell=True).communicate() 	
	
	BEELINE_cmd = '/mnt/large/jfiorentino/INTERACTomics/Beeline/BEELINE/bin/python ../Evaluation/BLEvaluator_adapted.py --config '+conf_file +' --eprtgt'
	(out, err) = subprocess.Popen(BEELINE_cmd, stdout=subprocess.PIPE, shell=True).communicate() 	

if __name__ == '__main__':
	main()
