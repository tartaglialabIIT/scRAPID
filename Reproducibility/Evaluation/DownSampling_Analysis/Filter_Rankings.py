import pandas as pd
import numpy as np
import os
import argparse
import csv

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Make benchmark heatmap.')
        
    parser.add_argument('-f','--folder', default='',
        help="Benchmark measure to be used")
        
    parser.add_argument('-d','--dataset', default='',
        help="Benchmark measure to be used")
    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def get_delimiter(file_path, bytes = 4096):
	sniffer = csv.Sniffer()
	data = open(file_path, "r").read(bytes)
	delimiter = sniffer.sniff(data).delimiter
	return delimiter

opts = parse_arguments()
    
    
ct=opts.dataset
myfolder=opts.folder

base_folder='../Inference/outputs/'+myfolder+'/'
gt_folder='../Ground_Truth_Networks/'
input_folder='../Inference/inputs/'+myfolder+'/'+ct+'/'
methods=['PIDC','GRNBOOST2','SINCERITIES','TENET', 'TENET_A','TENET_B','DeePSEM','ARACNe_FDR']

# Load the eCLIP ground truth
if 'HepG2' in ct:
	FullNet=pd.read_csv(gt_folder+'eCLIP_HepG2_single_FC_1_P_3_canonical-network.csv')
	DownNet=pd.read_csv(input_folder+'eCLIP_HepG2_single_FC_1_P_3_canonical-network.csv')
elif 'K562' in ct:
	FullNet=pd.read_csv(gt_folder+'eCLIP_K562_single_FC_1_P_3_canonical-network.csv')
	DownNet=pd.read_csv(input_folder+'eCLIP_K562_single_FC_1_P_3_canonical-network.csv')

FullNet['Edges']=FullNet.reg+'|'+FullNet.target
DownNet['Edges']=DownNet.Gene1+'|'+DownNet.Gene2

# Find intersection between FullNet and DownNet edges
intersEdges=list(set(FullNet['Edges']).intersection(set(DownNet['Edges'])))

# Remove those edges from FullNet
FullNet=FullNet[~FullNet.Edges.isin(intersEdges)]


for met in methods:
	out_folder=base_folder+ct+'/'+met+'/'
	if os.path.isfile(out_folder+'rankedEdges.csv')==True:
		delim=get_delimiter(out_folder+'rankedEdges.csv')
		rankedEdges=pd.read_csv(out_folder+'rankedEdges.csv',delimiter=delim)
		rankedEdges.sort_values('EdgeWeight',inplace=True, ascending=False)
		rankedEdges['Edges']=rankedEdges.Gene1+'|'+rankedEdges.Gene2
		
		# Intersection between 
		intersEdges2=list(set(FullNet['Edges']).intersection(set(rankedEdges['Edges'])))
		rankedEdges=rankedEdges[~rankedEdges.Edges.isin(intersEdges2)]
		
		rankedEdges.loc[:,['Gene1','Gene2','EdgeWeight']].to_csv(out_folder+'rankedEdges.csv',sep='\t')
