import pandas as pd
import numpy as np
import os
import argparse
import csv

def get_delimiter(file_path, bytes = 4096):
	sniffer = csv.Sniffer()
	data = open(file_path, "r").read(bytes)
	delimiter = sniffer.sniff(data).delimiter
	return delimiter


def get_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description='Post-process inferred GRN to obtain BEELINE-like formatted output files')

    parser.add_argument('--dataID', default='MyData',help='Name of the dataset')
    
    parser.add_argument('--folder', default='ExampleFolder',help='Folder containing the input data')

    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts

base_folder=opts.folder
dataID=opts.dataID


#regs=['TF_RNA500','TF_RNA1000','RBP_RNA500', 'RBP_RNA1000']
regs=['RBP_RNA500']

ctypes=['HepG2_Smartseq2','HepG2_DNBelab','HepG2_9CL_SCAN','K562_CELseq','K562_STORMseq','K562_Smartseq3','K562_9CL_SCAN', 'K562_UMI200_SCAN']

# Post process ARACNe and TENET output

if 'RBP' in dataID:
	reg2='RBPs'
elif 'TF' in reg:
	reg2='TFs'
		# Post process TENET inferred edges
TENET_folder=base_folder+dataID+'/TENET/'
df=pd.read_csv(TENET_folder+"TE_result_matrix.txt",delimiter='\t',index_col=0)
new_df=df.stack().reset_index()
new_df.columns=['Gene1','Gene2','EdgeWeight']
new_df.to_csv(TENET_folder+'rankedEdges.csv',index=False)

TENET_A=base_folder+dataID+'/TENET_A/'
TENET_B=base_folder+dataID+'/TENET_B/'
num_lines = sum(1 for line in open(TENET_A+'TE_result_matrix.fdr0.01.trimIndirect-0.1.sif'))
if os.path.getsize(TENET_A+'TE_result_matrix.fdr0.01.trimIndirect-0.1.sif') > 0 and num_lines>1:
	tenet_resA=np.loadtxt(TENET_A+'TE_result_matrix.fdr0.01.trimIndirect-0.1.sif',dtype=str)
	tenet_df=pd.DataFrame(data=tenet_resA,columns=['Gene1', 'EdgeWeight','Gene2'])
	tenet_df.sort_values('EdgeWeight',inplace=True,ascending=False)
	tenet_df=tenet_df[['Gene1','Gene2','EdgeWeight']]
	# Save the results in a csv file
	tenet_df.to_csv(TENET_A+'rankedEdges.csv',index=False,sep='\t')
num_lines = sum(1 for line in open(TENET_B+'TE_result_matrix.fdr0.5.trimIndirect-0.1.sif'))
if os.path.getsize(TENET_B+'TE_result_matrix.fdr0.5.trimIndirect-0.1.sif') > 0 and num_lines>1:
	tenet_resB=np.loadtxt(TENET_B+'TE_result_matrix.fdr0.5.trimIndirect-0.1.sif',dtype=str)
	tenet_df=pd.DataFrame(data=tenet_resB,columns=['Gene1', 'EdgeWeight','Gene2'])
	tenet_df.sort_values('EdgeWeight',inplace=True,ascending=False)
	tenet_df=tenet_df[['Gene1','Gene2','EdgeWeight']]
	# Save the results in a csv file
	tenet_df.to_csv(TENET_B+'rankedEdges.csv',index=False,sep='\t')
ARACNe_folder=base_folder+dataID+'/ARACNe/'
ARACNe_FDR_folder=base_folder+dataID+'/ARACNe_FDR/'

if os.path.isdir(ARACNe_FDR_folder)==False:
	os.mkdir(ARACNe_FDR_folder)
edges=pd.read_csv(ARACNe_folder+'ARACNE-RUN-'+dataID+'_'+reg2+'_4col.tsv',header=None,delimiter='\t')
edges.columns=['Gene1','Gene2','EdgeWeight','p-value']
edges2=edges.copy()
edges2=edges2.loc[:,['Gene1','Gene2','p-value']]
edges2.columns=['Gene1','Gene2','EdgeWeight']
edges2['EdgeWeight']= -np.log10(edges2['EdgeWeight'])
edges2.to_csv(ARACNe_FDR_folder+'rankedEdges.csv',index=False,sep='\t')

methods=['PIDC','GRNBOOST2','SINCERITIES','TENET', 'TENET_A','TENET_B','DeePSEM','ARACNe_FDR']

i=0
for met in methods:
	out_folder=base_folder+dataID+'/'+met+'/'
	if os.path.isfile(out_folder+'rankedEdges.csv')==True:
		delim=get_delimiter(out_folder+'rankedEdges.csv')
		rankedEdges=pd.read_csv(out_folder+'rankedEdges.csv',delimiter=delim)
		rankedEdges.sort_values('EdgeWeight',inplace=True, ascending=False)
		rankedEdges.to_csv(out_folder+'rankedEdges.csv',sep='\t')
	i+=1
