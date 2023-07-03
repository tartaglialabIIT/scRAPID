import ast
import numpy as np
import pandas as pd
import pickle
import argparse
import os.path

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Run Gene Regulatory Network inference pipeline.')
    
    parser.add_argument('--catflag', default='no',help='Path to the ground truth network file')
    
    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts

opts = parse_arguments()
catflag=opts.catflag


base_folder='/mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/COINTER_RBPs/'

datasets=['HEK293T','HEK293T_smartseq3','HCT116']
strings=['eCLIP_HepG2_single_FC_1_P_3_canonical','eCLIP_HepG2_single_FC_1_P_3_canonical','eCLIP_HepG2_single_FC_1_P_3_canonical']
sets=['RBP_RNA1000','RBP_RNA2000', 'RBP_RNA3000']

for (d,s) in zip(datasets,strings):
	print(d)
	res_folder=base_folder+s+d+'_'
	
	for se in sets:
		print(se)
		if catflag=='no':
			res_folder2=res_folder+se+'/'
			res_folder0=base_folder+d+'_'+se+'/'
		elif catflag=='yes':
			res_folder2=res_folder+se+'_cfilt/'
			res_folder0=base_folder+d+'_'+se+'_cfilt/'
		
		complexes=pd.read_csv(res_folder2+"Complexes.csv",index_col=0)
		
		
		methods=complexes.index
		
		for met in methods:
				newfilename=res_folder2+'Complexes_'+met+'.rnk'
				pp=ast.literal_eval(complexes.loc[met][0])
				print(met)
				if isinstance(pp,int)==False:
					df=pd.DataFrame(pp)
					df.columns= ['RBP1', 'RBP2','tgt1','tgt2','Jaccard']
					df['Edges']=[''.join(sorted(filter(None, x))) for x in df.loc[:,['RBP1','RBP2']].to_numpy()]
					df=df.drop_duplicates(subset=['Edges'])
					
					df.sort_values('Jaccard',ascending=False,inplace=True)
					tmpdf=df.loc[:,['Edges','Jaccard']]
					tmpdf.to_csv(newfilename,sep='\t',index=False)
