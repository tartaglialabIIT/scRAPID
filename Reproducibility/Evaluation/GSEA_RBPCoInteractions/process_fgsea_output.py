import numpy as np
import pandas as pd
import os
import argparse
#import subprocess

def get_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(description='Run Gene Regulatory Network inference pipeline.')
    
    parser.add_argument('--catflag', default='no',help='Path to the ground truth network file')
    return parser

def parse_arguments():
	parser = get_parser()
	opts = parser.parse_args()
	
	return opts


opts = parse_arguments()
base_folder='./'

datasets=['HEK293T','HEK293T_smartseq3','HCT116']
strings=['eCLIP_HepG2_single_FC_1_P_3_canonical','eCLIP_HepG2_single_FC_1_P_3_canonical','eCLIP_HepG2_single_FC_1_P_3_canonical']
sets=['RBP_RNA1000','RBP_RNA2000', 'RBP_RNA3000']

methods =['PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR']

data_list=[]

catflag=opts.catflag

#if catflag=='no':
#	plots_folder='/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_RBPs/fgsea_plots/'
#elif catflag=='yes':
#	plots_folder='/mnt/large/jfiorentino/INTERACTomics/Beeline/fgsea_RBPs/fgsea_plots_cfilt/'

#if os.path.isdir(plots_folder)==False:
#	os.mkdir(plots_folder)

for (d,s) in zip(datasets,strings):
	print(d)
	res_folder=base_folder+s+d+'_'
	for se in sets:
		print(se)
		if catflag=='yes':
			res_folder2=res_folder+se+'_cfilt/'
		elif catflag=='no':
			res_folder2=res_folder+se+'/'
		for met in methods:
			filename= res_folder2+'fgsea_res'+met+'.tsv'
			print(filename)	
			if os.path.isfile(filename):
				data=pd.read_csv(filename,delimiter='\t')
				data['dataset']=d
				data['reg']=se
				data['method']=met
				data_list.append(data)
				#fname= d+'_'+se+'_'+met+'.pdf'
				#mv_cmd='mv '+fname+' '+plots_folder
     				#BEELINE_cmd = '/mnt/large/jfiorentino/INTERACTomics/Beeline/BEELINE/bin/python /mnt/large/jfiorentino/INTERACTomics/Beeline/BLEvaluator_adapted.py --config '+conf_file +' --complex'
				#(out, err) = subprocess.Popen(mv_cmd, stdout=subprocess.PIPE, shell=True).communicate() 	 
			
			
final_df=pd.concat(data_list)
if catflag=='no':
	final_df.to_csv('fgsea_final_results.csv')
elif catflag=='yes':
        final_df.to_csv('fgsea_final_results_cfilt.csv')

