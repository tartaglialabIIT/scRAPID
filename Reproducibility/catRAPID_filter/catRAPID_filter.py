import pandas as pd
import numpy as np
import os
import argparse

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Make benchmark heatmap.')

    parser.add_argument('-t','--threshold', default=30,
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
    
def main():
    
	opts = parse_arguments()
	threshold=opts.threshold
    
	base_folder='../Inference/outputs/'
	catRAPID_folder='./catRAPID/'
	
	catRAPID_data=pd.read_csv(catRAPID_folder+'your_catrapid_table.csv',index_col=0)
	catRAPID_data['Edges']=catRAPID_data.protein_name+'|'+catRAPID_data.gene_name
	catRAPID_data=catRAPID_data.set_index('Edges')
	
#	regs='RBPs_1000'
	regs=['RBP_RNA500','RBP_RNA1000']	
	ctypes=[ 'K562_9CL_SCAN','K562_9CL_SCAN','K562_UMI200_SCAN','K562_CELseq','K562_STORMseq','K562_Smartseq3','HepG2_Smartseq2','HepG2_DNBelab','HepG2_9CL_SCAN']
	methods=['PIDC','GRNBOOST2','SINCERITIES', 'TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR']
	
	# Post process ARACNe output
	for ct in ctypes:
		for reg in regs:
			out_folder1=base_folder+ct+'_'+reg+'_cfilt/'
			if os.path.isdir(out_folder1)==False:
				os.mkdir(out_folder1)
			for met in methods:
				met_out_folder1=out_folder1+met+'/'
				if os.path.isdir(met_out_folder1)==False:
					os.mkdir(met_out_folder1)
				method_folder=base_folder+ct+'_'+reg+'/'+met+'/'
				if os.path.isfile(method_folder+'rankedEdges.csv')==True:
					data=pd.read_csv(method_folder+'rankedEdges.csv',delimiter='\t')
					data['Edges']=data['Gene1']+'|'+data['Gene2']
					data=data.set_index('Edges')
					
					# Filter the inferred edges with catRAPID
					#keys =['Gene1','Gene2']
					i1 = catRAPID_data[catRAPID_data['max']>float(threshold)].index
					i2 = data.index
					data1=data[i2.isin(i1)]
					data1=data1.loc[:,['Gene1','Gene2','EdgeWeight']]
					data1.sort_values('EdgeWeight',inplace=True,ascending=False)
					data1.to_csv(met_out_folder1+'rankedEdges.csv',sep='\t')
				
if __name__ == '__main__':
	main()

