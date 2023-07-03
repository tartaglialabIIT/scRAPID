import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import dataframe_image as dfi
import os
import argparse
from sklearn import preprocessing
import rpy2.robjects as robjects
import pickle

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Make benchmark heatmap.')

    parser.add_argument('-m','--measure', default='EPr',
        help="Benchmark measure to be used")
        
    parser.add_argument('-r','--res_folder', default='/mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/MAIN/',
        help="Folder with the results")
        
    parser.add_argument('-i','--input_folder', default='/mnt/large/jfiorentino/INTERACTomics/Beeline/inputs/MAIN/',
        help="Folder with the results")
    
    parser.add_argument('-g','--gtlab', default='eCLIP',
        help="Folder with the results")
#    parser.add_argument('-o','--output_folder', default='/mnt/large/jfiorentino/INTERACTomics/Beeline/inputs/example/heatmaps/',
#        help="Folder with the results")
        
#    parser.add_argument('-c','--catrapid', default='False',
#        help="Use catRAPID filtered results")
    
    #parser.add_argument('-t','--catrapid_threshold', default='0',
    #    help="Threshold on the catRAPID interaction propensity for filtering (30, 50, 70 available)")
    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def PlotRankings(res_folder,input_folder,celltype,set_lab,ground_truth_labels,methods,fsize,reg,output_folder,at,meas):
	colors=['C'+str(i) for i in range(len(methods))]
	for ct in celltype:
		ctn=ct.split('_')[0]
		i=0
		print(i,len(ground_truth_labels))
		fig,ax=plt.subplots(1,2,figsize=fsize,sharey=True)
		fig.suptitle(ct,fontsize=16)
		for gtlab in ground_truth_labels:
			old_gt_lab=gtlab[:]
			gtlab='eCLIP_'+ctn+'_single_FC_1_P_3_'+at
			my_foldercat=res_folder+gtlab+ct+'_'+reg+'_cfilt/'
			my_folder=res_folder+gtlab+ct+'_'+reg+'/'
			if meas=='Ranking':
				net_fts=np.loadtxt(input_folder+ct+'_'+reg+'/'+gtlab+'_stat.txt',dtype=str)
			elif meas=='RankingTgt':
				net_fts=np.loadtxt(input_folder+ct+'_'+reg+'/'+gtlab+'_stat_tgt.txt',dtype=str)
			ax[0].set_xscale('log')
			#ax[0].set_title(gtlab,fontsize=14)
			ax[0].set_xlabel('Ranking',fontsize=14)
			ax[0].set_ylabel('% true positives',fontsize=14)
			ax[1].set_xscale('log')
			#ax[1].set_title(gtlab + ' catRAPID filter',fontsize=12)
			ax[1].set_xlabel('Ranking',fontsize=14)
			ax[1].set_ylabel('% true positives',fontsize=14)
			# Import the data with the rankings for all the algorithms
			with open(my_folder+meas+'.pickle', 'rb') as handle:
				fulldata = pickle.load(handle)
			with open(my_foldercat+meas+'.pickle', 'rb') as handle:
				fulldata2 = pickle.load(handle)
			j=0
			for met in methods:
				
				data=fulldata.loc[met,:][0]
				data2=fulldata2.loc[met,:][0]
				#print(gtlab,met)
				#print(type(data),data.shape)
				if isinstance(data, np.ndarray) and data.shape[0]>4:
					print(gtlab,met)
					width=data[1,0]-data[0,0]
					xy = np.column_stack(( data[:,0] + width , 100.*(data[:,1]/width) ))
					x_range = np.linspace(min(xy[:,0]),max(xy[:,0]), num=100, endpoint=True)
					r_y = robjects.FloatVector(xy[:,1])
					r_x = robjects.FloatVector(xy[:,0])
		
					r_smooth_spline = robjects.r['smooth.spline'] #extract R function# run smoothing function
					spline1 = r_smooth_spline(x=r_x, y=r_y, spar=.7)
					ySpline=np.array(robjects.r['predict'](spline1,robjects.FloatVector(x_range)).rx2('y'))
					RSS=spline1[9]
					ax[0].plot(x_range,ySpline,label=met,linewidth=3,color=colors[j])
				if isinstance(data2, np.ndarray) and data2.shape[0]>4:
					width=data2[1,0]-data2[0,0]
					xy = np.column_stack(( data2[:,0] + width , 100.*(data2[:,1]/width) ))
					x_range = np.linspace(min(xy[:,0]),max(xy[:,0]), num=100, endpoint=True)
					r_y = robjects.FloatVector(xy[:,1])
					r_x = robjects.FloatVector(xy[:,0])
		
					r_smooth_spline = robjects.r['smooth.spline'] #extract R function# run smoothing function
					spline1 = r_smooth_spline(x=r_x, y=r_y, spar=.7)
					ySpline=np.array(robjects.r['predict'](spline1,robjects.FloatVector(x_range)).rx2('y'))
					RSS=spline1[9]
					ax[1].plot(x_range,ySpline,label=met,linewidth=3,color=colors[j])
				j+=1
			ax[0].axhline(100.0*float(net_fts[3]),linestyle='dashed',color='black',linewidth=3,label='Random predictor')
			ax[1].axhline(100.0*float(net_fts[3]),linestyle='dashed',color='black',linewidth=3,label='Random predictor')
			#if i==len(ground_truth_labels)-1:
			#ax[0].legend(frameon=False,bbox_to_anchor=(1,1))
			#ax[1].legend(frameon=False,bbox_to_anchor=(1,1))
			#ax[0].set_ylim(ax[1].get_ylim())
			ax[1].set_xlim(ax[0].get_xlim())
			ax[0].tick_params(axis='both', which='major', labelsize=14)
			ax[0].tick_params(axis='both', which='minor', labelsize=14)
			ax[1].tick_params(axis='both', which='major', labelsize=14)
			ax[1].tick_params(axis='both', which='minor', labelsize=14)
			i+=1
		fig.tight_layout()
		
		fig.savefig(output_folder+'rank_'+'_'+ct+'_'+set_lab+'_'+at+'_'+meas+'.pdf',bbox_inches='tight'),plt.close()

def main():
    opts = parse_arguments()
    
    res_folder=opts.res_folder
    input_folder=opts.input_folder
    output_folder=res_folder+'ranking_'+opts.gtlab+'/'
    
    if os.path.isdir(output_folder)==False:
        os.mkdir(output_folder)
    
    catflag=False
	
    methods=['PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR']

    celltype=['HepG2_Smartseq2','HepG2_DNBelab','HepG2_9CL_SCAN','K562_CELseq','K562_STORMseq','K562_Smartseq3','K562_9CL_SCAN','K562_UMI200_SCAN']
    set_labs=['RBPs_RNA500_spec','RBPs_RNA1000_spec']
    set_labs2=['RBP_RNA500','RBP_RNA1000']
    
    ground_truth_labels=[opts.gtlab]
    my_index = [[i] * len(celltype) for i in ground_truth_labels]
    my_index=[item for sublist in my_index for item in sublist]
    for (set_lab,set_lab2) in zip(set_labs,set_labs2):
        for at in ['canonical']:
            PlotRankings(res_folder,input_folder,celltype,set_lab,ground_truth_labels,methods,(10,4*len(ground_truth_labels)),set_lab2,output_folder,at,'Ranking')
if __name__ == '__main__':
  main()

