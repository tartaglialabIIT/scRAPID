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

    parser.add_argument('-m','--measure', default='EPRTgt',
        help="Benchmark measure to be used")
        
    parser.add_argument('-r','--res_folder', default='/mnt/large/jfiorentino/INTERACTomics/Beeline/outputs/mRNA_lncRNA/',
        help="Folder with the results")
        
    parser.add_argument('-i','--input_folder', default='/mnt/large/jfiorentino/INTERACTomics/Beeline/inputs/mRNA_lncRNA/',
        help="Folder with the results")
        
    parser.add_argument('-g','--gtlab', default='eCLIP',
        help="Folder with the results")
        
    parser.add_argument('-d','--downf', default='mRNA_downsampling_results',
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

def add_cosmetics(ax,title='mRNA-lncRNA comparison', 
                  xlabel='', ylabel='EPR',showxticks=False):
    ax.set_title('', fontsize=16)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.tick_params(axis="x", labelsize=15,rotation=45) 
    ax.tick_params(axis="y", labelsize=15) 
    if showxticks==False:
        ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    sns.despine()

def MakeBoxPlot(res_folder,input_folder,reg,gtlab,ctypes,mets,meas,catRAPID,at,output_folder,downsampling_folder):
	#df_ratios= pd.DataFrame(np.nan, index=ctypes, columns=mets)
	old_gt_lab=gtlab[:]
	fig,ax=plt.subplots(5,1,figsize=(15, 25))
	i=0

	for ct in ctypes:
		ctn=ct.split('_')[0]
		gt_string='eCLIP_'+ctn+'_single_FC_1_P_3_canonical_stat_'+ct+'_'+reg+'_'
		gtlab='eCLIP_'+ctn+'_single_FC_1_P_3_'+at
		AUC_data=[]
		if catRAPID=='catrapid':
			for k in range(100):
				stat=np.loadtxt(downsampling_folder+gt_string+str(k+1)+'.txt')
				density=stat[3]
				EPr=pd.read_csv(downsampling_folder+ct+'_'+reg+'_EPrTgt_cfilt_'+str(k+1)+'.csv',index_col=0)
				EPr.columns=['EPR']
				AUC_data.append(EPr.T/float(density))
			pvals=pd.read_csv(res_folder+gtlab+ct+'_'+reg+'_cfilt/p_'+meas+'Down.csv',index_col=0)
			AUC_lnc_data=pd.read_csv(res_folder+gtlab+ct+'_RBP_lnc400_cfilt/EPrTgt.csv',index_col=0)
			AUC_lnc_data.columns=['EPR']
			AUC_mRNA_data=pd.read_csv(res_folder+gtlab+ct+'_RBP_mRNA400_cfilt/EPr.csv',index_col=0)
			AUC_mRNA_data.columns=['EPR']
		else:
			for k in range(100):
				stat=np.loadtxt(downsampling_folder+gt_string+str(k+1)+'.txt')
				density=stat[3]
				EPr=pd.read_csv(downsampling_folder+ct+'_'+reg+'_EPrTgt_'+str(k+1)+'.csv',index_col=0)
				EPr.columns=['EPR']
				AUC_data.append(EPr.T/float(density))
			pvals=pd.read_csv(res_folder+gtlab+ct+'_'+reg+'/p_'+meas+'Down.csv',index_col=0)
			AUC_lnc_data=pd.read_csv(res_folder+gtlab+ct+'_RBP_lnc400/EPrTgt.csv',index_col=0)
			AUC_lnc_data.columns=['EPR']
			AUC_mRNA_data=pd.read_csv(res_folder+gtlab+ct+'_RBP_mRNA400/EPr.csv',index_col=0)
			AUC_mRNA_data.columns=['EPR']
		
		
		net_fts_lnc=np.loadtxt(input_folder+ct+'_RBP_lnc400/'+gtlab+'_stat_tgt.txt',dtype=str)
		AUC_lnc_data=AUC_lnc_data/float(net_fts_lnc[3])
		AUC_lnc_data=AUC_lnc_data.loc[mets]
		
		net_fts_mrna=np.loadtxt(input_folder+ct+'_RBP_mRNA400/'+gtlab+'_stat.txt',dtype=str)
		AUC_mRNA_data=AUC_mRNA_data/float(net_fts_mrna[3])
		AUC_mRNA_data=AUC_mRNA_data.loc[mets]
		
		pvals=pvals.loc[mets]
		pvals['p-value']=[i.replace('[','') for i in list(pvals['p-value'])]
		pvals['p-value']=[i.replace(']','') for i in list(pvals['p-value'])]
		pvals['p-value']=[float(i) for i in list(pvals['p-value'])]

		AUC_data=pd.concat(AUC_data)
		
		print(pvals)
		
		algo=[]
		EPR=[]
		for met in mets:
			if met in list(set(AUC_data.columns)):
				if met=='DeePSEM_CT':
					meto='DeePSEM'
				elif met=='ARACNe_FDR':
					meto='ARACNe'
				else:
					meto=met
			algo.append([meto]*len(AUC_data.loc[:,met]))
			EPR.append(list(AUC_data.loc[:,met]))
		
		df=pd.DataFrame(data={'method':[item for sublist in algo for item in sublist], 'EPR':[item for sublist in EPR for item in sublist]})
				
		# Create violin plots without mini-boxplots inside.
		#sns.violinplot(y='EPR', x='method', data=df,color='mediumslateblue', cut=0, inner=None,ax=ax[i])
		# Clip the right half of each violin.
		#for item in ax[i].collections:
		#	x0, y0, width, height = item.get_paths()[0].get_extents().bounds
		#	item.set_clip_path(plt.Rectangle((x0, y0), width/2, height,transform=ax[i].transData))
		# Create strip plots with partially transparent points of different colors depending on the group.
		#num_items = len(ax[i].collections)
		#sns.stripplot(y='EPR', x='method', color='black', data=df,alpha=0.4, size=7,ax=ax[i])
		# Shift each strip plot strictly below the correponding volin.
		#for item in ax[i].collections[num_items:]:
	#		item.set_offsets(item.get_offsets() + 0.15)
		# Create narrow boxplots on top of the corresponding violin and strip plots, with thick lines, the mean values, without the outliers.
		sns.boxplot(y='EPR', x='method', data=df,ax=ax[i], width=0.25,showfliers=False, showmeans=True, meanprops=dict(marker='o', markerfacecolor='darkorange',markersize=10, zorder=3),boxprops=dict(facecolor=(0,0,0,0), linewidth=3, zorder=3),whiskerprops=dict(linewidth=3),capprops=dict(linewidth=3),medianprops=dict(linewidth=3))
		if (i==len(ctypes)-1):
			add_cosmetics(ax[i],xlabel='', ylabel='EPR',showxticks=True)
		else:
			add_cosmetics(ax[i],xlabel='', ylabel='EPR')
		ax[i].axhline(1,linestyle='dashed',color='black',lw=3)
		
		# Add the EPR for long non-coding RNAs
		ax[i].scatter(np.arange(0,len(mets)),list(AUC_lnc_data['EPR']),marker='s',color='red',s=100,label='EPR lncRNA',zorder=100)
		ax[i].scatter(np.arange(0,len(mets)),list(AUC_mRNA_data['EPR']),marker="^",color='cyan',s=100,label='EPR mRNA (original)',zorder=101)
		
		# Add significance level
		r=0
		for l in mets:
			# Find max between the EPR of the lncRNAs and all the EPRs
			if met=='DeePSEM_CT':
				meto='DeePSEM'
			elif met=='ARACNe_FDR':
				meto='ARACNe'
			else:
				meto=met
			max_mRNA=df[df.method==meto]['EPR'].max()
			totMax=max(max_mRNA,AUC_lnc_data.loc[l][0])
			print(pvals.loc[l][0],type(pvals.loc[l][0]))
			if pvals.loc[l][0] < 0.05:
				ax[i].scatter(r,totMax+0.2,marker='*',color='black',s=70)
			r+=1			
		ax[i].set_title(ct, fontsize=16)
		if i==0:
			ax[i].legend(frameon=False, fontsize=15,loc=0)
			
		
		i+=1
	fig.tight_layout()
	if catRAPID=='catrapid':
		plt.savefig(output_folder+meas+'_'+at+'_cfilt_boxplot.pdf',bbox_inches='tight'),plt.close()
	else:
		plt.savefig(output_folder+meas+'_'+at+'_boxplot.pdf',bbox_inches='tight'),plt.close()
    
def main():
    opts = parse_arguments()
    
    res_folder=opts.res_folder
    input_folder=opts.input_folder
    output_folder1=res_folder+'boxplot_'+opts.gtlab+'Down/'
    downsampling_folder='../Inference/outputs/'+opts.downf+'/'    
    
    if os.path.isdir(output_folder1)==False:
        os.mkdir(output_folder1)

    
    # # For each ground truth dataset and each type (biased or high_exp)
	# # create a dataframe with the methods on the columns, the cell lines on the rows.
	# # The values are EP, AUPRC and AUROC or their ratios with that of a random predictor
	# # i.e. the ground truth network density

	# # Mark with grey boxes the missing values
	# # Mark with black boxes the worse than random values
    methods=['PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_B','DeePSEM','ARACNe_FDR']

    celltype=['HepG2_Smartseq2','HepG2_DNBelab','K562_CELseq','K562_STORMseq','K562_Smartseq3']
    set_labs=['RBPs_mRNA400_spec']
    set_labs2=['RBP_mRNA400']
    
    ground_truth_labels=[opts.gtlab]
    #ground_truth_labels=['eCLIP_loose','shRNA_seq_indirect_strict']
    my_index = [[i] * len(celltype) for i in ground_truth_labels]
    my_index=[item for sublist in my_index for item in sublist]
    for (set_lab,set_lab2) in zip(set_labs,set_labs2):
        for at in ['canonical']:
            MakeBoxPlot(res_folder,input_folder,set_lab2,ground_truth_labels[0],celltype,methods,opts.measure,False,at,output_folder1,downsampling_folder)
            MakeBoxPlot(res_folder,input_folder,set_lab2,ground_truth_labels[0],celltype,methods,opts.measure,'catrapid',at,output_folder1,downsampling_folder)
               	
if __name__ == '__main__':
  main()

