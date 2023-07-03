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
import math

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

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts


def MakeDf(res_folder,input_folder,reg,gtlab,ctypes,mets,meas,catRAPID,at):
	df_ratios= pd.DataFrame(np.nan, index=ctypes, columns=mets)
	df_labels= pd.DataFrame(np.nan, index=ctypes, columns=mets)
	net_fts_lst=[]
	old_gt_lab=gtlab[:]
	for ct in ctypes:
		ctn=ct.split('_')[0]
		gtlab='eCLIP_'+ctn+'_single_FC_1_P_3_'+at  
		AUC_data=pd.read_csv(res_folder+gtlab+ct+'_'+reg+'_cfilt/'+meas+'.csv')
		
		AUC_data_orig=pd.read_csv(res_folder+gtlab+ct+'_'+reg+'/'+meas+'.csv')
		AUC_data.columns=['Method',meas]
		AUC_data_orig.columns=['Method',meas]
		#print(gtlab)
		if meas=='EPr':
			net_fts=np.loadtxt(input_folder+ct+'_'+reg+'/'+gtlab+'_stat.txt',dtype=str)
		elif meas=='EPrReg':
			net_fts=np.loadtxt(input_folder+ct+'_'+reg+'/'+gtlab+'_stat_reg.txt',dtype=str)
		elif meas=='EPrTgt':
			net_fts=np.loadtxt(input_folder+ct+'_'+reg+'/'+gtlab+'_stat_tgt.txt',dtype=str)
		net_fts_lst.append(net_fts)
		for met in mets:
			if met in list(AUC_data.Method):
				EPR1=AUC_data_orig[AUC_data_orig.Method==met][meas].iloc[0]/float(net_fts[3])
				EPR2=AUC_data[AUC_data.Method==met][meas].iloc[0]/float(net_fts[3])
				df_ratios.loc[ct,met]=AUC_data[AUC_data.Method==met][meas].iloc[0]/float(net_fts[3])
				df_labels.loc[ct,met]='('+str("{0:.2g}".format(100*(EPR2-EPR1)/EPR1))+'%'+')'
				# Worse than random predictor: assign -2
				# Missing results filled with NaN
				if df_ratios.loc[ct,met]< 1.0:
					df_ratios.loc[ct,met]= -2
	return df_ratios,df_labels,net_fts_lst;

def round_to_nsf(value, nsf=2):
    """
    Rounds the number to the provided number of significant figures.
    """
    if value!=np.nan: 
        integer_part = math.floor(value)
        if integer_part > 0:
            integer_part_len = len(str(integer_part))
            return round(value, nsf-integer_part_len)
        else:
            str_value = str(value)
            #if of the form "8e-05"
            if '-' in str_value:
                index = int(str_value[str_value.find('-')+1:]) - 1
            else:
               st = {'1', '2', '3', '4', '5', '6', '7', '8', '9'}
               index = next((i for i, ch in enumerate(str(value)) if ch in st), None) - 2
            return round(value, index+nsf)
    else:
        return np.nan;
                   
def MakeHeatmap(res_folder,df,meas,set_type,mysize,save,catflag,output_folder,at,labels_df):
    cmap = sns.cm.flare_r
    df2=df.copy()
    df2[df2 == -2] = np.nan
    x = df2.values #returns a numpy array
    x_min = np.nanmin(x,axis=1)
    x_max = np.nanmax(x,axis=1)
    x_scaled = (x-x_min[:,np.newaxis])/(x_max[:,np.newaxis]-x_min[:,np.newaxis])
    nans_in_rows=df2.isnull().sum(axis=1)
    for i in range(len(nans_in_rows)):
        # If there is only one value that is not nan,
        # assign 0.5 to the scaled values
        if x_min[i]==x_max[i]:
            cols=df2.iloc[i,:].notnull()
            x_scaled[i,cols]=0.5
        if nans_in_rows.iloc[i] == len(df2.columns)-1:
            col=df2.iloc[i,:].notnull()
            x_scaled[i,col]=0.5

    df2 = pd.DataFrame(x_scaled,index=df.index,columns=df.columns)
    df2[df== -2]=-2
    # combining text with values
    #annot_df=df.copy()
    #for i in range(len(labels_df)):
    #    for j in range(len(labels_df)):
    #        annot_df.iloc[i,j]=str(df.iloc[i,j])+'\n'+str(labels_df.iloc[i,j])
    formatted_text = (np.asarray(["{0}\n{1}".format(text, data) for (text, data) in zip(df.round(1).values.flatten(), labels_df.values.flatten())])).reshape(df.values.shape[0], df.values.shape[1])
    f,ax=plt.subplots(figsize=mysize)
    sns.heatmap(df2,ax=ax,mask=df.isnull() | (df == -2.0)
                ,cmap=cmap,cbar=True,annot=formatted_text,fmt='',
                square=True,linecolor='white', linewidths=0.5,
                cbar_kws={"orientation": "horizontal","shrink": .4})
              
    ax.tick_params(axis='both', which='major', labelsize=10,rotation=90,labelbottom = False, bottom=False, top = False, labeltop=True)
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(np.array([0.0, 1.0]))
    colorbar.set_ticklabels(['Low/Poor', 'High/Good'])

    sns.heatmap(df, cmap=ListedColormap(['black']), linecolor='white', linewidths=0.5,square=True, cbar=False, mask=(df != -2.0), ax=ax)
#     f.tight_layout()
    if save:
        if catflag=='catrapid':
            plt.savefig(output_folder+meas+'_'+set_type+'_'+at+'_cfilt_heatmap.pdf',bbox_inches='tight')
        else:
            plt.savefig(output_folder+meas+'_'+set_type+'_'+at+'_heatmap.pdf',bbox_inches='tight')
    plt.show(),plt.close()
    
def CreateDFandPlot(celltype, methods,set_lab,reg, ground_truth_labels,meas,my_index,res_folder,input_folder,catRAPID,output_folder,at):
    df=pd.DataFrame(columns=methods)
    net_fts_lst=[]
    for gt_lab in ground_truth_labels:
        tmp_df,labels_df,net_fts=MakeDf(res_folder,input_folder, reg ,gt_lab,celltype,methods,meas,catRAPID,at)
        df=pd.concat([df,tmp_df])
        net_fts_lst.append(net_fts)
    df.set_index([my_index,df.index],inplace=True)
    Newdf=pd.DataFrame(data=[item for sublist in net_fts_lst for item in sublist]).astype(float)
    Newdf[0]=Newdf[0].astype(int)
    Newdf[1]=Newdf[1].astype(int)
    Newdf[2]=Newdf[2].astype(int)
    Newdf['ct']=celltype*int(len(df)/len(celltype))
    Newdf=Newdf.loc[:,['ct',0,1,2,3]]
    Newdf=Newdf.style.format({  3: '{: .2}'})
    dfi.export(Newdf, output_folder+'df_'+meas+'_'+set_lab+'_'+reg+'_'+'_styled.png')
    MakeHeatmap(res_folder,df,meas,set_lab,(16,len(my_index)),True,catRAPID,output_folder,at,labels_df)
    return df;
    
def main():
    opts = parse_arguments()
    
    
    res_folder=opts.res_folder
    input_folder=opts.input_folder
    output_folder=res_folder+'heatmaps_catrapid/'
    
    if os.path.isdir(output_folder)==False:
        os.mkdir(output_folder)

    catflag='catrapid'
    
    methods=['PIDC','GRNBOOST2','SINCERITIES','TENET','TENET_A','TENET_B','DeePSEM','ARACNe_FDR']

    celltype=['HepG2_Smartseq2','HepG2_DNBelab','HepG2_9CL_SCAN','K562_CELseq','K562_STORMseq','K562_Smartseq3','K562_9CL_SCAN','K562_UMI200_SCAN']
    set_labs=['RBPs_RNA500_spec','RBPs_RNA1000_spec']
    set_labs2=['RBP_RNA500','RBP_RNA1000']
    
    ground_truth_labels=[opts.gtlab]
    my_index = [[i] * len(celltype) for i in ground_truth_labels]
    my_index=[item for sublist in my_index for item in sublist]
    for (set_lab,set_lab2) in zip(set_labs,set_labs2):
        for at in ['canonical']:
            df=CreateDFandPlot(celltype,methods,set_lab,set_lab2, ground_truth_labels,opts.measure,my_index,res_folder,input_folder,catflag,output_folder,at)
            #PlotRankings(res_folder,input_folder,celltype,set_lab,ground_truth_labels,methods,(14,4*len(ground_truth_labels)),set_lab2,output_folder3,at,'Ranking')

	
	

				
if __name__ == '__main__':
  main()

