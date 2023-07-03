import numpy as np
import pandas as pd
import os
from statsmodels.distributions.empirical_distribution import ECDF

import argparse

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Make benchmark heatmap.')
        
    parser.add_argument('-d','--downf', default='mRNA_downsampling_results',
        help="Benchmark measure to be used")
        
    parser.add_argument('-f','--folder', default='mRNA_lncRNA',
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


opts = parse_arguments()
myfolder=opts.folder
downfolder=opts.downf
base_out_folder='../Inference/outputs/'+myfolder+'/'
base_in_folder='../Inference/inputs/'+myfolder+'/'
downsampling_folder='../Inference/outputs/'+downfolder+'/'
    
datasets=['HepG2_Smartseq2','HepG2_DNBelab','K562_CELseq','K562_STORMseq','K562_Smartseq3']
reg='RBP_mRNA400'

for d in datasets:
    ct=d.split('_')[0]
    gt_string='eCLIP_'+ct+'_single_FC_1_P_3_canonical_stat_'+d+'_'+reg+'_'
    
    mRNA_out_folder=base_out_folder+'eCLIP_'+ct+'_single_FC_1_P_3_canonical'+d+'_RBP_mRNA400'
    mRNA_in_folder=base_in_folder+d+'_RBP_mRNA400/'
    
    # Read the EPr and stat of the corresponding lnc400 dataset
    lnc_out_folder=base_out_folder+'eCLIP_'+ct+'_single_FC_1_P_3_canonical'+d+'_RBP_lnc400'
    lnc_in_folder=base_in_folder+d+'_RBP_lnc400/'
    stat_lnc=np.loadtxt(lnc_in_folder+'eCLIP_'+ct+'_single_FC_1_P_3_canonical_stat_tgt.txt')
    density_lnc=stat_lnc[3]
    EPr_lnc=pd.read_csv(lnc_out_folder+'/EPrTgt.csv',index_col=0)
    EPR_lnc=EPr_lnc/float(density_lnc)
    EPr_cfilt_lnc=pd.read_csv(lnc_out_folder+'_cfilt/EPrTgt.csv',index_col=0)
    EPR_cfilt_lnc=EPr_cfilt_lnc/float(density_lnc)
    
    # Now read the EPr and stat for each downsampled mRNA net
    EPr_list=[]
    EPr_cfilt_list=[]
    EPR_list=[]
    EPR_cfilt_list=[]
    stats_list=[]
    for k in range(50):
        stat=np.loadtxt(downsampling_folder+gt_string+str(k+1)+'.txt')
        density=stat[3]
        EPr=pd.read_csv(downsampling_folder+d+'_'+reg+'_EPrTgt_'+str(k+1)+'.csv',index_col=0)
        EPR=EPr/float(density)
        EPr_cfilt=pd.read_csv(downsampling_folder+d+'_'+reg+'_EPrTgt_cfilt_'+str(k+1)+'.csv',index_col=0)
        EPR_cfilt=EPr_cfilt/float(density)
        
        EPr_list.append(EPr.T)
        EPr_cfilt_list.append(EPr_cfilt.T)
        EPR_list.append(EPR.T)
        EPR_cfilt_list.append(EPR_cfilt.T)

        stats_list.append(stat)
    
    EPr_df=pd.concat(EPr_list)
    EPr_cfilt_df=pd.concat(EPr_cfilt_list)
    EPR_df=pd.concat(EPR_list)
    EPR_cfilt_df=pd.concat(EPR_cfilt_list)
    
    # Average each column of the dataframes
    EPr_mean=EPr_df.mean(axis=0).T
    EPr_cfilt_mean=EPr_cfilt_df.mean(axis=0).T
    EPR_mean=EPR_df.mean(axis=0).T
    EPR_cfilt_mean=EPR_cfilt_df.mean(axis=0).T
    
    EPr_mean.to_csv(mRNA_out_folder+'/EPrTgtDown.csv')
    EPR_mean.to_csv(mRNA_out_folder+'/EPRTgtDown.csv')
    EPr_cfilt_mean.to_csv(mRNA_out_folder+'_cfilt/EPrTgtDown.csv')
    EPR_cfilt_mean.to_csv(mRNA_out_folder+'_cfilt/EPRTgtDown.csv')
    
    # Compute the ECDF for each algorithm
    EPr_ECDF_df=EPr_df.apply(ECDF)
    EPr_cfilt_ECDF_df=EPr_cfilt_df.apply(ECDF)
    EPR_ECDF_df=EPR_df.apply(ECDF)
    EPR_cfilt_ECDF_df=EPR_cfilt_df.apply(ECDF)
    
    p_EPr=[1-EPr_ECDF_df.loc[i](EPr_lnc.loc[i]) for i in EPr_ECDF_df.index]
    p_EPr_df=pd.DataFrame(data={'p-value': p_EPr},index=EPr_ECDF_df.index)
    
    p_EPR=[1-EPR_ECDF_df.loc[i](EPR_lnc.loc[i]) for i in EPR_ECDF_df.index]
    p_EPR_df=pd.DataFrame(data={'p-value': p_EPR},index=EPR_ECDF_df.index)
    
    p_EPr_cfilt=[1-EPr_cfilt_ECDF_df.loc[i](EPr_cfilt_lnc.loc[i]) for i in EPr_cfilt_ECDF_df.index]
    p_EPr_cfilt_df=pd.DataFrame(data={'p-value': p_EPr_cfilt},index=EPr_cfilt_ECDF_df.index)
    
    p_EPR_cfilt=[1-EPR_cfilt_ECDF_df.loc[i](EPR_cfilt_lnc.loc[i]) for i in EPR_cfilt_ECDF_df.index]
    p_EPR_cfilt_df=pd.DataFrame(data={'p-value': p_EPR_cfilt},index=EPR_cfilt_ECDF_df.index)
    
    p_EPr_df.to_csv(mRNA_out_folder+'/p_EPrTgtDown.csv')
    p_EPR_df.to_csv(mRNA_out_folder+'/p_EPRTgtDown.csv')
    p_EPr_cfilt_df.to_csv(mRNA_out_folder+'_cfilt/p_EPrTgtDown.csv')
    p_EPR_cfilt_df.to_csv(mRNA_out_folder+'_cfilt/p_EPRTgtDown.csv')
    
    stats_df=pd.DataFrame(list(map(np.ravel, stats_list)))
    stats_df.columns=['reg','tgt','edges','d']
    stats_mean=stats_df.mean(axis=0)
    stats_mean.values
    np.savetxt(mRNA_in_folder+'eCLIP_'+ct+'_single_FC_1_P_3_canonical_stat_tgt_down.txt',
               np.c_[stats_mean.values],fmt='%f')
