#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 14})


# This uses the output file produced by the R script that makes the bar plots of the EPR
data500=pd.read_csv("./TF_RBP_500.csv")

data500_orig=data500[data500['Reg_Type']=='TF'].copy()
data500_cat=data500[data500['Reg_Type']=='RBP'].copy()

from scipy.stats import kstest
print(kstest(data500_orig[data500_orig.EPR>0].EPR,data500_cat[data500_cat.EPR>0].EPR))


fig,ax=plt.subplots(figsize=(6,4))
ax.set_xlabel('EPR',fontsize=12)
ax.set_ylabel('Probability Density',fontsize=12)
sns.kdeplot(data500_orig[data500_orig.EPR>0].EPR, color='#a1e9f0', shade=True, label='TF')
sns.kdeplot(data500_cat[data500_cat.EPR>0].EPR, color='#d9b1f0', shade=True, label='RBP')
ax.axvline(data500_orig[data500_orig.EPR>0].EPR.mean(),color='#a1e9f0',lw=3)
ax.axvline(data500_cat[data500_cat.EPR>0].EPR.mean(),color='#d9b1f0',lw=3)
ax.axvline(1,color='black',lw=3,linestyle='dashed',label='Random pred.')
ax.legend(loc=0,frameon=False)
fig.tight_layout()
plt.savefig("TF_RBP_500.pdf",bbox_inches='tight')


data1000=pd.read_csv("./TF_RBP_1000.csv")

data1000_orig=data1000[data1000['Reg_Type']=='TF'].copy()
data1000_cat=data1000[data1000['Reg_Type']=='RBP'].copy()


from scipy.stats import kstest
print(kstest(data1000_orig[data1000_orig.EPR>0].EPR,data1000_cat[data1000_cat.EPR>0].EPR))


fig,ax=plt.subplots(figsize=(6,4))
ax.set_xlabel('EPR',fontsize=12)
ax.set_ylabel('Probability Density',fontsize=12)
sns.kdeplot(data1000_orig[data1000_orig.EPR>0].EPR, color='#a1e9f0', shade=True, label='TF')
sns.kdeplot(data1000_cat[data1000_cat.EPR>0].EPR, color='#d9b1f0', shade=True, label='RBP')
ax.axvline(data1000_orig[data1000_orig.EPR>0].EPR.mean(),color='#a1e9f0',lw=3)
ax.axvline(data1000_cat[data1000_cat.EPR>0].EPR.mean(),color='#d9b1f0',lw=3)
ax.axvline(1,color='black',lw=3,linestyle='dashed',label='Random pred.')
ax.legend(loc=0,frameon=False)
fig.tight_layout()
plt.savefig("TF_RBP_1000.pdf",bbox_inches='tight')



