import pandas as pd
import argparse
import numpy as np
from os import walk
import re
import time

def get_sec(time_str):
    if len(time_str.split(':'))==2:
        m,s=time_str.split(':')
        return int(m)*60 + float(s)
    elif len(time_str.split(':'))==3:
        h,m,s=time_str.split(':')
        return int(h)*60*60+int(m)*60+float(s)


parser = argparse.ArgumentParser()
parser.add_argument('--out_folder', type=str, default='/out')
opt = parser.parse_args()

res = []
for i in range(10):
    res.append(pd.read_csv(opt.out_folder+'/res_' + str(i+1)+ '.tsv',sep='\t'))
res = pd.concat(res)
res['EdgeWeight'] = abs(res['EdgeWeight'])
res=res.groupby(['TF','Target']).mean()

res=res.reset_index()


print('Renaming columns')
res.columns=['Gene1','Gene2','EdgeWeight']
print('Sorting according to EdgeWeight')
res.sort_values('EdgeWeight',ascending=False,inplace=True)
print('Saving')
res.to_csv(opt.out_folder+'/rankedEdges.csv',index=False)
