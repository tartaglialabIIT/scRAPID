import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--folder', type=str)
parser.add_argument('--filename', type=str)
opt = parser.parse_args()


# In[ ]:

folder=opt.folder
filename=opt.filename

data=pd.read_csv(folder+filename,index_col=0)
data=data.T
data.to_csv(folder+filename)
