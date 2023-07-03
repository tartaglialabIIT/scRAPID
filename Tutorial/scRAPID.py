import pandas as pd
import numpy as np
import itertools
from itertools import product,permutations
import networkx as nx

def preprocess(folder):
    predDF = pd.read_csv(folder+"rankedEdges.csv", sep = '\t', header =  0, index_col = None)
    
    # Remove self-loops
    predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
    
    # Remove NaN (if present)
    predDF = predDF.dropna(subset='EdgeWeight')

    # Sort edges according to the EdgeWeight
    predDF.sort_values('EdgeWeight', ascending=False,inplace=True)
    
    # Remove duplicated edges (if present)
    predDF.drop_duplicates(keep = 'first', inplace=True)
    
    predDF.reset_index(drop = True,  inplace= True)
    
    return predDF.loc[:,['Gene1','Gene2','EdgeWeight']];

def CoInter(predDF):
	gene_groups = predDF.groupby('Gene1')
	gene_targets = {gene: set(gene_group['Gene2']) for gene, gene_group in gene_groups}
	
	regList = list(set(predDF.Gene1))
	regCombinations = itertools.combinations(regList, 2)
	algJ = []
	
	for reg1, reg2 in regCombinations:
		targets1 = gene_targets[reg1]
		targets2 = gene_targets[reg2]
		j = len(targets1 & targets2) / len(targets1 | targets2)
		algJ.append([reg1, reg2, len(targets1), len(targets2), j])
		
	df=pd.DataFrame(algJ)
	df.columns=['RBP1','RBP2','nr_Tgt1','nr_Tgt2','Jaccard']
	return df;
	
def FilterRanking(predDF,catRAPIDDF, threshold = 30):
	predDF['Edges']=predDF['Gene1']+'|'+predDF['Gene2']
	predDF=predDF.set_index('Edges')
	
	# Filter the inferred edges with catRAPID
	i1 = catRAPIDDF[catRAPIDDF['score']>float(threshold)].index
	i2 = predDF.index
	predDF_filtered=predDF[i2.isin(i1)]
	predDF_filtered=predDF_filtered.loc[:,['Gene1','Gene2','EdgeWeight']]
	predDF_filtered.sort_values('EdgeWeight',inplace=True,ascending=False)
	
	return predDF_filtered;	
	
	
def HubReg(predDF,regs,algo):
	predGraph = nx.from_pandas_edgelist(predDF,source='Gene1',
                                   target='Gene2',
                                   create_using=nx.DiGraph())
	
	if algo=='PPCOR' or algo=='PIDC':
		inGraphU = nx.DiGraph.to_undirected(predGraph)
	else:
		inGraphU = predGraph.copy()
	deg_centr = nx.degree_centrality(inGraphU)
	deg_centr_s = sorted(deg_centr.items(), key=lambda kv: kv[1],reverse=True)
	#print(deg_centr_s)
	deg_centr_s = [i for i in deg_centr_s if i[0] in list(regs)]
	df = pd.DataFrame(deg_centr_s, columns=['RBP', 'outdeg_centr'])

	return df;
    
def HubTarget(predDF,tgts,algo):
	predGraph = nx.from_pandas_edgelist(predDF,source='Gene1',
                                   target='Gene2',
                                   create_using=nx.DiGraph())
	if algo=='PPCOR' or algo=='PIDC':
		inGraphU = nx.DiGraph.to_undirected(predGraph)
		ideg_centr = nx.degree_centrality(predGraph)
	else:
		ideg_centr = nx.in_degree_centrality(predGraph)
    
	ideg_centr_s = sorted(ideg_centr.items(), key=lambda kv: kv[1],reverse=True)
	ideg_centr_s = [i for i in ideg_centr_s if i[0] in list(tgts)]
	df = pd.DataFrame(ideg_centr_s, columns=['RNA', 'indeg_centr'])

	return df;
