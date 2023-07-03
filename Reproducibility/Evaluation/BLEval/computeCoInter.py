import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
import itertools
from tqdm import tqdm

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

def CoInter(datasetDict, inputSettings,TFEdges=True):
    '''
    Computes protein co-interactions based on shared RNA targets
    :param datasetDict:   A dictionary containing the dataset name, path to reference network.
    :type datasetDict: dict
    :param inputSettings: An object of class :class:`BLEval.InputSettings`.
    :type inputSettings: :class:`BLEval.InputSettings`
    :returns:
    '''
    
    Norm_Data=pd.read_csv(str(inputSettings.datadir)+'/'+ datasetDict['name'] +'/' +datasetDict['exprData'],sep = ',',header = 0, index_col = 0)
    all_genes=list(Norm_Data.index)
    possibleEdges=(len(all_genes)*len(all_genes))-len(all_genes)                            

    # provide a list of proteins for which you want to compute co-interactions
    # CHANGE PATH HERE
    bioplex_proteins=list(np.loadtxt("/mnt/large/jfiorentino/INTERACTomics/Beeline/Bioplex/Bioplex_HEK293T_RBPs.txt",dtype=str))
    
    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' + datasetDict['name']
    dataDict = {}

    dataDict['CO'] = {}
    
    for algo in tqdm(inputSettings.algorithms,
                     total = len(inputSettings.algorithms), unit = " Algorithms"):
        # if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
        #     continue
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', sep = '\t', header =  0, index_col = None)


            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            predDF.sort_values('EdgeWeight', ascending=False,inplace=True)
            
            if len(list(predDF.index))>int(0.05*possibleEdges):
                predDF=predDF.iloc[:int(0.05*possibleEdges)]
            
            # Select edges going out of Bioplex proteins
            predDF=predDF.loc[predDF['Gene1'].isin(bioplex_proteins)]
            

            if not predDF.shape[0] == 0:

                
                
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
                
               
                dataDict['CO'][algo[0]]= algJ
            
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')

            dataDict['CO'][algo[0]] = 0

    dataDF = pd.DataFrame(dataDict)
    print("success")
    return dataDF['CO'];
