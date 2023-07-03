import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from itertools import product, permutations, combinations, combinations_with_replacement
from tqdm import tqdm
import networkx as nx
#import multiprocessing as mp

#def RandomJaccard(N,M,m):
#    rGr = nx.gnm_random_graph(N, M, directed=False)
#    rGrd = nx.gnm_random_graph(N, M, directed=True)
#    rGr = nx.relabel_nodes(rGr, m)
#    rGrd = nx.relabel_nodes(rGrd, m)
#    r_vals = np.array(Hubs_undirected(rGr,refGraph)[:6])
#    rd_vals = np.array(Hubs_directed(rGrd,refGraph)[:6])
#    
#    return r_vals, rd_vals;

def HubsTargets(datasetDict, inputSettings,TFEdges=True):
    '''
    Computes hub similarity. I  use in- and out- degree centrality, see STREAMLINE (https://www.biorxiv.org/content/10.1101/2022.10.31.514493v3) 
    for the original functions computing all measures of hub similarities. We exclude RNAs coding for eCLIP RBPs. 
    :param datasetDict:   A dictionary containing the dataset name, path to reference network.
    :type datasetDict: dict
    :param inputSettings: An object of class :class:`BLEval.InputSettings`.
    :type inputSettings: :class:`BLEval.InputSettings`
    :returns:
    '''

    # Read file for trueEdges
    trueEdgesDF = pd.read_csv(str(inputSettings.datadir)+'/'+ datasetDict['name'] +
                                '/' +datasetDict['trueEdges'],
                                sep = ',',
                                header = 0, index_col = None)
    
    # Load the RBPs from the two cell lines
    RBPs_1=np.loadtxt('RBPs_HepG2.txt',dtype=str)
    RBPs_2=np.loadtxt('RBPs_K562.txt',dtype=str)
        
    RBPs=list(set(list(RBPs_1)+list(RBPs_2)))

    trueEdgesDF=trueEdgesDF.loc[~trueEdgesDF.Gene2.isin(set(RBPs))]
    uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
    possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)-set(RBPs)))
    #possibleEdges = list(permutations(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]),
    #                         r = 2))
    # Get a list of all possible interactions 
    possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))
    # Find intersection of above lists to ignore self edges
    # TODO: is there a better way of doing this?
    possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)
    EdgeDict = {'|'.join(p):0 for p in possibleEdges}

    refGraph = nx.DiGraph()

    for key in EdgeDict.keys():
        u = key.split('|')[0]
        v = key.split('|')[1]
        if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == u) &
               (trueEdgesDF['Gene2'] == v)])>0:
                refGraph.add_edge(u,v)

    numEdges = len(refGraph.edges())
    numNodes = refGraph.number_of_nodes()

    refPR, refIC, refOC, refEIG, refBT, refRD,REFprhubs1,REFprhubs2,REFidhubs1,REFidhubs2,REFodhubs1,REFodhubs2,REFbthubs1,REFbthubs2,REFradhubs1,REFradhubs2 = Hubs_directed(refGraph, refGraph,set(trueEdgesDF.Gene1),set(trueEdgesDF.Gene2)-set(RBPs))

    # set-up outDir that stores output directory name
    outDir = "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/' + datasetDict['name']
    dataDict = {}

    dataDict['PR'] = {}
    dataDict['IC'] = {}
    dataDict['OC'] = {}
    dataDict['EIG'] = {}
    dataDict['BT'] = {}
    dataDict['RD'] = {}

    dataDict['PR']["ground truth"] = refPR
    dataDict['IC']["ground truth"] = refIC
    dataDict['OC']["ground truth"] = refOC
    dataDict['EIG']["ground truth"] = refEIG
    dataDict['BT']["ground truth"] = refBT
    dataDict['RD']["ground truth"] = refRD

    hubsDict = {}

    hubsDict['PR'] = {}
    hubsDict['IC'] = {}
    hubsDict['OC'] = {}
    hubsDict['BT'] = {}
    hubsDict['RD'] = {}

    hubsDict['PR']["ground truth"] = REFprhubs1
    hubsDict['IC']["ground truth"] = REFidhubs1
    hubsDict['OC']["ground truth"] = REFodhubs1
    hubsDict['BT']["ground truth"] = REFbthubs1
    hubsDict['RD']["ground truth"] = REFradhubs1

    for algo in tqdm(inputSettings.algorithms,
                     total = len(inputSettings.algorithms), unit = " Algorithms"):
        # if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
        #     continue
        # check if the output rankedEdges file exists
        if Path(outDir + '/' +algo[0]+'/rankedEdges.csv').exists():
             # Initialize Precsion

            predDF = pd.read_csv(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                                        sep = '\t', header =  0, index_col = None)


            predDF = predDF.loc[(predDF['Gene1'] != predDF['Gene2'])]
            predDF.drop_duplicates(keep = 'first', inplace=True)
            predDF.reset_index(drop = True,  inplace= True)
            if TFEdges:
                # Consider only edges going out of TFs
                # Get a list of all possible TF to gene interactions 
                #trueEdgesDF=trueEdgesDF.loc[~trueEdgesDF.Gene2.isin(set(trueEdgesDF.Gene1))]
                uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
                possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)-set(RBPs)))
                
                #uniqueNodes = np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']])
                #possibleEdges_TF = set(product(set(trueEdgesDF.Gene1),set(uniqueNodes)))
                # Get a list of all possible interactions 
                possibleEdges_noSelf = set(permutations(uniqueNodes, r = 2))
                # Find intersection of above lists to ignore self edges
                # TODO: is there a better way of doing this?
                possibleEdges = possibleEdges_TF.intersection(possibleEdges_noSelf)
                TrueEdgeDict = {'|'.join(p):0 for p in possibleEdges}
                trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
                trueEdges = trueEdges[trueEdges.isin(TrueEdgeDict)]
                #print("\nEdges considered ", len(trueEdges))
                #print("\nPossible edges ", len(possibleEdges))
                numEdges = len(trueEdges)
                
                predDF['Edges'] = predDF['Gene1'] + "|" + predDF['Gene2']
                # limit the predicted edges to the genes that are in the ground truth
                predDF = predDF[predDF['Edges'].isin(TrueEdgeDict)]
            else:
                trueEdges = trueEdgesDF['Gene1'] + "|" + trueEdgesDF['Gene2']
                trueEdges = set(trueEdges.values)
                numEdges = len(trueEdges)
        
            # check if ranked edges list is empty
            # if so, it is just set to an empty set

            if not predDF.shape[0] == 0:

                # we want to ensure that we do not include
                # edges without any edge weight
                # so check if the non-zero minimum is
                # greater than the edge weight of the top-kth
                # node, else use the non-zero minimum value.
                predDF.EdgeWeight = predDF.EdgeWeight.round(6)
                predDF.EdgeWeight = predDF.EdgeWeight.abs()

                # Use num True edges or the number of
                # edges in the dataframe, which ever is lower
                maxk = min(predDF.shape[0], numEdges)

                if algo[0] == 'PPCOR':
                    edgeWeightTopk = predDF.iloc[2*maxk-1+25].EdgeWeight
                elif algo[0] == 'PIDC':
                    edgeWeightTopk = predDF.iloc[2*maxk-1].EdgeWeight
                else:
                    edgeWeightTopk = predDF.iloc[maxk-1].EdgeWeight

                nonZeroMin = np.nanmin(predDF.EdgeWeight.replace(0, np.nan).values)
                bestVal = max(nonZeroMin, edgeWeightTopk)

                if algo[0] == 'PPCOR':
                    newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal) & (predDF['EdgeWeight'] < 1)]
                else:
                    newDF = predDF.loc[(predDF['EdgeWeight'] >= bestVal)]


                predGraph = nx.DiGraph()


                for key in EdgeDict.keys():
                    u = key.split('|')[0]
                    v = key.split('|')[1]
                    if len(newDF.loc[(newDF['Gene1'] == u) &
                           (newDF['Gene2'] == v)])>0:
                            predGraph.add_edge(u,v)

                if algo[0] == 'PPCOR' or algo[0] == 'PIDC':
                    dataDict['PR'][algo[0]], dataDict['IC'][algo[0]], dataDict['OC'][algo[0]], dataDict['EIG'][algo[0]], dataDict['BT'][algo[0]], dataDict['RD'][algo[0]],hubsDict['PR'][algo[0]],a, hubsDict['IC'][algo[0]],b, hubsDict['OC'][algo[0]],c, hubsDict['BT'][algo[0]],d, hubsDict['RD'][algo[0]],e = Hubs_undirected(predGraph,refGraph,set(trueEdgesDF.Gene1),set(trueEdgesDF.Gene2)-set(RBPs))
                else:
                    dataDict['PR'][algo[0]], dataDict['IC'][algo[0]], dataDict['OC'][algo[0]], dataDict['EIG'][algo[0]], dataDict['BT'][algo[0]], dataDict['RD'][algo[0]],hubsDict['PR'][algo[0]],a, hubsDict['IC'][algo[0]],b, hubsDict['OC'][algo[0]],c, hubsDict['BT'][algo[0]],d, hubsDict['RD'][algo[0]],e = Hubs_directed(predGraph,refGraph,set(trueEdgesDF.Gene1),set(trueEdgesDF.Gene2)-set(RBPs))
            else:
                # no edges are predicted, set to 0!
                dataDict['PR'][algo[0]] = 0
                dataDict['IC'][algo[0]] = 0
                dataDict['OC'][algo[0]] = 0
                dataDict['EIG'][algo[0]] = 0
                dataDict['BT'][algo[0]] = 0
                dataDict['RD'][algo[0]] = 0
                hubsDict['PR'][algo[0]] = 0
                hubsDict['IC'][algo[0]] = 0
                hubsDict['OC'][algo[0]] = 0
                hubsDict['BT'][algo[0]] = 0
                hubsDict['RD'][algo[0]] = 0
        else:
            print(outDir + '/' +algo[0]+'/rankedEdges.csv', \
                  ' does not exist. Skipping...')

            dataDict['PR'][algo[0]] = 0
            dataDict['IC'][algo[0]] = 0
            dataDict['OC'][algo[0]] = 0
            dataDict['EIG'][algo[0]] = 0
            dataDict['BT'][algo[0]] = 0
            dataDict['RD'][algo[0]] = 0
            hubsDict['PR'][algo[0]] = 0
            hubsDict['IC'][algo[0]] = 0
            hubsDict['OC'][algo[0]] = 0
            hubsDict['BT'][algo[0]] = 0
            hubsDict['RD'][algo[0]] = 0

    dataDF = pd.DataFrame(dataDict)
    hubsDF = pd.DataFrame(hubsDict)
    print("success")
    return dataDF['PR'], dataDF['IC'], dataDF['OC'], dataDF['EIG'], dataDF['BT'], dataDF['RD'],hubsDF['PR'], hubsDF['IC'], hubsDF['OC'], hubsDF['BT'], hubsDF['RD'];


def Hubs_undirected(inGraph, refGraph,regs,tgts):
    # A helper function to compute
    # hub similarity for undirected output
    inGraphU = nx.DiGraph.to_undirected(inGraph)
    refGraphU = nx.DiGraph.to_undirected(refGraph)
    # pranks_1 = nx.pagerank(inGraphU, max_iter=500)
    # pranks_1_s = sorted(pranks_1.items(), key=lambda kv: kv[1],reverse=True)
    # pranks_1_s = [i for i in pranks_1_s if i[0] in list(regs)]
    # pranks_2 = nx.pagerank(refGraphU, max_iter=500)
    # pranks_2_s = sorted(pranks_2.items(), key=lambda kv: kv[1],reverse=True)
    # pranks_2_s = [i for i in pranks_2_s if i[0] in list(regs)]
    # prhubs1 = [a_tuple1[0] for a_tuple1 in pranks_1_s[0:int(0.1*len(pranks_2_s))]]
    # prhubs2 = [a_tuple2[0] for a_tuple2 in pranks_2_s[0:int(0.1*len(pranks_2_s))]]
    # pr = jaccard(prhubs1, prhubs2)

    deg_centr1 = nx.degree_centrality(inGraphU)
    deg_centr1_s = sorted(deg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr1_s = [i for i in deg_centr1_s if i[0] in list(tgts)]
    deg_centr2 = nx.degree_centrality(refGraphU)
    deg_centr2_s = sorted(deg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr2_s = [i for i in deg_centr2_s if i[0] in list(tgts)]
    idhubs1 = [d_tuple1[0] for d_tuple1 in list(deg_centr1_s)[0:int(0.1*len(deg_centr2_s))]]
    idhubs2 = [d_tuple2[0] for d_tuple2 in list(deg_centr2_s)[0:int(0.1*len(deg_centr2_s))]]
    icentr = jaccard(idhubs1, idhubs2)
    
    deg_centr1 = nx.degree_centrality(inGraphU)
    deg_centr1_s = sorted(deg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr1_s = [i for i in deg_centr1_s if i[0] in list(regs)]
    deg_centr2 = nx.degree_centrality(refGraphU)
    deg_centr2_s = sorted(deg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    deg_centr2_s = [i for i in deg_centr2_s if i[0] in list(regs)]
    odhubs1 = [d_tuple1[0] for d_tuple1 in list(deg_centr1_s)[0:int(0.1*len(deg_centr2_s))]]
    odhubs2 = [d_tuple2[0] for d_tuple2 in list(deg_centr2_s)[0:int(0.1*len(deg_centr2_s))]]
    ocentr = jaccard(odhubs1, odhubs2)

    # eig_centr1 = nx.eigenvector_centrality(inGraphU, max_iter=500)
    # eig_centr1_s = sorted(eig_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # eig_centr1_s = [i for i in eig_centr1_s if i[0] in list(regs)]
    # eig_centr2 = nx.eigenvector_centrality(refGraphU, max_iter=500)
    # eig_centr2_s = sorted(eig_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # eig_centr2_s = [i for i in eig_centr2_s if i[0] in list(regs)]
    # eighubs1 = [d_tuple1[0] for d_tuple1 in list(eig_centr1_s)[0:int(0.1*len(eig_centr2_s))]]
    # eighubs2 = [d_tuple2[0] for d_tuple2 in list(eig_centr2_s)[0:int(0.1*len(eig_centr2_s))]]
    # eig = jaccard(eighubs1, eighubs2)

    # bt_centr1 = nx.betweenness_centrality(inGraphU)
    # bt_centr1_s = sorted(bt_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # bt_centr1_s = [i for i in bt_centr1_s if i[0] in list(regs)]
    # bt_centr2 = nx.betweenness_centrality(refGraphU)
    # bt_centr2_s = sorted(bt_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # bt_centr2_s = [i for i in bt_centr2_s if i[0] in list(regs)]
    # bthubs1 = [d_tuple1[0] for d_tuple1 in list(bt_centr1_s)[0:int(0.1*len(bt_centr2_s))]]
    # bthubs2 = [d_tuple2[0] for d_tuple2 in list(bt_centr2_s)[0:int(0.1*len(bt_centr2_s))]]
    # bt = jaccard(bthubs1, bthubs2)

    # rad_centr1 = radiality(inGraph)
    # rad_centr1_s = sorted(rad_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # rad_centr1_s = [i for i in rad_centr1_s if i[0] in list(regs)]
    # rad_centr2 = radiality(refGraph)
    # rad_centr2_s = sorted(rad_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # rad_centr2_s = [i for i in rad_centr2_s if i[0] in list(regs)]
    # radhubs1 = [d_tuple1[0] for d_tuple1 in list(rad_centr1_s)[0:int(0.1*len(rad_centr2_s))]]
    # radhubs2 = [d_tuple2[0] for d_tuple2 in list(rad_centr2_s)[0:int(0.1*len(rad_centr2_s))]]
    # rad = jaccard(radhubs1, radhubs2)
    #rad = 0

    return icentr, icentr, ocentr, icentr, icentr, icentr,idhubs1,idhubs2,idhubs1,idhubs2,odhubs1,odhubs2,idhubs1,idhubs2,idhubs1,idhubs2;

def Hubs_directed(inGraph, refGraph,regs,tgts):
    # A helper function to compute
    # hub similarity for directed output
    # pranks_1 = nx.pagerank(inGraph, max_iter=500)
    # pranks_1_s = sorted(pranks_1.items(), key=lambda kv: kv[1],reverse=True)
    # pranks_1_s = [i for i in pranks_1_s if i[0] in list(regs)]
    # pranks_2 = nx.pagerank(refGraph, max_iter=500)
    # pranks_2_s = sorted(pranks_2.items(), key=lambda kv: kv[1],reverse=True)
    # pranks_2_s = [i for i in pranks_2_s if i[0] in list(regs)]
    # prhubs1 = [a_tuple1[0] for a_tuple1 in pranks_1_s[0:int(0.1*len(pranks_2_s))]]
    # prhubs2 = [a_tuple2[0] for a_tuple2 in pranks_2_s[0:int(0.1*len(pranks_2_s))]]
    # pr = jaccard(prhubs1, prhubs2)

    ideg_centr1 = nx.in_degree_centrality(inGraph)
    ideg_centr1_s = sorted(ideg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    ideg_centr1_s = [i for i in ideg_centr1_s if i[0] in list(tgts)]
    ideg_centr2 = nx.in_degree_centrality(refGraph)
    ideg_centr2_s = sorted(ideg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    ideg_centr2_s = [i for i in ideg_centr2_s if i[0] in list(tgts)]
    idhubs1 = [d_tuple1[0] for d_tuple1 in list(ideg_centr1_s)[0:int(0.1*len(ideg_centr2_s))]]
    idhubs2 = [d_tuple2[0] for d_tuple2 in list(ideg_centr2_s)[0:int(0.1*len(ideg_centr2_s))]]
    icentr = jaccard(idhubs1, idhubs2)

    odeg_centr1 = nx.out_degree_centrality(inGraph)
    odeg_centr1_s = sorted(odeg_centr1.items(), key=lambda kv: kv[1],reverse=True)
    odeg_centr1_s = [i for i in odeg_centr1_s if i[0] in list(regs)]    
    odeg_centr2 = nx.out_degree_centrality(refGraph)
    odeg_centr2_s = sorted(odeg_centr2.items(), key=lambda kv: kv[1],reverse=True)
    odeg_centr2_s = [i for i in odeg_centr2_s if i[0] in list(regs)]    
    odhubs1 = [d_tuple1[0] for d_tuple1 in list(odeg_centr1_s)[0:int(0.1*len(odeg_centr2_s))]]
    odhubs2 = [d_tuple2[0] for d_tuple2 in list(odeg_centr2_s)[0:int(0.1*len(odeg_centr2_s))]]
    ocentr = jaccard(odhubs1, odhubs2)

    # eig_centr1 = nx.eigenvector_centrality(inGraph)
    # eig_centr1_s = sorted(eig_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # eig_centr2 = nx.eigenvector_centrality(refGraph)
    # eig_centr2_s = sorted(eig_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # eighubs1 = [d_tuple1[0] for d_tuple1 in list(eig_centr1_s)[0:4]]
    # eighubs2 = [d_tuple2[0] for d_tuple2 in list(eig_centr2_s)[0:4]]
    # eig = jaccard(eighubs1, eighubs2)
    # eig = 0

    # bt_centr1 = nx.betweenness_centrality(inGraph)
    # bt_centr1_s = sorted(bt_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # bt_centr1_s = [i for i in bt_centr1_s if i[0] in list(regs)]        
    # bt_centr2 = nx.betweenness_centrality(refGraph)
    # bt_centr2_s = sorted(bt_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # bt_centr2_s = [i for i in bt_centr2_s if i[0] in list(regs)]            
    # bthubs1 = [d_tuple1[0] for d_tuple1 in list(bt_centr1_s)[0:int(0.1*len(bt_centr2_s))]]
    # bthubs2 = [d_tuple2[0] for d_tuple2 in list(bt_centr2_s)[0:int(0.1*len(bt_centr2_s))]]
    # bt = jaccard(bthubs1, bthubs2)

    # rad_centr1 = radiality(inGraph)
    # rad_centr1_s = sorted(rad_centr1.items(), key=lambda kv: kv[1],reverse=True)
    # rad_centr1_s = [i for i in rad_centr1_s if i[0] in list(regs)]        
    # rad_centr2 = radiality(refGraph)
    # rad_centr2_s = sorted(rad_centr2.items(), key=lambda kv: kv[1],reverse=True)
    # rad_centr2_s = [i for i in rad_centr2_s if i[0] in list(regs)]        
    # radhubs1 = [d_tuple1[0] for d_tuple1 in list(rad_centr1_s)[0:int(0.1*len(rad_centr2_s))]]
    # radhubs2 = [d_tuple2[0] for d_tuple2 in list(rad_centr2_s)[0:int(0.1*len(rad_centr2_s))]]
    # rad = jaccard(radhubs1, radhubs2)
    #rad = 0

    return icentr, icentr, ocentr, icentr, icentr, icentr,idhubs1,idhubs2,idhubs1,idhubs2,odhubs1,odhubs2,idhubs1,idhubs2,idhubs1,idhubs2;


### Helper functions

def kl_divergence(d1, d2):
    max_l=max(d1.shape[0],d2.shape[0])

    p = np.zeros(max_l)
    q = np.zeros(max_l)

    p[:d1.shape[0]] = d1
    q[:d2.shape[0]] = d2

    return np.sum(np.where(np.bitwise_and(p != 0,q != 0), p * np.log2(p / q), 0))

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

def radiality(Graph):
	fw = nx.floyd_warshall(Graph)
	results = {a:dict(b) for a,b in fw.items()}
	N = Graph.order()
	DistanceMatrix = np.zeros((N,N))
	for i in range(0,N):
		for j in range(0,N):
			key1 = list(results.keys())[i]
			key2 = list(results.keys())[j]
			# DistanceMatrix[i,j]=results[i][j]
			DistanceMatrix[i,j]=results[key1][key2]
	if nx.is_connected(Graph.to_undirected()):
		RD=nx.diameter(Graph.to_undirected())*np.ones((N,N))+1-DistanceMatrix
	else:
		#S = max(nx.connected_component_subgraphs(Graph.to_undirected()), key=len)
		S = max((Graph.to_undirected().subgraph(c).copy() for c in nx.connected_components(Graph.to_undirected())), key=len)
		RD=nx.diameter(S)*np.ones((N,N))+1-DistanceMatrix
	radialities=(RD.sum(axis=1)-RD.diagonal())/(len(RD.diagonal())-1)
	#rad_dict = dict((k,radialities[k]) for k in range(0,N))
	rad_dict=dict((n,radialities[k]) for (n,k) in zip(Graph.nodes(),range(0,N)))
	return rad_dict
