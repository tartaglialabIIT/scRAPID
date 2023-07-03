#!/usr/bin/env python
# coding: utf-8

import os
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency
import pickle

# local imports
import BLEval as ev

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('-c','--config', default='config.yaml',
        help="Configuration file containing list of datasets "
              "algorithms and output specifications.\n")

    parser.add_argument('-a', '--auc', action="store_true", default=False,
        help="Compute median of areas under Precision-Recall and ROC curves.\n")

    parser.add_argument('-j', '--jaccard', action="store_true", default=False,
      help="Compute median Jaccard index of predicted top-k networks "
      "for each algorithm for a given set of datasets generated "
      "from the same ground truth network.\n")

    parser.add_argument('-r', '--spearman', action="store_true", default=False,
      help="Compute median Spearman Corr. of predicted edges "
      "for each algorithm  for a given set of datasets generated "
      " from the same ground truth network.\n")

    parser.add_argument('-t', '--time', action="store_true", default=False,
      help="Analyze time taken by each algorithm for a.\n")

    parser.add_argument('-e', '--epr', action="store_true", default=False,
      help="Compute median early precision.")

    parser.add_argument('-et', '--eprtgt', action="store_true", default=False,
      help="Compute median early precision.")
            
    parser.add_argument('-s','--sepr', action="store_true", default=False,
      help="Analyze median (signed) early precision for activation and inhibitory edges.")

    parser.add_argument('-m','--motifs', action="store_true", default=False,
      help="Compute network motifs in the predicted top-k networks.")

    parser.add_argument('-p','--paths', action="store_true", default=False,
      help="Compute path length statistics on the predicted top-k networks.")

    parser.add_argument('-b','--borda', action="store_true", default=False,
      help="Compute edge ranked list using the various Borda aggregatio methods.")

    parser.add_argument('-top','--topology', action="store_true", default=False,
      help="Analyse the network topology.")

    parser.add_argument('-ExpTop','--experimental_topology', action="store_true", default=False,
      help="Analyse the experimental network topology.")

    parser.add_argument('-hu','--hubs', action="store_true", default=False,
      help="Analyse the hub detection.")
      
    parser.add_argument('-hut','--hubstargets', action="store_true", default=False,
      help="Analyse the hub detection.")
    
    parser.add_argument('-co','--cointer', action="store_true", default=False,
      help="Compute protein co-interactions based on shared RNA targets.")
    
    parser.add_argument('-mat','--matrix', action="store_true", default=False,
      help="Compute adjacency matrices")
      
    parser.add_argument('-topk','--topk_ranking', action="store_true", default=False,
      help="Compute ranking of true positive fraction")
    
    parser.add_argument('-topktgt','--topk_ranking_tgt', action="store_true", default=False,
      help="Compute ranking of true positive fraction")

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    opts = parse_arguments()
    config_file = opts.config

    evalConfig = None

    with open(config_file, 'r') as conf:
        evalConfig = ev.ConfigParser.parse(conf)

    print('\nPost-run evaluation started...')
    evalSummarizer = ev.BLEval(evalConfig.input_settings, evalConfig.output_settings)
    print(str(evalSummarizer.output_settings.base_dir))
    print(str(evalSummarizer.input_settings.datadir))
    print(str(evalSummarizer.output_settings.output_prefix))
    outDir = str(evalSummarizer.output_settings.base_dir) + \
            str(evalSummarizer.input_settings.datadir).split("inputs")[1] + "/"+\
            str(evalSummarizer.output_settings.output_prefix) + "/"
    
    if os.path.isdir(outDir)==False:
        os.mkdir(outDir)

    # Compute and plot ROC, PRC and report median AUROC, AUPRC
    if (opts.auc):
        print('\n\nComputing areas under ROC and PR curves...')

        AUPRC, AUROC = evalSummarizer.computeAUC()
        AUPRC.to_csv(outDir+'AUPRC.csv')
        AUROC.to_csv(outDir+'AUROC.csv')

    # Compute Jaccard index
    if (opts.jaccard):
        print('\n\nComputing Jaccard index...')

        jaccDict = evalSummarizer.computeJaccard()
        jaccDict.to_csv(outDir + "Jaccard.csv")

    # Compute Spearman correlation scores
    if (opts.spearman):
        print('\n\nComputing Spearman\'s correlation...')

        corrDict = evalSummarizer.computeSpearman()
        corrDict.to_csv(outDir + "Spearman.csv")

    # Compute median time taken
    if (opts.time):
        print('\n\nComputing time taken...')

        TimeDict = evalSummarizer.parseTime()
        pd.DataFrame(TimeDict).to_csv(outDir+'Times.csv')

    # Compute early precision
    if (opts.epr):
        print('\n\nComputing early precision values...')
        ePRDF,FracDF = evalSummarizer.computeEarlyPrec()
        # REMOVE THIS FOR STANDARD ANALYSIS
        #Borda=True
        Borda=False
        if Borda:
            ePRDF.to_csv(outDir + "BordaEPr.csv")
        else:
            ePRDF.to_csv(outDir + "EPr.csv")
            FracDF.to_csv(outDir + "FracLeftEdges.csv")

    if (opts.eprtgt):
        print('\n\nComputing early precision values...')
        ePRDFTgt,FracDFTgt = evalSummarizer.computeEarlyPrecTgt()
        ePRDFTgt.to_csv(outDir + "EPrTgt.csv")
        FracDFTgt.to_csv(outDir + "FracLeftEdgesTgt.csv")    
                
    # Compute early precision for activation and inhibitory edges
    if (opts.sepr):
        print('\n\nComputing early precision values for activation and inhibitory edges...')

        actDF, inhDF = evalSummarizer.computeSignedEPrec()
        actDF.to_csv(outDir + "EPr-Activation.csv")
        inhDF.to_csv(outDir + "EPr-Inhibitory.csv")

    # Compute median time taken
    if (opts.motifs):
        print('\n\nComputing network motifs...')

        FBL, FFL, MI = evalSummarizer.computeNetMotifs()
        FBL.to_csv(outDir+'NetworkMotifs-FBL.csv')
        FFL.to_csv(outDir+'NetworkMotifs-FFL.csv')
        MI.to_csv(outDir+'NetworkMotifs-MI.csv')

    # Compute path statistics such as number of TP, FP,
    # and path lengths among TP in the top-k networks.
    if (opts.paths):
        print('\n\nComputing path length statistics on predicted networks...')
        evalSummarizer.computePaths()

   # Compute edge ranked list using the borda method
    if (opts.borda):
        print('\n\nComputing edge ranked list using the borda method')
        evalSummarizer.computeBorda()

    if (opts.topology):
        print('\n\nPerforming topological analysis')
        CC, GEFF, LEFF, NC, SP, ASS, iDS, oDS, iCT, oCT, HJ, DJ = evalSummarizer.computeTopology()
        CC.to_csv(outDir+'ClusteringCoefficient.csv')
        GEFF.to_csv(outDir+'GlobalEfficiency.csv')
        LEFF.to_csv(outDir+'LocalEfficiency.csv')
        NC.to_csv(outDir+'ConnectedComponents.csv')
        SP.to_csv(outDir+'ShortestPathLength.csv')
        ASS.to_csv(outDir+'Assortativity.csv')
        iDS.to_csv(outDir+'InDegreeSimilarity.csv')
        oDS.to_csv(outDir+'OutDegreeSimilarity.csv')
        iCT.to_csv(outDir+'InCentralization.csv')
        oCT.to_csv(outDir+'OutCentralization.csv')
        HJ.to_csv(outDir+'HubSimilarity.csv')
        DJ.to_csv(outDir+'DegreeJaccard.csv')

    if (opts.experimental_topology):
        print('\n\nPerforming topological analysis on experimental data')
        CC, GEFF, LEFF, NC, SP, ASS, CT = evalSummarizer.computeExperimentalTopology()
        CC.to_csv(outDir+'ClusteringCoefficient.csv')
        GEFF.to_csv(outDir+'GlobalEfficiency.csv')
        LEFF.to_csv(outDir+'LocalEfficiency.csv')
        NC.to_csv(outDir+'ConnectedComponents.csv')
        SP.to_csv(outDir+'ShortestPathLength.csv')
        ASS.to_csv(outDir+'Assortativity.csv')
        CT.to_csv(outDir+'Centralization.csv')

    if (opts.hubs):
        print('\n\nPerforming topological analysis on experimental data')
        PR, IC, OC, EIG, BT, RD, PRHubs, ICHubs, OCHubs, BTHubs, RDHubs = evalSummarizer.computeHubs()
        PR.to_csv(outDir+'PageRank.csv')
        IC.to_csv(outDir+'InCentrality.csv')
        OC.to_csv(outDir+'OutCentrality.csv')
        EIG.to_csv(outDir+'Eigenvector.csv')
        BT.to_csv(outDir+'Betweenness.csv')
        RD.to_csv(outDir+'Radiality.csv')
        PRHubs.to_csv(outDir+'PageRankHubs.csv')
        ICHubs.to_csv(outDir+'InCentralityHubs.csv')
        OCHubs.to_csv(outDir+'OutCentralityHubs.csv')
        BTHubs.to_csv(outDir+'BetweennessHubs.csv')
        RDHubs.to_csv(outDir+'RadialityHubs.csv')
        
    if (opts.hubstargets):
        print('\n\nPerforming topological analysis on experimental data')
        PR, IC, OC, EIG, BT, RD, PRHubs, ICHubs, OCHubs, BTHubs, RDHubs = evalSummarizer.computeHubsTargets()
        PR.to_csv(outDir+'PageRank_T.csv')
        IC.to_csv(outDir+'InCentrality_T.csv')
        OC.to_csv(outDir+'OutCentrality_T.csv')
        EIG.to_csv(outDir+'Eigenvector_T.csv')
        BT.to_csv(outDir+'Betweenness_T.csv')
        RD.to_csv(outDir+'Radiality_T.csv')
        PRHubs.to_csv(outDir+'PageRankHubs_T.csv')
        ICHubs.to_csv(outDir+'InCentralityHubs_T.csv')
        OCHubs.to_csv(outDir+'OutCentralityHubs_T.csv')
        BTHubs.to_csv(outDir+'BetweennessHubs_T.csv')
        RDHubs.to_csv(outDir+'RadialityHubs_T.csv')
    
    if (opts.cointer):
        print('\n\nPerforming topological analysis on experimental data')
        PR = evalSummarizer.computeCoInter()
        PR.to_csv(outDir+'CoInterFull.csv')

    if (opts.matrix):
        print('\n\nCalculating adjacency matrices')
        R = evalSummarizer.computeTopologyMatrixbased()
        R.to_csv(outDir+'AdjacencyMatrixStatus.csv')
        
    if (opts.topk_ranking):
        print('\n\nCalculating rankings of true positives')
        top = evalSummarizer.computeTopRank()
        with open(outDir+'Ranking.pickle', 'wb') as handle:
            pickle.dump(top, handle, protocol=pickle.HIGHEST_PROTOCOL)

    if (opts.topk_ranking_tgt):
        print('\n\nCalculating rankings of true positives')
        top = evalSummarizer.computeTopRankTgt()
        with open(outDir+'RankingTgt.pickle', 'wb') as handle:
            pickle.dump(top, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('\n\nEvaluation complete...\n')


if __name__ == '__main__':
  main()
