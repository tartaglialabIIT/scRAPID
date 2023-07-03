
"""
BEELINE Evaluation (:mod:`BLEval`) module contains the following
:class:`BLEval.BLEval` and three additional classes used in the
definition of BLEval class

- :class:`BLEval.ConfigParser`
- :class:`BLEval.InputSettings`
- :class:`BLEval.OutputSettings`


"""
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


# local imports
from BLEval.parseTime import getTime
from BLEval.computeDGAUC import PRROC
from BLEval.computeBorda import Borda
from BLEval.computeJaccard import Jaccard
from BLEval.computeSpearman import Spearman
from BLEval.computeNetMotifs import Motifs
from BLEval.computeEarlyPrec import EarlyPrec
from BLEval.computeEarlyPrecTgt import EarlyPrecTgt
from BLEval.computePathStats import pathAnalysis
from BLEval.computeSignedEPrec import signedEPrec
from BLEval.computeTopology import Topologies
#from BLEval.computeExperimentalTopology import ExpTopologies
from BLEval.computeHubs import Hubs
from BLEval.computeCoInter import CoInter
from BLEval.computeHubsTargets import HubsTargets
#from BLEval.computeTopologyMatrixbased import AdjacencyMatrices
from BLEval.computeTopRank import RankTopK
from BLEval.computeTopRankTgt import RankTopKTgt

class InputSettings(object):
    '''
    The class for storing the names of input files.
    This initilizes an InputSettings object based on the
    following three parameters.


    :param datadir:   input dataset root directory, typically 'inputs/'
    :type datadir: str

    :param datasets:   List of dataset names
    :type datasets: list

    :param algorithms:   List of algorithm names
    :type algorithms: list
    '''

    def __init__(self,
            datadir, datasets, algorithms) -> None:

        self.datadir = datadir
        self.datasets = datasets
        self.algorithms = algorithms


class OutputSettings(object):
    '''
    The class for storing the names of directories that output should
    be written to. This initilizes an OutputSettings object based on the
    following two parameters.


    :param base_dir: output root directory, typically 'outputs/'
    :type base_dir: str
    :param output_prefix: A prefix added to the final output files.
    :type str:
    '''

    def __init__(self, base_dir, output_prefix: Path) -> None:
        self.base_dir = base_dir
        self.output_prefix = output_prefix




class BLEval(object):
    '''
    The BEELINE Evaluation object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs.
    '''

    def __init__(self,
            input_settings: InputSettings,
            output_settings: OutputSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings


    def computeAUC(self, directed = True):

        '''
        Computes areas under the precision-recall (PR) and
        and ROC plots for each algorithm-dataset combination.

        Parameters
        ----------
        directedFlag: bool
            A flag to specifiy whether to treat predictions
            as directed edges (directed = True) or
            undirected edges (directed = False).

        :returns:
            - AUPRC: A dataframe containing AUPRC values for each algorithm-dataset combination
            - AUROC: A dataframe containing AUROC values for each algorithm-dataset combination
        '''
        AUPRCDict = {}
        AUROCDict = {}

        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):

            AUPRC, AUROC = PRROC(dataset, self.input_settings,
                                    directed = directed, selfEdges = False, plotFlag = False)
            AUPRCDict[dataset['name']] = AUPRC
            AUROCDict[dataset['name']] = AUROC
        AUPRC = pd.DataFrame(AUPRCDict)
        AUROC = pd.DataFrame(AUROCDict)
        return AUPRC, AUROC


    def parseTime(self):
        """
        Parse time output for each
        algorithm-dataset combination.

        :returns:
            A dictionary of times for all dataset-algorithm combinations
        """
        TimeDict = dict()

        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            timevals  = getTime(self, dataset)
            TimeDict[dataset["name"]] = timevals

        return TimeDict

    def computeJaccard(self):

        '''
        Computes Jaccard Index between top-k edge predictions
        of the same algorithm.

        :returns:
            A dataframe containing the median and median absolute
            deviation of the Jaccard Index values of each algorithm
            on the given set of datasets.
        '''
        JaccDF = {}
        JaccDF['Jaccard Median'] = {}
        JaccDF['Jaccard MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                JaccDF['Jaccard Median'][algo[0]], JaccDF['Jaccard MAD'][algo[0]] = Jaccard(self, algo[0])

        return pd.DataFrame(JaccDF)


    def computeSpearman(self):

        '''
        Finds the Spearman's correlation coefficient between
        the ranked edges of the same algorithm on the given
        set of datasets.

        :returns:
            A dataframe containing the median and median absolute
            deviation of the Separman's correlation values of each algorithm.
        '''
        corrDF = {}
        corrDF['Spearman Median'] = {}
        corrDF['Spearman MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                corrDF['Spearman Median'][algo[0]],corrDF['Spearman MAD'][algo[0]] = Spearman(self, algo[0])

        return pd.DataFrame(corrDF)



    def computeNetMotifs(self):

        '''
        For each algorithm-dataset combination, this function computes the network motifs such as
        Feedforward loops, Feedback loops and
        Mutual interactions in the predicted top-k network. It returns the ratio of network motif counts
        compared to their respective values in the reference network.

        :returns:
            - FBL: A dataframe containing ratios of number of Feedback loops

            - FFL: A dataframe containing ratios of number of Feedforward loops

            - MI: A dataframe containing ratios of number of Mutual Interactions
        '''
        FFLDict = {}
        FBLDict = {}
        MIDict = {}
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            FBLDict[dataset["name"]], FFLDict[dataset["name"]], MIDict[dataset["name"]] = Motifs(dataset, self.input_settings)

        FBL = pd.DataFrame(FBLDict)
        FFL = pd.DataFrame(FFLDict)
        MI = pd.DataFrame(MIDict)

        return FBL, FFL, MI

    def computeTopology(self):

        '''
        For each algorithm-dataset combination, this function computes the network topologic properties

        :returns:ns
        '''
        # FFLDict = {}
        # FBLDict = {}
        # MIDict = {}
        print("Reached init.py")

        CCDict = {}
        GEFFDict = {}
        LEFFDict = {}
        NCDict = {}
        SPDict = {}
        ASSDict = {}
        iDSDict = {}
        oDSDict = {}
        iCTDict = {}
        oCTDict = {}
        HJDict = {}
        DJDict = {}
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            CCDict[dataset["name"]], GEFFDict[dataset["name"]], LEFFDict[dataset["name"]], NCDict[dataset["name"]], SPDict[dataset["name"]], ASSDict[dataset["name"]], iDSDict[dataset["name"]], oDSDict[dataset["name"]], iCTDict[dataset["name"]], oCTDict[dataset["name"]], HJDict[dataset["name"]], DJDict[dataset["name"]] = Topologies(dataset, self.input_settings)

        CC = pd.DataFrame(CCDict)
        GEFF = pd.DataFrame(GEFFDict)
        LEFF = pd.DataFrame(LEFFDict)
        NC = pd.DataFrame(NCDict)
        SP = pd.DataFrame(SPDict)
        ASS = pd.DataFrame(ASSDict)
        iDS = pd.DataFrame(iDSDict)
        oDS = pd.DataFrame(oDSDict)
        iCT = pd.DataFrame(iCTDict)
        oCT = pd.DataFrame(oCTDict)
        HJ = pd.DataFrame(HJDict)
        DJ = pd.DataFrame(DJDict)
        #
        # FBL = pd.DataFrame(FBLDict)
        # FFL = pd.DataFrame(FFLDict)
        # MI = pd.DataFrame(MIDict)
        print("Finished init.py")

        return CC, GEFF, LEFF, NC, SP, ASS, iDS, oDS, iCT, oCT, HJ, DJ

    def computeExperimentalTopology(self):

        '''
        For each algorithm-dataset combination, this function computes the network topologic properties

        :returns:ns
        '''
        # FFLDict = {}
        # FBLDict = {}
        # MIDict = {}
        print("Reached init.py")

        CCDict = {}
        GEFFDict = {}
        LEFFDict = {}
        NCDict = {}
        SPDict = {}
        ASSDict = {}
        CTDict = {}
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            CCDict[dataset["name"]], GEFFDict[dataset["name"]], LEFFDict[dataset["name"]], NCDict[dataset["name"]], SPDict[dataset["name"]], ASSDict[dataset["name"]], DSDict[dataset["name"]], CTDict[dataset["name"]], HJDict[dataset["name"]], DJDict[dataset["name"]] = ExpTopologies(dataset, self.input_settings)

        CC = pd.DataFrame(CCDict)
        GEFF = pd.DataFrame(GEFFDict)
        LEFF = pd.DataFrame(LEFFDict)
        NC = pd.DataFrame(NCDict)
        SP = pd.DataFrame(SPDict)
        ASS = pd.DataFrame(ASSDict)
        CT = pd.DataFrame(CTDict)
        print("Finished init.py")

        return CC, GEFF, LEFF, NC, SP, ASS, CT

    def computeHubs(self):

        '''
        For each algorithm-dataset combination, this function computes the hub regulators

        :returns:ns
        '''
        # FFLDict = {}
        # FBLDict = {}
        # MIDict = {}
        print("Reached init.py")

        PRict = {}
        ICDict = {}
        OCDict = {}
        EIGDict = {}
        BTDict = {}
        RDDict = {}
        PRHubs={}
        ICHubs={}
        OCHubs={}
        BTHubs={}
        RDHubs={}
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            PRict[dataset["name"]], ICDict[dataset["name"]], OCDict[dataset["name"]], EIGDict[dataset["name"]], BTDict[dataset["name"]], RDDict[dataset["name"]],PRHubs[dataset["name"]], ICHubs[dataset["name"]], OCHubs[dataset["name"]], BTHubs[dataset["name"]], RDHubs[dataset["name"]]  = Hubs(dataset, self.input_settings)

        PR = pd.DataFrame(PRict)
        IC = pd.DataFrame(ICDict)
        OC = pd.DataFrame(OCDict)
        EIG = pd.DataFrame(EIGDict)
        BT = pd.DataFrame(BTDict)
        RD = pd.DataFrame(RDDict)
        PRHubsDF = pd.DataFrame(PRHubs)
        ICHubsDF = pd.DataFrame(ICHubs)
        OCHubsDF = pd.DataFrame(OCHubs)
        BTHubsDF = pd.DataFrame(BTHubs)
        RDHubsDF = pd.DataFrame(RDHubs)
        print("Finished init.py")

        return PR, IC, OC, EIG, BT, RD,PRHubsDF,ICHubsDF,OCHubsDF,BTHubsDF,RDHubsDF;


    def computeCoInter(self):

        '''
        For each algorithm-dataset combination, this function computes the RBP co-interactions

        :returns:ns
        '''
        
        print("Reached init.py")

        PRict = {}
    
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            PRict[dataset["name"]]  = Complexes(dataset, self.input_settings)

        PR = pd.DataFrame(PRict)
        print("Finished init.py")

        return PR;

        
    def computeHubsTargets(self):

        '''
        For each algorithm-dataset combination, this function computes the hub targets

        :returns:ns
        '''
        # FFLDict = {}
        # FBLDict = {}
        # MIDict = {}
        print("Reached init.py")

        PRict = {}
        ICDict = {}
        OCDict = {}
        EIGDict = {}
        BTDict = {}
        RDDict = {}
        PRHubs={}
        ICHubs={}
        OCHubs={}
        BTHubs={}
        RDHubs={}
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            PRict[dataset["name"]], ICDict[dataset["name"]], OCDict[dataset["name"]], EIGDict[dataset["name"]], BTDict[dataset["name"]], RDDict[dataset["name"]],PRHubs[dataset["name"]], ICHubs[dataset["name"]], OCHubs[dataset["name"]], BTHubs[dataset["name"]], RDHubs[dataset["name"]]  = HubsTargets(dataset, self.input_settings)

        PR = pd.DataFrame(PRict)
        IC = pd.DataFrame(ICDict)
        OC = pd.DataFrame(OCDict)
        EIG = pd.DataFrame(EIGDict)
        BT = pd.DataFrame(BTDict)
        RD = pd.DataFrame(RDDict)
        PRHubsDF = pd.DataFrame(PRHubs)
        ICHubsDF = pd.DataFrame(ICHubs)
        OCHubsDF = pd.DataFrame(OCHubs)
        BTHubsDF = pd.DataFrame(BTHubs)
        RDHubsDF = pd.DataFrame(RDHubs)
        print("Finished init.py")

        return PR, IC, OC, EIG, BT, RD,PRHubsDF,ICHubsDF,OCHubsDF,BTHubsDF,RDHubsDF;


    def computeTopologyMatrixbased(self):

        '''
        For each algorithm-dataset combination, this function computes the network topologic properties

        :returns:ns
        '''
        print("Reached init.py")

        RDict = {}
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            RDict[dataset["name"]] = AdjacencyMatrices(dataset, self.input_settings)

        R = pd.DataFrame(RDict)
        print("Finished init.py")

        return R



    def computePaths(self):
        '''
        For each algorithm-dataset combination, this function computes path lengths
        through TP edges and FP edges, returns statistics on path lengths.

        :returns:
            - pathStats: A dataframe path lengths in predicted network
        '''
        for dataset in tqdm(self.input_settings.datasets,
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            pathAnalysis(dataset, self.input_settings, self.output_settings)


    def computeEarlyPrec(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values of the
        network formed using the predicted top-k edges.

        :returns:
            A dataframe containing the early precision values
            for each algorithm-dataset combination.

        '''
        Eprec = {}
        FracEdge = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                Eprec[algo[0]],FracEdge[algo[0]] = EarlyPrec(self, algo[0])
        return pd.DataFrame(Eprec).T,pd.DataFrame(FracEdge).T;
        
    
    def computeEarlyPrecTgt(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values of the
        network formed using the predicted top-k edges.

        :returns:
            A dataframe containing the early precision values
            for each algorithm-dataset combination.

        '''
        Eprectgt = {}
        FracEdgetgt = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                Eprectgt[algo[0]],FracEdgetgt[algo[0]] = EarlyPrecTgt(self, algo[0])
        return pd.DataFrame(Eprectgt).T,pd.DataFrame(FracEdgetgt).T;
                
    def computeTopRank(self):

        '''
        For each algorithm-dataset combination,
        this function computes the true positive in the top k edges.

        :returns:
            A dataframe containing the early precision values
            for each algorithm-dataset combination.

        '''
        rank = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                rank[algo[0]] = RankTopK(self, algo[0])
        return pd.DataFrame(rank).T

    def computeTopRankTgt(self):

        '''
        For each algorithm-dataset combination,
        this function computes the true positive in the top k edges.

        :returns:
            A dataframe containing the early precision values
            for each algorithm-dataset combination.

        '''
        rank = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                rank[algo[0]] = RankTopKTgt(self, algo[0])
        return pd.DataFrame(rank).T
        
    def computeSignedEPrec(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values separately
        for the activation and inhibitory edges.

        :returns:
            - A dataframe containing early precision for activation edges
            - A dataframe containing early precision for inhibitory edges
        '''
        sEPRDict = {}
        sEPRDict['EPrec Activation'] = {}
        sEPRDict['EPrec Inhibition'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                sEPrecDF = signedEPrec(self, algo[0])
                sEPRDict['EPrec Activation'][algo[0]] = sEPrecDF['+']
                sEPRDict['EPrec Inhibition'][algo[0]] = sEPrecDF['-']
        return(pd.DataFrame(sEPRDict['EPrec Activation']).T, pd.DataFrame(sEPRDict['EPrec Inhibition']).T)

    def computeBorda(self, selectedAlgorithms=None, aggregationMethod="average"):

        '''
        Computes edge ranked list using the Borda method for each dataset.
        Parameters
        ----------
        selectedAlgorithms: [str]
          List of algorithm names used to run borda method on selected
          algorithms. If nothing is provided, the function runs borda on
          all available algorithms.
        aggregationMethod: str
          Method used to aggregate rank in borda method. Available options are
          {‘average’, ‘min’, ‘max’, ‘first’}, default ‘average’
        :returns:
            None
        '''
        feasibleAlgorithmOptions = [algorithmName for algorithmName, _ in self.input_settings.algorithms]
        feasibleaggregationMethodOptions = ['average', 'min', 'max', 'first']

        selectedAlgorithms = feasibleAlgorithmOptions if selectedAlgorithms is None else selectedAlgorithms

        for a in selectedAlgorithms:
            if a not in feasibleAlgorithmOptions:
                print("\nERROR: No data available on algorithm %s. Please choose an algorithm from the following options: %s" % (a, feasibleAlgorithmOptions))
                return

        if aggregationMethod not in feasibleaggregationMethodOptions:
                print("\nERROR: Please choose an aggregation method algorithm from following options: " % feasibleaggregationMethodOptions)
                return

        Borda(self, selectedAlgorithms, aggregationMethod)

class ConfigParser(object):
    '''
    The class define static methods for parsing and storing the contents
    of the config file that sets a that sets a large number of parameters
    used in the BLEval.
    '''
    @staticmethod
    def parse(config_file_handle) -> BLEval:
        '''
        A method for parsing the input .yaml file.

        :param config_file_handle: Name of the .yaml file to be parsed
        :type config_file_handle: str

        :returns:
            An object of class :class:`BLEval.BLEval`.

        '''
        config_map = yaml.load(config_file_handle)
        return BLEval(
            ConfigParser.__parse_input_settings(
                config_map['input_settings']),
            ConfigParser.__parse_output_settings(
                config_map['output_settings']))

    @staticmethod
    def __parse_input_settings(input_settings_map) -> InputSettings:
        '''
        A method for parsing and initializing
        InputSettings object.
        '''
        input_dir = input_settings_map['input_dir']
        dataset_dir = input_settings_map['dataset_dir']
        datasets = input_settings_map['datasets']

        return InputSettings(
                Path(input_dir, dataset_dir),
                datasets,
                ConfigParser.__parse_algorithms(
                input_settings_map['algorithms']))


    @staticmethod
    def __parse_algorithms(algorithms_list):
        '''
        A method for parsing the list of algorithms
        that are being evaluated, along with
        any parameters being passed.

        Note that these parameters may not be
        used in the current evaluation, but can
        be used at a later point.
        '''

        # Initilalize the list of algorithms
        algorithms = []

        # Parse contents of algorithms_list
        for algorithm in algorithms_list:
                combos = [dict(zip(algorithm['params'], val))
                    for val in itertools.product(
                        *(algorithm['params'][param]
                            for param in algorithm['params']))]
                for combo in combos:
                    algorithms.append([algorithm['name'],combo])


        return algorithms

    @staticmethod
    def __parse_output_settings(output_settings_map) -> OutputSettings:
        '''
        A method for parsing and initializing
        Output object.
        '''
        output_dir = Path(output_settings_map['output_dir'])
        output_prefix = Path(output_settings_map['output_prefix'])

        return OutputSettings(output_dir,
                             output_prefix)
