For the inference of GRNs from single-cell RNA-seq data, we chose the three top performing methods from BEELINE (PIDC, GRNBOOST2, SINCERITIES). We also added two more recent methods that have been shown to outperform previous ones (TENET and DeePSEM), and a method not specifically designed for single-cell RNA-seq, but that has good performance and it is widely used (ARACNe).  
GRN inference is based on [BEELINE](https://github.com/Murali-group/Beeline). 

PIDC, GRNBOOST2 and SINCERITIES were already available in BEELINE, thus we used the Docker containers provided in it. For TENET, we followed the instructions provided at [https://github.com/neocaleb/TENET](https://github.com/neocaleb/TENET). The script [InferencePipeline.py](InferencePipeline.py) is a wrapper to automatically generate the yaml configuration file used in BEELINE and adding TENET inference to the pipeline. 

DeePSEM instead runs on a GPU architecture and it was implemented following the instructions provided by the authors [https://github.com/HantaoShu/DeepSEM](https://github.com/HantaoShu/DeepSEM). In the folder [DeePSEM](/DeePSEM/) we provide the codes used to run GRN inference using DeePSEM on GPUs.


ARACNe-AP was run on a HPC cluster following the instructions provided by the authors [https://github.com/califano-lab/PISCES/tree/master/data](https://github.com/califano-lab/PISCES/tree/master/data). In the folder [ARACNe](/ARACNe/) we provide the codes used to run GRN inference using ARACNe on a HPC cluster.

The GRNs inferred by the different algorithms are then formatted as in BEELINE using our script [PostProcessing.py](PostProcessing.py), which can be run as

```
python PostProcessRankings.py --folder ./ --dataID HepG2_Smartseq2_RBP_RNA500
```
for the HepG2 Smart-seq2 dataset with RBPs and 500 HVGs.
