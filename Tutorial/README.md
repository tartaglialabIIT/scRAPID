Here we provide data and codes to run scRAPID on a GRN inferred from a scRNA-seq dataset, without knowing the ground truth network. A python or conda virtual environment can be installed based on the provided file [requirements.txt](requirements.txt).

We use as an example the network inferred by DeePSEM on a short-read SPLIT-seq dataset from the murine C2C12 cell line. We provide the input data in [inputs](./inputs/) and the GRNs inferred by each algorithm in [outputs](./outputs/).

The [tutorial](scRAPID_tutorial.ipynb) includes:
* Loading and pre-processing of a GRN inferred from single-cell transcriptomic data
* Prediction of RBP co-interactions based on shared RNA targets
* catRAPID-based filtering of RBP-RNA interactions
* Selection of RBP-lncRNA interactions
* Prediction of hub RBPs and hub RNAs (mRNAs and lncRNAs)
  
The script [scRAPID.py](scRAPID.py) contains the functions used in the tutorial.

We provide also:
* Lists of RBPs for human and mouse
* Lists of human and mouse long non-coding RNAs (obtained from GENCODE V41 and M30, respectively)
* The table with catRAPID interaction propensities for RBPs and RNAs contained in the dataset used in scRAPID tutorial.

To facilitate the usage of scRAPID with new scRNA-seq datasets in different organisms, we provide a SQL database containing the maximum interaction propensity scores from catRAPID for: 
* 3131 RBPs vs 62055 RNAs (all the canonical isoforms for the full transcriptome) in human, for a total of 194.3 millions interactions;
* 2900 RBPs and 53087 RNAs (all the canonical isoforms for the full transcriptome) in mouse, for a total of 154.0 millions interactions.

The lists of human and mouse RBPs were compiled by combining the RBPs from the RBP2GO database having score larger than 10 with those that make up the [catRAPID omics v2.0](http://service.tartaglialab.com/page/catrapid_omics2_group) RBP libraries; the latter sets were further expanded by including, for human and mouse, proteins that are orthologous to the RBPs identified in mouse and human, respectively.

The csv file can be obtained from our SQL database via curl. Example queries for a single RNA, a single RBP or a RBP-RNA pair:

```
curl "http://scrapid.tartaglialab.com/database?rna=Malat1&protein=Caprin1" > my_MOUSE_interactions.csv

curl "http://scrapid.tartaglialab.com/database?protein=FUS" > FUS_HUMAN_interactions.csv

curl "http://scrapid.tartaglialab.com/database?rna=XIST" > XIST_HUMAN_interactions.csv

```

We also provide a bash script ([query_multiple_pairs.sh](query_multiple_pairs.sh)) to perform queries for multiple RBP-RNA pairs stored in a text file (see [protein_rna_pairs.txt](protein_rna_pairs.txt)); the example output is in [concatenated_output.csv](concatenated_output.csv).
