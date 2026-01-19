# scRAPID tutorial

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
* Lists of RBPs for 8 model organisms (Arabidopsis Thaliana, Caenorhabditis Elegans, Danio Rerio, Drosophila Melanogaster, Homo sapiens, Mus musculus, Rattus norvegicus and Xenopus tropicalis)
* RNA libraries should be downloaded from the [Zenodo record](https://zenodo.org/records/18155494) and uncompressed.

# SQL databases for multiple model organisms

To facilitate the use of scRAPID with new scRNA-seq datasets and across different species, we now provide SQL databases for 8 model organisms, containing pre-computed catRAPID interaction propensities, distributed via the associated Zenodo record: [https://zenodo.org/records/18155494](https://zenodo.org/records/18155494). These databases can be locally downloaded and queried for the RBP–RNA interactions of interest (see the tutorial notebook), and directly used to reproduce the scRAPID tutorial on organism-specific datasets.

For each organism, the SQL database contains the maximum catRAPID interaction propensity scores between RBPs and RNAs, considering canonical transcript isoforms for the full transcriptome.
