

First install a conda virtual environment using the [requirements file](requirements.txt) provided in this folder.

* The folder [PreProcessing](/PreProcessing/) contains the Jupyter Notebooks for the pre-processing and gene selection of all the single-cell transcriptomic dataset used in our manuscript. We also provide the processed data for the HepG2 and K562 single-cell RNA-seq datasets to help the reproducibility of the inference results.

* The folder [Inference](/Inference/) contains the codes to run the Gene Regulatory Network inference using the algorithms PIDC, GRNBOOST2, SINCERITIES, TENET, DeePSEM and ARACNe.

* The folder [Evaluation](/Evaluation/) contains the codes associated to our additions to the [BEELINE](https://github.com/Murali-group/Beeline) and [STREAMLINE](https://github.com/ScialdoneLab/STREAMLINE) frameworks.

* The folder [Plotting](/Plotting/) contains the codes to reproduce the plots in the manuscript.

* The folder [Ground_Truth_Networks](/Ground_Truth_Networks/) contains the processed eCLIP, ChIP-seq and shRNA RNA-seq data from the ENCODE project used to evalute the inference performance for RBP-target and TF-target interactions for the HepG2 and K562 cell lines.

* The folder [catRAPID_filter](/catRAPID_filter/) contains the codes for the catRAPID-based filter of RBP-RNA interactions, and instructions for downloading catRAPID interaction propensity tables for human and mouse from our database.
