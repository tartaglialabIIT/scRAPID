Here we provide the codes to run the evaluation of the inference performance. The [BLEval](./BLEval/) folder and the script [BLEvaluator_adapted.py](BLEvaluator_adapted.py) provide several extensions to the BEELINE and STREAMLINE frameworks, including:
* Computation of the early Precision restricted to RBP-RNA interactions
* Computation of hub RBPs and hub RNAs
* Computation of the fraction of true positives in the inferred ranked GRN
* Computation of RBP co-interactions based on shared RNA targets

The script [EvaluationPipeline.py](EvaluationPipeline.py) is a wrapper for running the evaluation using all the GRN inference algorithms, including the new functionalities that we provide. The bash script [RunEvaluation.sh](RunEvaluation.sh) runs the evaluation pipeline on different scRNA-seq datasets and gene sets.

The folder [Downsampling_Analysis](./Downsampling_Analysis/) contains the codes to perform the downsampling of eCLIP RBP-mRNA networks to match the statistics of RBP-lncRNA networks.

The folder [GSEA_RBPCoInteractions](./GSEA_RBPCoInteractions/) contains the codes to perform the Gene Set Enrichment Analysis for RBP co-interactions using the Bioplex Interactome database.
