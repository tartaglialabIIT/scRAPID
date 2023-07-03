Here we provide the codes to run ARACNE-AP on a HPC cluster, as described in [https://github.com/califano-lab/PISCES](https://github.com/califano-lab/PISCES).

We provide a def file to build the Singularity image for running ARACNe.

The R scripts [convert_h5_seurat.R](convert_h5_seurat.R) and [data_prep_HPC.R](data_prep_HPC.R) are used to format the input required by ARACNe. Starting from the processed data that we provide, rds objects are generated and the tsv files required by ARACNe are then generated ([input_tables](./input_tables/)).

We also provide the code to generate metacells for large datasets ([HepG2_DNBelab_METACELLS](./HepG2_DNBelab_METACELLS/)).

The folder [HPC_scripts](./HPC_scripts/) provides example bash scripts to run ARACNe on a HPC cluster.
