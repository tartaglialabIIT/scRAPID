#!/bin/sh
#PBS -l select=1:ncpus=24
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -N ARACNe_FULL
#PBS -q workq

cd /work/jfiorentino/ARACNe
chmod +x aracne_run_all_datasets.sh
singularity exec ./ARACNe-AP.sif bash ./aracne_run_all_datasets.sh
