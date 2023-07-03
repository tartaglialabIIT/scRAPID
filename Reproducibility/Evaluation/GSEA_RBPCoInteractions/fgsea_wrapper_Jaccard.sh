#!/bin/bash

python prepare_rnk_files.py

Rscript fgsea.R

python process_fgsea_output.py
