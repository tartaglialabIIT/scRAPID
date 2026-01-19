# scRAPID
single-cell transcriptomic-based RnA-Protein Interaction Detection

[![DOI](https://zenodo.org/badge/653177334.svg)](https://zenodo.org/doi/10.5281/zenodo.10210488)

## 🚀 New!

A web-based implementation of **scRAPID** is now available at  
👉 **https://tools.tartaglialab.com/scrapid**

The web server provides an interactive interface to run the scRAPID pipeline without local installation, enabling protein–RNA interaction prediction directly from single-cell transcriptomic data.

To align the local scRAPID code to the web server functionalities, the Zenodo record now provides also the SQL databases with pre-computed catRAPID interaction propensities for 8 model organisms. See the updated tutorial for details.

--- 

This repository includes data and codes associated with the manuscript [**Fiorentino et al, Prediction of protein-RNA interactions from single-cell transcriptomic data, Nucleic Acids Research, 2024**](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkae076/7607879?searchresult=1). 

* The folder [Reproducibility](/Reproducibility/) contains the codes needed to reproduce the results shown in the manuscript, including data pre-processing, Gene Regulatory Network inference, evaluation of protein-RNA interaction inference and plotting the results.
* The folder [Tutorial](/Tutorial/) contains a guided tutorial to run the scRAPID pipeline on a new single-cell transcriptomic dataset.

Detailed explanations and instructions are provided in each folder.

If you use our pipeline please cite:

- Jonathan Fiorentino, Alexandros Armaos, Alessio Colantoni, Gian Gaetano Tartaglia, Prediction of protein-RNA interactions from single-cell transcriptomic data, Nucleic Acids Research, 2024, gkae076, https://doi.org/10.1093/nar/gkae076
- Jonathan Fiorentino*, Alexandros Armaos*, Chiara Montrone, Alessio Colantoni, Gian Gaetano Tartaglia,  scRAPID-web: a web server for predicting protein–RNA interactions from single-cell transcriptomics, bioRxiv, 2025, https://doi.org/10.1101/2025.03.12.642785
