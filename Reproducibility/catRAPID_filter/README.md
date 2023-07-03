Here we provide the python script used to perform the catRAPID-based filter of inferred GRN and a bash script that runs it with varying thresholds. 

This analysis requires a csv file containing the catRAPID interaction propensities for RBP-RNA pairs potentially included in the GRN inferred from a scRNA-seq dataset. See an example table [here](/scRAPID/Tutorial/catRAPID_table_C2C12SR9KMBSC.csv).

To facilitate the usage of scRAPID with new scRNA-seq datasets in different organisms, we provide a SQL database containing the maximum interaction propensity scores from catRAPID for: 
* 3131 RBPs vs 62055 RNAs (all the canonical isoforms for the full transcriptome) in human, for a total of 194.3 millions interactions;
* 2900 RBPs and 53087 RNAs (all the canonical isoforms for the full transcriptome) in mouse, for a total of 154.0 millions interactions.

The lists of human and mouse RBPs were compiled by combining the RBPs from the RBP2GO database having score larger than 10 with those that make up the [catRAPID omics v2.0](http://service.tartaglialab.com/page/catrapid_omics2_group) RBP libraries; the latter sets were further expanded by including, for human and mouse, proteins that are orthologous to the RBPs identified in mouse and human, respectively.

We make the RBP lists available [here](scRAPID/Tutorial/).

Conversion between protein names and Uniprot IDs can be done using [Uniprot](https://www.uniprot.org/id-mapping); tables obtained from Uniprot for human and mouse are also provided [here](scRAPID/Tutorial/)
The csv file can be obtained from our SQL database via curl. Example queries for a single RNA, a single RBP or RBP-RNA pairs:

CHANGE address!!
```
http://127.0.0.1:5000/database?protein=P35637&rna=GAPDH" > my_interactions.csv

http://127.0.0.1:5000/database?rna=GAPDH" > my_interactions.csv

http://127.0.0.1:5000/database?protein=P35637" > my_interactions.csv

```

