Here use an example ranking and compute (without ground truth):

1) RBP co-interactions
2) Filter with catRAPID
3) RBP-lncRNA interactions
4) Hub RBPs and hub mRNAs and lncRNAs

These should be options provided to a script main.py

main.py takes:

- config file: algorithm for which the analysis should be run 
- optional: set of RBPs for RBP co-interactions (if this analysis is chosen)
- optional: catRAPID table file (if catRAPID filter is included in the analysis)

Then it 
- cleans the file with ranked edges
- runs function with RBP co-interactions

We can base the script on BEELINE evaluation, that takes the inferred ranking, the catRAPID file and the set of RBPs


Put in this folder example data:
- Ranking from an algorithm
- catRAPID table
- set of human and mouse RBPs