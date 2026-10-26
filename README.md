Building PANDA Regulatory Networks from GTEx Gene Expression Data in R

Lopes-Ramos et al. [https://pubmed.ncbi.nlm.nih.gov/28899340/] previously applied PANDA to expression data from GTEx to estimate four gene regulatory networks, including networks for two cell lines (lymphoblastoids and fibroblasts) together with their associated tissue-of-origin (whole blood and skin). Subsequent analyses of these networks identified important differences, including increased, but repressive, regulatory targeting of cell-cycle pathways in tissue compared to cell-line networks.

Code demonstrating how to model and analyze PANDA networks based on this analysis is available on netBooks (https://github.com/netZoo/netbooks). A copy of that code, together with the input data necessary to execute it, has been copied to this repository.

This repository contains two files and one directory:

(1) panda_gtex_tutorial_code.R -- R code demonstrating how to run PANDA and visualize the results when performing comparative network analysis

(2) panda_gtex_tutorial_visNetwork.html -- an html file containing the output from running the tutorial

(3) data -- a directory containing all the input files needed to run the code in panda_gtex_tutorial_code.R
