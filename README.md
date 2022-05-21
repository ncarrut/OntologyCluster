# ontologyCluster

## Motivation  

Gene ontology (GO) enrichment analysis results are plagued by redundant categories bloating the results and making them hard to interpret.  It is desirable to group redundant categories together and display simplified results.  

## Methodology  

Various approaches are implemented for this task, here I calculate the overlap of protein membership between selected categories.  I use the overlap matrix to cluster the categories and finally summarize them according to the most strongly enriched category in a cluster.  This approach is particularity fit for proteomics analyses where coverage of categories can be limited.  It only uses identified proteins to group categories.  

## Implementation  
This is implemented in a R-package composed of 3 functions, (1) findCategoryOverlap.R, (2) clusterByOverlap.R and (3) summarizeCategoryCluster.R.  This package is not on CRAN etc.  You have to install it from GitHub.  
