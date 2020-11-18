findCategoryOverlap <- function(membershipList, similarityCut = 0.5){
  ## get overlap matrix
  OG<-mat.or.vec(nr=length(membershipList), nc=length(membershipList))

  namesMap <- make.names(names(membershipList))
  names(namesMap) <- names(membershipList)

  rownames(OG) <- make.names(names(membershipList)); colnames(OG) <- make.names(names(membershipList))
  for (i in 1:nrow(OG)) for(j in 1:ncol(OG)){
    OG[i,j]<-length(intersect(membershipList[[i]], membershipList[[j]]))/
      max(c(length(membershipList[[i]]), length(membershipList[[j]])))
  }
  OG
}

clusterByOverlap <- function(OG, similarityCut = 0.5){
  ## Default to only cluster points with >50% overlap.
  ## make a diagnostic plot to see if that makes sense.
  ## look for a plateau around similarityCut followed by rapid drop

  require(kernlab)
  require(Spectrum)

  cutVec <- seq(0, 1, 0.05)
  cutResVec <- c()
  for(i in 1:length(cutVec)){
    cutResVec[i] <- length(unique(cutree(tree = hclust(as.dist(1-OG), method = "ward.D"), h = cutVec[i])))
  }
  diagnosticPlot <- data.frame(cutVec, cutResVec) %>%
    ggplot(aes(x = cutVec, y = cutResVec)) +
    geom_point() +
    labs(x = "Minimum similarity", y = "Clusters")

  ## decide which categories should be clustered
  doCluster <- apply(OG - diag(nrow = nrow(OG), ncol = ncol(OG)), 1, max) > similarityCut

  ## cluster using spectral clustering
  ## get K from Spectrum package
  myK <- Spectrum::estimate_k(OG[doCluster, doCluster], maxk = sum(doCluster)-1, showplots = FALSE) %>%
    filter(K > sqrt(sum(doCluster))) %>%
    filter(z == max(z)) %>%
    pull(K)

  ## then assign clusters from specc
  ## except where similarity was too small to cluster then
  ## just fill in with a sequence
  cluster1 <- specc(OG[doCluster, doCluster], centers = myK)@.Data
  cluster <- rep(0, nrow(OG))
  names(cluster) = rownames(OG)
  cluster[doCluster] <- cluster1
  cluster[!doCluster] <- seq(from = myK+1, to = myK+sum(cluster==0), by = 1)

  ## order according to euclidean distance.  This is for plotting and
  ## will be returned to help with evaluating other clustering methods if needed
  ord <- hclust(as.dist(1-OG), method = "ward.D")$order

  #plot to visualize clusters and assess clustering
  tilePlot <- data.frame(OG) %>%
    mutate(myRowNames = ordered(row.names(OG), levels = rownames(OG)[ord])) %>%
    mutate(cluster = cluster) %>%
    pivot_longer(-c(myRowNames, cluster), names_to = "myColNames") %>%
    #mutate(myColNames = gsub("\\.", ":", myColNames)) %>% ## because data.frame replaces ":"
    mutate(myColNames = ordered(myColNames, levels = rownames(OG)[ord])) %>%
    ggplot(aes(x = myRowNames, y = myColNames, fill = value)) +
    geom_tile() +
    theme(axis.text.y = element_text(size=6)) +
    scale_fill_gradient(name = "fraction",
                        low = "#FFFFFF",
                        high = "#012345") +
    labs(x = "cluster", y = "pathway") +
    scale_x_discrete(labels = cluster[ord]) +
    scale_y_discrete(label = function(x) strtrim(x, 35))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=6))

  return(list(overlapMatrix = OG, plot = tilePlot,
              plotOrder = ord, cluster = cluster, diagnosticPlot = diagnosticPlot))
}

## Summarize clustered pathways producing a table and graph of the top pathway for each cluster
## input is a series of vectors with essential category enrichment details
## clusters can come from the findCategoryOverlap function or from your favorite algorithm.

summarizeCategoryClusters <- function(catID, catName, catScore, catFDR, catSize, cluster){
  summaryTable <- data.frame(catID, catName, catScore, catFDR, catSize, cluster, stringsAsFactors = FALSE) %>%
    group_by(cluster) %>%
    mutate(clusterSize = n()) %>%
    dplyr::arrange(catFDR, desc(abs(catScore))) %>%
    mutate(myRank = 1:n()) %>%
    ungroup() %>%
    dplyr::filter(myRank == 1) %>%
    dplyr::select(-myRank)

  myLabeller <- function(x, size = 40){
    output <- c()
    for(i in 1:length(x)){
      output[i] <- paste(strwrap(x[i], size), collapse = "\n")
    }
    output
  }

  summaryPlot <- summaryTable %>%
    mutate(catFDR = -log10(catFDR)) %>%
    ggplot(aes(x = reorder(catName, catScore), y = catScore, fill = catSize)) +
    geom_col() +
    scale_fill_gradient(low = grey(0.8), high = grey(0.2), name = "Number of \n proteins") +
    labs(x = "", y = "Score") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme(axis.text.y = element_text(size = 8)) +
    labs(y = "-log10 p-value") +
    coord_flip()

  return(list(table = summaryTable, plot = summaryPlot))
}
