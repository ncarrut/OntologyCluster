#' Cluster Overlapping Categories
#'
#' First, decides which categories are similar enough to be clustered.  The usual similarity structure is that a few categories are almost identical and then the rest are mainly unique.  Just cluster nearly identical ones
#'
#' @param OG matrix of overlap
#'
#' @param similarityCut minimum similarity required for clustering
#'
#' @param clusteringMethod a string, either "hierarchical" or "spectral" to indicate clustering method.
#'
#' @return A list containing
#'      cluster: cluster assignments,
#'      plot: ordered plot to visualize clustering,
#'      plotOrder: ordering for the visualization plot in case you want to try other clusterings,
#'      overlapMatrix: the original matrix
#'      diagnosticPlot: A diagnostic plot for the cluster/no cluster decision
#'
#' @examples
#' data("membershipList")
#'
#' OG <- findCategoryOverlap(membershipList = membershipList)
#'
#' clusterByOverlap(OG, similarityCut = 0.5, clusteringMethod = "hierarchical")
#'
#' @export
#' @import stats
#' @importFrom ggplot2 ggplot aes
#' @importFrom dplyr %>% arrange filter pull
#' @importFrom tidyr pivot_longer
#' @importFrom Spectrum estimate_k
#' @importFrom kernlab specc
clusterByOverlap <- function(OG, similarityCut = 0.5, clusteringMethod = "hierarchical"){
  cutVec <- seq(0, 1, 0.05)
  cutResVec <- c()
  for(i in 1:length(cutVec)){
    cutResVec[i] <- length(unique(cutree(tree = hclust(as.dist(1-OG), method = "ward.D"), h = cutVec[i])))
  }
  diagnosticPlot <- data.frame(cutVec, cutResVec) %>%
    ggplot2::ggplot(ggplot2::aes(x = cutVec, y = cutResVec)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Minimum similarity", y = "Clusters")

  ## decide which categories should be clustered
  hcRes <- cutree(tree = hclust(as.dist(1-OG), method = "ward.D"), h = similarityCut)
  doCluster <- as.character(hcRes) %in% names(table(hcRes))[table(hcRes) > 1]

  if(clusteringMethod == "spectral"){
    ## cluster using spectral clustering
    ## get K from Spectrum package
    myK <- Spectrum::estimate_k(OG[doCluster, doCluster], maxk = sum(doCluster) - 2, showplots = FALSE) %>%
      dplyr::filter(K > sqrt(sum(doCluster))) %>%
      dplyr::filter(z == max(z)) %>%
      dplyr::pull(K)

    ## then assign clusters from specc
    ## except where distance was too great to cluster then
    ## just fill in with a sequence
    cluster1 <- kernlab::specc(OG[doCluster, doCluster], centers = myK)@.Data
    cluster <- rep(0, nrow(OG))
    names(cluster) = rownames(OG)
    cluster[doCluster] <- cluster1
    cluster[!doCluster] <- seq(from = myK+1, to = myK+sum(cluster==0), by = 1)
  }

  if(clusteringMethod == "hierarchical"){
    cluster <- rep(0, nrow(OG))
    names(cluster) = rownames(OG)
    cluster[doCluster] <- hcRes[doCluster]
    cluster[!doCluster] <- seq(from = max(cluster) + 1, to = max(cluster) + sum(!doCluster), by = 1)
  }

  ## order according to euclidean distance.  This is for plotting and
  ## will be returned to help with evaluating other clustering methods if needed
  ord <- hclust(as.dist(1-OG), method = "ward.D")$order

  #plot to visualize clusters and assess clustering
  tilePlot <- data.frame(OG) %>%
    dplyr::mutate(myRowNames = ordered(row.names(OG), levels = rownames(OG)[ord])) %>%
    dplyr::mutate(cluster = cluster) %>%
    tidyr::pivot_longer(-c(myRowNames, cluster), names_to = "myColNames") %>%
    dplyr::mutate(myColNames = ordered(myColNames, levels = rownames(OG)[ord])) %>%
    ggplot2::ggplot(ggplot2::aes(x = myRowNames, y = myColNames, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=6)) +
    ggplot2::scale_fill_gradient(name = "fraction",
                        low = "#FFFFFF",
                        high = "#012345") +
    ggplot2::labs(x = "cluster", y = "pathway") +
    ggplot2::scale_x_discrete(labels = cluster[ord]) +
    ggplot2::scale_y_discrete(label = function(x) strtrim(x, 35))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=6))

  return(list(cluster = cluster, plot = tilePlot, plotOrder = ord,
              overlapMatrix = OG, diagnosticPlot = diagnosticPlot))
}
