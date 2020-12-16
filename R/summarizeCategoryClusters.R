#' Summarize categories
#'
#' produce a table and a graph of the top pathway for each cluster
#' Input is a series of vectors with essential category enrichment details
#'
#' @param catName vector of category names
#'
#' @param catScore vector of test statistics for each category
#'
#' @param catFDR vector of FDR corrected p-values
#'
#' @param catSize vector with number of proteins per category
#'
#' @param cluster vector indicating which cluster a category is assigned to
#'
#' @return A list containing
#'      summaryTable: a table with one entry per cluster,
#'      SummaryPlot: a bargraph with one bar per cluster,
#'
#'
#' @export
#' @import stats
#' @importFrom ggplot2 ggplot aes
#' @importFrom dplyr %>% arrange
summarizeCategoryClusters <- function(catName, catScore, catFDR, catSize, cluster){
  summaryTable <- data.frame(catName, catScore, catFDR, catSize, cluster, stringsAsFactors = FALSE) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(clusterSize = n()) %>%
    dplyr::arrange(catFDR, desc(abs(catScore))) %>%
    dplyr::mutate(myRank = 1:n()) %>%
    dplyr::ungroup() %>%
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
    dplyr::mutate(catFDR = -log10(catFDR)) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(catName, catScore), y = catScore, fill = catSize)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = grey(0.8), high = grey(0.2), name = "Number of \n proteins") +
    ggplot2::labs(x = "", y = "Score") +
    ggplot2::scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ggplot2::theme(axis.text.y = element_text(size = 8)) +
    ggplot2::labs(y = "Mean t-statistic") +
    coord_flip()

  return(list(table = summaryTable, plot = summaryPlot))
}
