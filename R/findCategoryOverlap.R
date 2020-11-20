#' Category Overlap
#'
#' Creates an overlap matrix from a list of category members
#' Overlap is calculated as overlap/(size of the smaller cluster)
#'
#' @param membershipList Named list where each element is a category and contains a vector of protein identifiers
#'
#' @return Matrix of overlaps, each between 0 - 1.  Matrix names come from make.names(names(membershipList))
#'
#' @examples
#' data("membershipList")
#'
#' findCategoryOverlap(membershipList = membershipList)
#'
#' @export
#'
#' @import stats
findCategoryOverlap <- function(membershipList){
  ## get overlap matrix
  OG <- mat.or.vec(nr=length(membershipList), nc=length(membershipList))

  namesMap <- make.names(names(membershipList))
  names(namesMap) <- names(membershipList)

  rownames(OG) <- make.names(names(membershipList)); colnames(OG) <- make.names(names(membershipList))
  for (i in 1:nrow(OG)) for(j in 1:ncol(OG)){
    OG[i,j]<-length(intersect(membershipList[[i]], membershipList[[j]]))/
      min(c(length(membershipList[[i]]), length(membershipList[[j]])))
  }
  OG
}
