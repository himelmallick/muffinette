#' @title measureNetwork
#'
#' @description function to compute a network property for each microbial taxa from the estimated association matrix.
#' @param asso.matinput an association matrix estimated using a specified network estimation method.
#' @return A vector of network measures of each taxa.
#' @export

### network centrality (total connectivity)
measureNetwork <- function(asso.matinput) {
    colSums(as.matrix(asso.matinput)) - diag(as.matrix(asso.matinput))
}


