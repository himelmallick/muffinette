#' @title measureNetwork
#'
#' @description A function to compute the network centrality (i.e. total connectivity)
#'         of each microbial taxa from the association matrix.
#'
#' @param asso.matinput An input is an association matrix that is estimated from
#'        the user-provided OTU data.
#'
#' @return A vector containing network centrality of each taxa.
#'
#' @export

## Total connectivity (similar to PRANA)
# measureNetwork <- function(asso.matinput) {
#     results <- vector()
#     for(j in 1:ncol(asso.matinput)) {
#         results[j] <- sum(asso.matinput[j, -j])
#     }
#     return(results)
# }

measureNetwork <- function(asso.matinput) {
    rowSums(asso.matinput) - diag(asso.matinput)
}


