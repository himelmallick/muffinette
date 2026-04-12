#' @title extractScore
#'
#' @description function to compute sample-by-feature pseudovalue matrix from sample-by-feature abundance matrix.
#' @param x a sample-by-feature abundance matrix.
#' @param cl cluster environment for resampled network estimation in parallel.
#' @param comp logical; TRUE for compositional data.
#' @param net.est.method network estimation method specified by the user according to the (non)compositionality of abundances.
#' @param net_ctrl named list of control parameters specific to the network estimation method specified.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param return_net logical. If TRUE, the estimated network based on all the subjects is returned as part of the output. Default value: TRUE.
#' @return a list with
#' \item{thetahat}{marginal network centrality measures based on all the subjects.}
#' \item{thetahat_loo}{resampled marginal network centrality measures leaving-one-subject at a time.}
#' \item{pseudoVal}{sample-by-feature matrix of pseudovalues.}
#' \item{estimatedNet}{the estimated network based on all the subjects; returned only when \code{returm_net} is TRUE.}
#' @export

extractScore <- function(x, cl, comp, net.est.method, net_ctrl,
                         verbose = TRUE, return_net = TRUE) {
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("'x' must be a matrix or data frame.")
    }

    x <- as.matrix(x)

    if (nrow(x) < 2) {
        stop("Need at least 2 observations in 'x' for leave-one-out estimation.")
    }

    # Full-data network
    estimatedNet <- do.call(networkEst, c(list(x = x,
                                               comp = comp,
                                               estimethod = net.est.method),
                                          net_ctrl))

    ### global network properties for whole study
    thetahat <- measureNetwork(estimatedNet)

    # Export current x to workers
    parallel::clusterExport(cl, varlist = "x", envir = environment())

    # Leave-one-out network estimation in parallel
    estimatedNet_loo <- parallel::parLapply(cl, X = seq_len(nrow(x)),
                                            fun = function(j) {tryCatch(do.call(networkEst,
                                                                                c(list(x = x[-j, , drop = FALSE],
                                                                                       comp = comp,
                                                                                       estimethod = net.est.method),
                                                                                  net_ctrl)),
                                                                        error = function(e) e)})


    failed <- vapply(estimatedNet_loo, inherits, logical(1), what = "error")
    if (any(failed)) {
        stop(sprintf("Leave-One-Out network estimation failed for %d out of %d observations. First error: %s",
                sum(failed), length(failed), estimatedNet_loo[[which(failed)[1]]]$message))
    }

    ### global network properties leaving-one-subject-out within the study
    thetahat_loo <- sapply(estimatedNet_loo, measureNetwork)

    ### jackknife pseudovalues from global network properties
    pseudoVal <- pseudoValue(thetahat, thetahat_loo, nrow(x))

    if (verbose) {
        message(sprintf("extractScore completed for n = %d", nrow(x)))
    }

    out <- list(thetahat = thetahat, thetahat_loo = thetahat_loo,
                pseudoVal = pseudoVal)

    if (return_net) {
        out$estimatedNet <- estimatedNet
        # out$estimatedNet_loo <- estimatedNet_loo
    }

    out
}
