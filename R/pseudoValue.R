#' @title pseudoValue
#'
#' @description function to compute the Jackknife pseudovalues of network measures.
#'
#' @param thetahatinput network measure of all subjects.
#' @param thetahatdropinput network measure excluding a particular subject.
#' @param sizegroup number of observations.
#' @return matrix of jackknife pseudovalues.
#' @export


# Calculate jackknife pseudovalues (\tilde{\theta}_{ik} using \hat{\theta}_{k} and \hat{\theta}_{k(i)}.
pseudoValue <- function(thetahatinput, thetahatdropinput, sizegroup) {
    thetatildeout <- matrix(NA, ncol = length(thetahatinput), nrow = sizegroup)
    thetatildeout <- sapply(1:nrow(thetahatdropinput), function(k) {
        sizegroup * thetahatinput[k] - (sizegroup - 1) * thetahatdropinput[k, ]
    })
    thetatildeout
}