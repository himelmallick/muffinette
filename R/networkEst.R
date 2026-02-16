#' @title networkEst
#'
#' @description function to estimate network based on feature abundances.
#'
#' @param x abundance matrix
#' @param count logical. If TRUE, abundances are count instead of proportions. Default value: FALSE.
#' @param estimethod network estimation method. Options include "SparCC" and "SpiecEasi". Default value: SparCC.
#' @param ... additional arguments appropriate for the network estimation method chosen.
#' @return an association matrix.
#' @export

networkEst <- function(x, count = FALSE, estimethod = "SparCC", ...){

    estimethod <- match.arg(estimethod,
                            choices = c("SparCC", "SpiecEasi"))

    if(estimethod == "SparCC") {
        estNetObj <- SpiecEasi::sparcc(data = x, ...)$Cor
    } else if(estimethod == "SpiecEasi") {
        fit <- SpiecEasi::spiec.easi(data = as.matrix(x), ...)
        beta.mat <- as.matrix(SpiecEasi::getOptBeta(fit))
        adj.mat <- ((beta.mat != 0) + t(beta.mat != 0)) > 0
        estNetObj <- as.numeric(adj.mat)
        dim(estNetObj) <- dim(adj.mat)
    }

    ### Implement GLasso method ###
    # if(estimethod == "glasso"){
    #   x <- as.matrix(x)
    #   lambda_seq <- exp(seq(log(1), log(0.01), length.out = 50))
    #   fit <- huge(x = x, method = method, lambda = lambda_seq)
    #   fit_select <- huge.select(est = fit, ...)
    #   estNetObj <- fit_select$opt.icov
    # }

    estNetObj
}

