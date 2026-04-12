#' @title networkEst
#'
#' @description function to estimate network based on feature abundances.
#'
#' @param x abundance matrix
#' @param comp logical. If FALSE, the abundances are non-compositional. Default value: TRUE.
#' @param estimethod network estimation method. Options include "SparCC" and "SpiecEasi". Default value: SparCC.
#' @param ... additional arguments appropriate for the network estimation method chosen.
#' @return an association matrix.
#' @export

networkEst <- function(x, comp = TRUE, estimethod = "SparCC", ...){

    estimethod <- match.arg(estimethod,
                            choices = c("SparCC", "SpiecEasi", "glasso"))

    if(estimethod == "SparCC") {
        estNetObj <- SpiecEasi::sparcc(data = x, ...)$Cor
    } else if(estimethod == "SpiecEasi") {
        fit <- SpiecEasi::spiec.easi(data = as.matrix(x), ...)
        beta.mat <- as.matrix(SpiecEasi::getOptBeta(fit))
        adj.mat <- ((beta.mat != 0) + t(beta.mat != 0)) > 0
        estNetObj <- as.numeric(adj.mat)
        dim(estNetObj) <- dim(adj.mat)
    }

    extra_args <- list(...)

    ### Implement GLasso method ###
    if(estimethod == "glasso"){
      x <- as.matrix(x)

      huge_args <- extra_args$huge
      select_args <- extra_args$huge.select

      if (is.null(huge_args)) huge_args <- list()
      if (is.null(select_args)) select_args <- list()

      if (is.null(huge_args$lambda)) {
          huge_args$lambda <- exp(seq(log(1), log(0.01), length.out = 25))
      }

      fit <- do.call(huge::huge, c(list(x = x), huge_args))
      fit_select <- do.call(huge::huge.select, c(list(est = fit), select_args))
      estNetObj <- fit_select$opt.icov
    }

    estNetObj
}

