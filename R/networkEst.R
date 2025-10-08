#' @title networkEst
#'
#' @description Meta-analysis of networks using a pseudo-value approach
#'
#' @param X abundance matrix
#' @param count logical. If TRUE, abundances are count instead of proportions. Default value: FALSE.
#' @param est.method network estimation method. Options include "SparCC" and "SpiecEasi".
#' @param ... additional arguments appropriate for the network estimation method chosen.
#' @return an association matrix.
#' @import gtools
#' @import SpiecEasi
#' @export

##################################################
#### Function to implement network estimation ####
##################################################
networkEst <- function(X, count = FALSE, est.method, ...){

  ### Implement SparCC method ###
  # if(method == "SparCC"){
  #   if(count == TRUE){
  #     estNetObj <- sparcc(x = X, ...)$cor.w
  #   } else {
  #     estNetObj <- SparCC.frac(x = X, ...)$cor.w
  #   }
  # }
  if(est.method == "SpiecEasi") {
      fit <- SpiecEasi::spiec.easi(data = as.matrix(X), ...)
      beta.mat <- as.matrix(getOptBeta(fit))
      adj.mat <- ((beta.mat != 0) + t(beta.mat != 0)) > 0
      estNetObj <- as.numeric(adj.mat)
      dim(estNetObj) <- dim(adj.mat)
  }
  # if(est.method == "SparCC"){
  #     estNetObj <- sparcc(x = X, ...)$cor.w
  # }
  if(est.method == "SparCC"){
      estNetObj <- SpiecEasi::sparcc(data = X, ...)$Cor
  }
  ### Implement CCLasso method ###
  # if(est.method == "CCLasso"){
  #   estNetObj <- cclasso(x = X, counts = count, ...)$cor_w
  # }

  ### Implement GLasso method ###
  # if(est.method == "glasso"){
  #   X <- as.matrix(X)
  #   lambda_seq <- exp(seq(log(1), log(0.01), length.out = 50))
  #   fit <- huge(x = X, method = method, lambda = lambda_seq)
  #   fit_select <- huge.select(est = fit, ...)
  #   estNetObj <- fit_select$opt.icov
  # }

  return(estNetObj)
}
