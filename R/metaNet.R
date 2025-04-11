#' @title metaNet
#'
#' @description Meta-analysis of networks using a pseudo-value approach
#'
#' @param feature_abd_list a list of OTU tables of different studies. Each OTU table is a sample-by-feature data frame.
#' @param batchvar name of the batch variable. This variable in data should be a factor variable
#' @param exposurevar name of the exposure variable for differential abundance testing.
#' @param metaData a data frame of metadata, columns must include exposure, batch as factor, covariates to adjust for.
#' @param abd_threshold threshold for abundance of features, cannot be NULL if filter_ft is TRUE.
#' @param prev_threshold threshold for prevalence of features, cannot be NULL if filter_ft is TRUE.
#' @param covariates optional covariates for batch effect correction - default is NULL.
#' @param method method for network estimation - default is sparcc.
#' @param seed.value seed for reproducibility
#' @return list of per-feature meta-analytical differential abundance results such as data frame with columns as effect sizes, p-values and q-values, heterogeneity statistics.
#' @import robustbase
#' @import parallel
#' @import dplyr
#' @import gtools
#' @import fdrtool
#' @import SOHPIE
#' @import caret
#' @import MMUPHin
#' @export


metaNet <- function(feature_abd_list, batchvar, exposurevar,
                    metaData,
                    abd_threshold = NULL, prev_threshold = NULL,
                    covariates = NULL, method, seed.value) {
    set.seed(seed.value)
    
    nstudy <- length(feature_abd_list) ## number of studies
    studies <- unique(batchvar) ## vector of names of studies
    ni <- unlist(lapply(feature_abd_list, nrow)) ## vector of sample sizes for all studies
    if(sum(ni) != nrow(metaData)) {
        stop("Total number of samples mismatch!")
    }
    
    for(i in 1:nstudy) {
        check_features_abd(feature_abd_list[[i]])
    }
    
    #########################################
    ########### filtering features ##########
    #########################################
    
    if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
        stop("Threshold values need to be provided")
    }
    for(i in 1:nstudy) {
        feature_abd_list[[i]] <- filter_features(feature_abd_list[[i]],
                                                 abd_threshold = abd_threshold,
                                                 prev_threshold = prev_threshold)
    }
    
    #########################################
    ############ Common features ############
    #########################################
    
    feature_names_list <- vector("list", nstudy)
    for(i in 1:nstudy) {
        feature_names_list[[i]] <- colnames(feature_abd_list[[i]])
    }
    common_features <- Reduce(intersect, feature_names_list)
    feature_abd_list <- lapply(feature_abd_list,
                               function(foo) foo[, common_features])
    
    #########################################
    ############ Batch correction ###########
    #########################################
    
    meta_abd_mat <- t(as.matrix(do.call(rbind, feature_abd_list))) ## feature_by_sample matrix
    
    data_meta <- data.frame(sampleID = colnames(meta_abd_mat),
                            study = as.factor(rep(studies, ni)),
                            exposure = exposurevar)
    rownames(data_meta) <- data_meta$sampleID
    batch_corrected_abd <- MMUPHin::adjust_batch(feature_abd = meta_abd_mat,
                                                 batch = "study",
                                                 data = data_meta)$feature_abd_adj
    batch_corrected_abd_df <- as.data.frame(t(batch_corrected_abd)); rm(batch_corrected_abd)
    feature_abd_list_batchcor <- vector("list", nstudy)
    for(i in 1:nstudy) {
        feature_abd_list_batchcor[[i]] <- batch_corrected_abd_df[rownames(feature_abd_list[[i]]), ]
    }
    
    message("batchcorrection done")
    
    jackknifes_std_list <- vector("list", nstudy)
    for(i in 1:nstudy) {
        jackknifes_std_list[[i]] <- sohpie_jk(OTUdat = feature_abd_list_batchcor[[i]],
                                              groupA = which(data_meta[data_meta$study == studies[i], 3] == 1),
                                              groupB = which(data_meta[data_meta$study == studies[i], 3] == 0),
                                              method = method, seed = seed.value)
    }
    
    jackknife_std_abd <- do.call("rbind", jackknifes_std_list)
    rownames(jackknife_std_abd) <- rownames(data_meta)
    fit_lmmeta <- MMUPHin::lm_meta(feature_abd = t(jackknife_std_abd),
                          exposure = "exposure",
                          batch = "study",
                          data = data_meta, control = list(forest_plot = "forest.pdf",
                                                           normalization = 'NONE',
                                                           transform = 'NONE'))
    metafits <- fit_lmmeta$meta_fits[order(abs(fit_lmmeta$meta_fits$qval.fdr), decreasing = FALSE), ]
    list(feature_abd_list_batchcor = feature_abd_list_batchcor,
         jackknife_std_abd = jackknife_std_abd,
         metafits = metafits)
    #feature_abd_list_batchcor
}
