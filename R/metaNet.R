#' @title metaNet
#'
#' @description Meta-analysis of networks using a pseudo-value approach
#'
#' @param feature_abd_list a list of OTU tables of different studies. Each OTU table is a sample-by-feature data frame.
#' @param batchvar name of the batch variable. This variable in data should be a factor variable.
#' @param exposurevar name of the exposure variable for differential abundance testing.
#' @param metaData a data frame of metadata, columns must include exposure, batch as factor, covariates to adjust for.
#' @param count logical. If TRUE, the abundances are count. Default value: FALSE.
#' @param net.est.method method for network estimation. Default value: "SparCC".
#' @param abd_threshold threshold for abundance of features.
#' @param prev_threshold threshold for prevalence of features.
#' @param covariates optional covariates for batch effect correction. Default value: NULL.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param ... additional arguments passed to \code{\link{networkEst}}.
#' @return list of per-feature meta-analytical differential abundance results such as data frame with columns as effect sizes, p-values and q-values, heterogeneity statistics.
#' @import parallel
#' @import dplyr
#' @import gtools
#' @import SpiecEasi
#' @import caret
#' @import MMUPHin
#' @export


metaNet <- function(feature_abd_list, batchvar, exposurevar,
                    metaData, count, net.est.method,
                    abd_threshold = NULL, prev_threshold = NULL,
                    covariates = NULL, verbose = TRUE, ...){

    #set.seed(seed.value)

    nstudy <- length(feature_abd_list) ## number of studies
    if(nstudy < 2) {
        stop("There should be at least two studies!")
    }
    studies <- unique(batchvar) ## vector of names of studies
    # Check if feature abundances are neither missing nor negative.
    feature_abd_list <- sapply(1:nstudy, function(i) check_features_abd(feature_abd_list[[i]])$feature_abd_cleaned)
    ni <- unlist(lapply(feature_abd_list, nrow)) ## vector of sample sizes for all studies

    # metaData <- metadata[!(metadata$samplid %in% result$removed_rows), , drop=FALSE]
    # Check if number of samples mismatch
    if(sum(ni) != nrow(metaData)) {
        stop("Total number of samples in feature_abd_list mismatch with metaData!")
    }

    #########################################
    ########### filtering features ##########
    #########################################
    if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
        stop("Threshold values need to be provided")
    }
    feature_names_list <- vector("list", nstudy)
    for(i in 1:nstudy) {
        feature_abd_list[[i]] <- filter_features(feature_abd_list[[i]],
                                                 abd_threshold = abd_threshold,
                                                 prev_threshold = prev_threshold)
        feature_names_list[[i]] <- colnames(feature_abd_list[[i]])
    }

    #########################################
    ####### Extracting common features ######
    #########################################
    common_features <- Reduce(intersect, feature_names_list)
    feature_abd_list <- lapply(feature_abd_list,
                               function(foo) foo[, common_features])
    if(verbose)
        message("Common features across studies extracted...")

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
    if(verbose)
        message("Batch correction done...")

    #################################################################
    #### Network estimation & pseudo-value calculation per study ####
    #################################################################
    feature_abd_list_batchcor <- vector("list", nstudy)
    pseudoVal_list <- vector("list", nstudy)
    estimatedNet <- vector("list", nstudy)

    uniqueExposure <- unique(exposurevar)

    for(i in 1:nstudy) {
        feature_abd_list_batchcor[[i]] <- batch_corrected_abd_df[rownames(feature_abd_list[[i]]), ]

        groupA <- which(data_meta[data_meta$study == studies[i], 3] == uniqueExposure[1])
        groupB <- which(data_meta[data_meta$study == studies[i], 3] == uniqueExposure[2])
        estimatedNet_groupA <- networkEst(X = feature_abd_list_batchcor[[i]][groupA, ], count = count, est.method = net.est.method, ...)
        estimatedNet_groupB <- networkEst(X = feature_abd_list_batchcor[[i]][groupB, ], count = count, est.method = net.est.method, ...)
        estimatedNet[[i]] <- list(estimatedNet_groupA = estimatedNet_groupA, estimatedNet_groupB = estimatedNet_groupB)
        thetahat_groupA <- measureNetwork(estimatedNet_groupA)
        thetahat_groupB <- measureNetwork(estimatedNet_groupB)
        estimatedNet_drop_groupA <- parallel::mclapply(1:nrow(feature_abd_list_batchcor[[i]][groupA, ]), function(j) networkEst(X = feature_abd_list_batchcor[[i]][groupA, ][-j, ], count = count, est.method = net.est.method, ...))
        estimatedNet_drop_groupB <- parallel::mclapply(1:nrow(feature_abd_list_batchcor[[i]][groupB, ]), function(j) networkEst(X = feature_abd_list_batchcor[[i]][groupB, ][-j, ], count = count, est.method = net.est.method, ...))
        thetahat_drop_groupA <- sapply(estimatedNet_drop_groupA, measureNetwork)
        thetahat_drop_groupB <- sapply(estimatedNet_drop_groupB, measureNetwork)
        pseudoVal_groupA <- pseudoValueJK(thetahat_groupA, thetahat_drop_groupA, nrow(feature_abd_list_batchcor[[i]][groupA, ]))
        pseudoVal_groupB <- pseudoValueJK(thetahat_groupB, thetahat_drop_groupB, nrow(feature_abd_list_batchcor[[i]][groupB, ]))
        pseudoVal <- matrix(NA, nrow(feature_abd_list_batchcor[[i]]), ncol(feature_abd_list_batchcor[[i]]))
        pseudoVal[groupA, ] <- pseudoVal_groupA
        pseudoVal[groupB, ] <- pseudoVal_groupB
        pseudoVal <- scale(pseudoVal, center = TRUE, scale = TRUE)
        colnames(pseudoVal) <- colnames(feature_abd_list_batchcor[[i]]) # Map the column names (feature names)

        pseudoVal_list[[i]] <- pseudoVal
        message("Study", i)
    }

    pseudoVal_abd <- do.call(rbind, pseudoVal_list) ## sample-by-feature matrix of pseudovalues
    rownames(pseudoVal_abd) <- rownames(data_meta)

    #########################################
    #### Meta-analysis of pseudo-values  ####
    #########################################
    fit_lmmeta <- MMUPHin::lm_meta(feature_abd = t(pseudoVal_abd),
                                   exposure = "exposure",
                                   batch = "study",
                                   data = data_meta, control = list(forest_plot = "forest.pdf",
                                                                    normalization = 'NONE',
                                                                    transform = 'NONE'))
    if(verbose){
        message("Meta-analysis done...")
    }
    metafits <- fit_lmmeta$meta_fits[order(abs(fit_lmmeta$meta_fits$qval.fdr), decreasing = FALSE), ]

    list(metafits = metafits, pseudoVal_abd = pseudoVal_abd,
         batch_corrected_abd_df = batch_corrected_abd_df,
         feature_abd_list_batchcor = feature_abd_list_batchcor,
         estimatedNet = estimatedNet)

}

