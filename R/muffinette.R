#' @title muffinette
#'
#' @description Meta-analysis of networks using a pseudo-value approach
#'
#' @param metaAbd a feature-by-sample matrix of abundances across all studies.
#' @param batchvar a character vector mapping each observation to its originating study.
#' @param exposurevar a character vector mapping each observation to its level of a binary exposure variable.
#' @param metaData a data frame of metadata, columns must include exposure, batch as factor, covariates to adjust for, if any.
#' @param filter logical. If TRUE, features are filtered according to specified abundance and prevalence thresholds. Default value: TRUE.
#' @param abd_threshold numeric; threshold for abundance of features, only used when \code{filter = TRUE}. Default value: 0.
#' @param prev_threshold numeric; threshold for prevalence of features, only used when \code{filter = TRUE}. Default value: 0.1.
#' @param topfeatures number of features with highest variability to retain for analysis. Default value is \code{NULL} for which all features are kept.
#' @param batchCorrect logical. If TRUE, batch effect correction is performed. Default value: TRUE.
#' @param count logical. If TRUE, the abundances are count. Default value: FALSE.
#' @param net.est.method method for network estimation. Default value: "SparCC".
#' @param covariates optional covariates for batch effect correction. Default value: NULL.
#' @param ncores number of cores to use for parallelized leave-one-subject-out network estimation.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param ... additional arguments passed to \code{\link{networkEst}} according to the choice specified for \code{net.est.method}.
#' @return a list with
#' \item{metafits}{data frame with columns as feature names, effect sizes, p-values and q-values, heterogeneity statistics.}
#' \item{pseudoValues}{sample-by-feature matrix of pseudovalues.}
#'
#' @import parallel
#' @import dplyr
#' @import gtools
#' @import SpiecEasi
#' @import caret
#' @import MMUPHin
#' @import stats
#' @export

muffinette <- function(metaAbd, batchvar, exposurevar, metaData,
                       filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
                       batchCorrect = TRUE, count, net.est.method,
                       covariates = NULL, ncores = 4, verbose = TRUE, ...){


    metaAbd <- check_features_abd(metaAbd)

    ni <- as.numeric(table(batchvar)) ## vector of sample sizes for different studies
    ni_ends <- cumsum(ni)
    ni_starts <- c(1, ni_ends[-length(ni)] + 1)
    if(sum(ni) != nrow(metaData)) {
        stop("Total number of samples in feature_abd_list mismatch with metaData!")
    }

    nstudy <- length(ni) ## number of studies
    if(nstudy < 2) {
        stop("At least two studies are required for meta-analysis!")
    }
    studies <- unique(batchvar) ## vector of names of studies


    ##################################################
    ########### filtering abundant features ##########
    ##################################################
    if(filter) {
        if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
            stop("Threshold values needed for filtering!")
        }
        filtered_featuretable <- filter_abdfeatures(x = t(metaAbd), abd_threshold = abd_threshold,
                                                    prev_threshold = prev_threshold)
    } else {
        filtered_featuretable <- t(metaAbd)
    }
    filtered_featuretable <- filter_varfeatures(x = filtered_featuretable, topV = topfeatures)
    ##################################################


    feature_abd_list_filtered <- lapply(seq_len(nstudy), function(i) {
        filtered_featuretable[ni_starts[i]:ni_ends[i], ]
    })
    feature_names_list <- lapply(feature_abd_list_filtered, colnames)


    #########################################
    ############ Batch correction ###########
    #########################################
    if(batchCorrect) {
        meta_abd_mat <- t(as.matrix(do.call(rbind, feature_abd_list_filtered))) ## feature_by_sample matrix
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
    } else {
        meta_abd_mat <- t(as.matrix(do.call(rbind, feature_abd_list_filtered)))
        data_meta <- data.frame(sampleID = colnames(meta_abd_mat),
                                study = as.factor(rep(studies, ni)),
                                exposure = exposurevar)
        rownames(data_meta) <- data_meta$sampleID
        batch_corrected_abd_df <- as.data.frame(t(meta_abd_mat))
    }
    #########################################


    #################################################################
    #### Network estimation & pseudo-value calculation per study ####
    #################################################################
    feature_abd_list_batchcor <- vector("list", nstudy)
    pseudoVal_list <- vector("list", nstudy)
    estimatedNet <- vector("list", nstudy)

    uniqueExposure <- unique(exposurevar)

    for(i in 1:nstudy) {
        feature_abd_list_batchcor[[i]] <- batch_corrected_abd_df[rownames(feature_abd_list_filtered[[i]]), ]

        groupA <- which(data_meta[data_meta$study == studies[i], 3] == uniqueExposure[1])
        groupB <- which(data_meta[data_meta$study == studies[i], 3] == uniqueExposure[2])
        estimatedNet_groupA <- networkEst(x = feature_abd_list_batchcor[[i]][groupA, ], count = count, estimethod = net.est.method, ...)
        estimatedNet_groupB <- networkEst(x = feature_abd_list_batchcor[[i]][groupB, ], count = count, estimethod = net.est.method, ...)
        estimatedNet[[i]] <- list(estimatedNet_groupA = estimatedNet_groupA, estimatedNet_groupB = estimatedNet_groupB)
        thetahat_groupA <- measureNetwork(estimatedNet_groupA)
        message("DONE")
        thetahat_groupB <- measureNetwork(estimatedNet_groupB)
        estimatedNet_drop_groupA <- parallel::mclapply(1:nrow(feature_abd_list_batchcor[[i]][groupA, ]), function(j) networkEst(x = feature_abd_list_batchcor[[i]][groupA, ][-j, ], count = count, estimethod = net.est.method, ...), mc.cores = ncores)
        estimatedNet_drop_groupB <- parallel::mclapply(1:nrow(feature_abd_list_batchcor[[i]][groupB, ]), function(j) networkEst(x = feature_abd_list_batchcor[[i]][groupB, ][-j, ], count = count, estimethod = net.est.method, ...), mc.cores = ncores)
        thetahat_drop_groupA <- sapply(estimatedNet_drop_groupA, measureNetwork)
        message("DONE")
        thetahat_drop_groupB <- sapply(estimatedNet_drop_groupB, measureNetwork)
        pseudoVal_groupA <- pseudoValue(thetahat_groupA, thetahat_drop_groupA, nrow(feature_abd_list_batchcor[[i]][groupA, ]))
        pseudoVal_groupB <- pseudoValue(thetahat_groupB, thetahat_drop_groupB, nrow(feature_abd_list_batchcor[[i]][groupB, ]))
        pseudoVal <- matrix(NA, nrow(feature_abd_list_batchcor[[i]]), ncol(feature_abd_list_batchcor[[i]]))
        pseudoVal[groupA, ] <- pseudoVal_groupA
        pseudoVal[groupB, ] <- pseudoVal_groupB
        pseudoVal <- scale(pseudoVal, center = TRUE, scale = TRUE)
        colnames(pseudoVal) <- colnames(feature_abd_list_batchcor[[i]]) # Map the column names (feature names)

        pseudoVal_list[[i]] <- pseudoVal
        if(verbose) {
            message(sprintf("Network estimated for study %d / %d", i, nstudy))
        }
    }

    pseudoVal_abd <- do.call(rbind, pseudoVal_list) ## sample-by-feature matrix of pseudovalues
    rownames(pseudoVal_abd) <- rownames(data_meta)
    #########################################


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

    list(metafits = metafits, pseudoValues = pseudoVal_abd)

}
