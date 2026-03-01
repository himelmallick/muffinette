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
#' @param analysismethod method of meta-analysis.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param fixseed seed for reproducibility
#' @param ... additional arguments passed to \code{\link{networkEst}} according to the choice specified for \code{net.est.method}.
#' @return a list with
#' \item{metafits}{data frame with columns as feature names, effect sizes, p-values and q-values, heterogeneity statistics.}
#' \item{pseudoValues}{sample-by-feature matrix of pseudovalues.}
#' @export

muffinette <- function(metaAbd, batchvar, exposurevar, metaData,
                       filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
                       batchCorrect = TRUE, count, net.est.method,
                       covariates = NULL, ncores = 4, analysismethod = 'LM', verbose = TRUE, fixseed = NULL, ...){

    if(!is.null(fixseed)) {
        RNGkind("L'Ecuyer-CMRG")
        set.seed(fixseed)
    }


    metaAbd <- check_features_abd(metaAbd)

    ni <- as.numeric(table(factor(batchvar, levels = unique(batchvar)))) ## vector of sample sizes for different studies
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
        filtered_featuretable <- filter_varfeatures(x = filtered_featuretable, topV = topfeatures) ## sample-by-feature data frame
    } else {
        filtered_featuretable <- t(metaAbd)
    }
    ##################################################


    #########################################
    ############ Batch correction ###########
    #########################################
    if(batchCorrect) {
        meta_abd_mat <- t(as.matrix(filtered_featuretable)) ## feature-by-sample matrix
        data_meta <- data.frame(sampleID = colnames(meta_abd_mat),
                                study = as.factor(batchvar),
                                exposure = exposurevar)
        rownames(data_meta) <- data_meta$sampleID
        batch_corrected_abd <- adjust_batch_muff(feature_abd = meta_abd_mat,
                                                     batch = "study",
                                                     data = data_meta)$feature_abd_adj
        batch_corrected_abd_df <- as.data.frame(t(batch_corrected_abd)) ## sample-by-feature data frame
        rm(batch_corrected_abd)
        if(verbose)
            message("Batch correction done...")
    } else {
        data_meta <- data.frame(sampleID = rownames(filtered_featuretable),
                                study = as.factor(batchvar),
                                exposure = exposurevar)
        rownames(data_meta) <- data_meta$sampleID
        batch_corrected_abd_df <- filtered_featuretable
    }
    #########################################


    #################################################################
    #### Network estimation & pseudovalue calculation per study ####
    #################################################################
    feature_abd_list_batchcor <- vector("list", nstudy)
    pseudoVal_list <- vector("list", nstudy)
    estimatedNet <- vector("list", nstudy)

    uniqueExposure <- unique(exposurevar)

    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if(!is.null(fixseed)) {
        # ensure master uses reproducible parallel RNG
        RNGkind("L'Ecuyer-CMRG")
        set.seed(fixseed)

        # ensure workers use same RNG engine
        parallel::clusterEvalQ(cl, {
            RNGkind("L'Ecuyer-CMRG")
            NULL
        })

        # generate independent streams for workers
        parallel::clusterSetRNGStream(cl, iseed = fixseed)
    }
    # if(!is.null(fixseed)) {
    #     parallel::clusterSetRNGStream(cl, iseed = fixseed)
    # }
    parallel::clusterEvalQ(cl, {
        if(!requireNamespace("SpiecEasi", quietly = TRUE)) {
            stop("Package 'SpiecEasi' is required.")
        }
        NULL
    })
    parallel::clusterExport(cl, varlist = c("networkEst", "count", "net.est.method"),
                            envir = environment())

    kept_studies <- character(0)
    for(i in 1:nstudy) {
        feature_abd_list_batchcor[[i]] <- batch_corrected_abd_df[ni_starts[i]:ni_ends[i], ]

        groupA <- which(data_meta$exposure[data_meta$study == studies[i]] == uniqueExposure[1])
        groupB <- which(data_meta$exposure[data_meta$study == studies[i]] == uniqueExposure[2])
        xA <- feature_abd_list_batchcor[[i]][groupA, , drop = FALSE]
        xB <- feature_abd_list_batchcor[[i]][groupB, , drop = FALSE]

        if(nrow(xA) < 5 | nrow(xB) < 5) {
            if(verbose) {
                message(sprintf("Skipping study %s: too few samples per exposure (nA = %d, nB = %d)", studies[i], nrow(xA), nrow(xB)))
            }
            estimatedNet[[i]] <- NULL
            pseudoVal_list[[i]] <- NULL
            next
        }

        estimatedNet_groupA <- networkEst(x = xA, count = count, estimethod = net.est.method, ...)
        estimatedNet_groupB <- networkEst(x = xB, count = count, estimethod = net.est.method, ...)
        ## remove later after checking reproducibility
        estimatedNet[[i]] <- list(estimatedNet_groupA = estimatedNet_groupA, estimatedNet_groupB = estimatedNet_groupB)

        thetahat_groupA <- measureNetwork(estimatedNet_groupA)
        thetahat_groupB <- measureNetwork(estimatedNet_groupB)

        parallel::clusterExport(cl, varlist = c("xA"),
                                envir = environment())

        estimatedNet_loo_groupA <- parallel::parLapply(cl, X = seq_len(nrow(xA)),
                                                   fun = function(j) {
                                                       tryCatch(networkEst(x = xA[-j, , drop = FALSE], count = count, estimethod = net.est.method, ...),
                                                                error = function(e) e)
                                                   })
        parallel::clusterExport(cl, varlist = c("xB"),
                                envir = environment())
        estimatedNet_loo_groupB <- parallel::parLapply(cl, X = seq_len(nrow(xB)),
                                                   fun = function(j) {
                                                       tryCatch(networkEst(x = xB[-j, , drop = FALSE], count = count, estimethod = net.est.method, ...),
                                                                error = function(e) e)
                                                   })

        thetahat_loo_groupA <- sapply(estimatedNet_loo_groupA, measureNetwork)
        thetahat_loo_groupB <- sapply(estimatedNet_loo_groupB, measureNetwork)
        pseudoVal_groupA <- pseudoValue(thetahat_groupA, thetahat_loo_groupA, nrow(xA))
        pseudoVal_groupB <- pseudoValue(thetahat_groupB, thetahat_loo_groupB, nrow(xB))
        pseudoVal <- matrix(NA, nrow(feature_abd_list_batchcor[[i]]), ncol(feature_abd_list_batchcor[[i]]))
        pseudoVal[groupA, ] <- pseudoVal_groupA
        pseudoVal[groupB, ] <- pseudoVal_groupB
        pseudoVal <- scale(pseudoVal, center = TRUE, scale = TRUE)
        colnames(pseudoVal) <- colnames(feature_abd_list_batchcor[[i]]) # Map the column names (feature names)

        pseudoVal_list[[i]] <- pseudoVal
        rm(pseudoVal)
        if(verbose) {
            cat(sprintf("Network estimated for study %d / %d\n", i, nstudy))
        }
        kept_studies <- c(kept_studies, studies[i])
    }

    pseudoVal_list <- Filter(Negate(is.null), pseudoVal_list)
    pseudoVal_abd <- do.call(rbind, pseudoVal_list) ## sample-by-feature matrix of pseudovalues
    rownames(pseudoVal_abd) <- rownames(data_meta)
    # data_meta_eff <- data_meta[rownames(pseudoVal_abd), , drop = FALSE]
    #########################################


    #########################################
    #### Meta-analysis of pseudo-values  ####
    #########################################
    fit_lmmeta <- lm_meta_muff(feature_abd = t(pseudoVal_abd),
                                   exposure = "exposure",
                                   batch = "study",
                                   data = data_meta, control = list(output = getwd(),
                                                                    analysis_method = analysismethod,
                                                                    forest_plot = "forest.pdf",
                                                                    normalization = 'NONE',
                                                                    transform = 'NONE'))
    if(verbose){
        cat("Meta-analysis done...")
    }
    metafits <- fit_lmmeta$meta_fits[order(abs(fit_lmmeta$meta_fits$qval), decreasing = FALSE), ]

    list(metafits = metafits, pseudoValues = pseudoVal_abd,
         feature_abd_list_batchcor = feature_abd_list_batchcor,
         estimatedNet = estimatedNet)

}
