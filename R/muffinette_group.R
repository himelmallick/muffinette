#' @title muffinette_group
#'
#' @description Meta-analysis of feature abundances using a network connectivity approach.
#' Normalized abundances must be provided for per-study network construction and connectivity score calculation.
#' Standardized connectivity scores are corrected for batch (study) effect for subsequent meta-analysis.
#'
#' @param metaAbd a feature-by-sample matrix of abundances across all studies. The abundances must be normalized.
#' @param feature_abd_list a list of feature-by-sample matrices of abundances per study; this is necessary only when \code{filter_perstudy = TRUE}.
#' @param batchvar a character vector mapping each observation to its originating study.
#' @param exposurevar a character vector mapping each observation to its level of a binary exposure variable, else a vector of values of a continuous exposure variable.
#' @param metaData a data frame of metadata, columns must include exposure, batch as factor, covariates to adjust for, if any.
#' @param prevfilter logical; if TRUE, features are filtered according to specified abundance and prevalence thresholds. Default value: TRUE.
#' @param filter_perstudy logical; if TRUE, prevalence and/or variance filtering is applied on each study instead of \code{metaAbd}.
#' @param abd_threshold numeric; threshold for abundance of features, used only when \code{prevfilter = TRUE}. Default value: 0.
#' @param prev_threshold numeric; threshold for prevalence of features, used only when \code{prevfilter = TRUE}. Default value: 0.1.
#' @param topfeatures number of features with highest variability to retain for analysis. Default value is \code{NULL} such that all the remaining features after removing the near zero predictors are kept.
#' @param comp logical. If FALSE, the abundances are non-compositional and network construction step is avoided. Default value: TRUE.
#' @param net.est.method method for network estimation. Default value: "SparCC".
#' @param covariates optional covariates to adjust for in Maaslin2 model. Default value: NULL.
#' @param covariates_random optional random effects to adjust for in Maaslin2 model. Default value: NULL.
#' @param ncores number of cores to use for leave-one-subject-out network estimation in parallel.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param fixseed seed for reproducibility
#' @param control A named list of control parameters. May contain the following
#' components:
#' \describe{
#'   \item{meta}{A named list of control arguments for meta-analysis.}
#'   \item{network}{A named list of additional arguments passed to the
#'    network estimation method specified by \code{net.est.method}.}
#' }
#'
#' Each component is optional and contains specific arguments, see Details.
#' Missing components are replaced by defaults.
#' @details
#' The \code{control} argument allows fine-grained control of different stages
#' of the pipeline:
#'
#' \itemize{
#'   \item \code{control$meta}: Controls the meta-analysis step for compositional
#'   data, including \code{analysis_method} parameter for Maaslin2 and meta-analysis method.
#'   See Details in [MMUPHin::lm_meta()]. Default values are taken from internal settings.
#'
#'   \item \code{control$network}: Additional arguments passed to the selected
#'   network estimation method. The valid arguments depend on
#'   \code{net.est.method}.
#' }
#'
#' Any unspecified arguments are set to their default values.
#' @return a list with
#' \item{metafits}{data frame with columns as feature names, effect sizes, p-values and q-values, heterogeneity statistics.}
#' \item{pseudoValues}{sample-by-feature matrix of pseudovalues.}
#' @export

muffinette_group <- function(metaAbd, feature_abd_list = NULL, batchvar, exposurevar, metaData,
                       prevfilter = TRUE, filter_perstudy = FALSE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
                       comp = TRUE, net.est.method = "SparCC",
                       covariates = NULL, covariates_random = NULL, ncores = 4, verbose = TRUE, fixseed = NULL, control = NULL){

    if(is.null(control))
        control <- list()
    if(is.null(control$meta))
        control$meta <- list()
    if(is.null(control$network))
        control$network <- list()
    meta_ctrl <- utils::modifyList(control_lm_meta, control$meta)
    net_ctrl <- control$network

    metaAbd <- check_features_abd(metaAbd, comp = comp)

    ni <- as.numeric(table(factor(batchvar, levels = unique(batchvar)))) ## vector of sample sizes for different studies
    ni_ends <- cumsum(ni)
    ni_starts <- c(1, ni_ends[-length(ni)] + 1)
    if(sum(ni) != nrow(metaData)) {
        stop("Total number of samples in feature_abd_list mismatch with `metaData`!")
    }
    nstudy <- length(ni) ## number of studies

    # if(comp) {
    #     if(!all(round(colSums(metaAbd)) == 1)) {
    #         stop("If comp is TRUE, abundances must be compositional.")
    #     }
    # }

    sample_names <- rownames(metaData)

    if(is.null(sample_names)) {
        sample_names <- paste0("sample", seq_len(nrow(metaData)))
        rownames(metaData) <- sample_names
    }

    names(batchvar) <- sample_names
    names(exposurevar) <- sample_names

    if(!is.null(covariates)) {
        rownames(covariates) <- sample_names
    }

    if(!is.null(covariates_random)) {
        rownames(covariates_random) <- sample_names
    }

    batchvar <- as.factor(batchvar)
    names(batchvar) <- sample_names
    studies <- unique(batchvar) ## vector of names of studies


    ##################################################
    ########### filtering abundant features ##########
    ##################################################
    if(prevfilter == TRUE & filter_perstudy == TRUE) {
        if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
            stop("Threshold values needed for prevalence filtering!")
        }
        if(is.null(feature_abd_list)) {
            stop("Abundances are needed as a list for study-wise filtering!")
        }
        for(i in 1:nstudy) {
            feature_abd_list[[i]] <- filter_abdfeatures(t(feature_abd_list[[i]]),
                                                        abd_threshold = abd_threshold,
                                                        prev_threshold = prev_threshold) ## sample-by-feature data frame
            feature_abd_list[[i]] <- filter_varfeatures(feature_abd_list[[i]],
                                                        topV = topfeatures) ## sample-by-feature data frame
        }
        feature_names_list <- lapply(feature_abd_list, colnames)
        common_features <- Reduce(intersect, feature_names_list)
        feature_abd_list <- lapply(feature_abd_list,
                                   function(foo) foo[, common_features])
    } else if(prevfilter == TRUE & filter_perstudy == FALSE) {
        if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
            stop("Threshold values needed for prevalence filtering!")
        }
        filtered_featuretable <- filter_abdfeatures(x = t(metaAbd), abd_threshold = abd_threshold,
                                                    prev_threshold = prev_threshold) ## sample-by-feature data frame
        filtered_featuretable <- filter_varfeatures(x = filtered_featuretable, topV = topfeatures) ## sample-by-feature data frame

    } else if(prevfilter == FALSE & filter_perstudy == TRUE) {
        if(is.null(feature_abd_list)) {
            stop("Abundances are needed as a list for study-wise filtering!")
        }
        message("Prevalence filtering is ignored; only variance filtering is applied per study...")
        for(i in 1:nstudy) {
            feature_abd_list[[i]] <- filter_varfeatures(t(feature_abd_list[[i]]),
                                                        topV = topfeatures) ## sample-by-feature data frame
        }
        feature_names_list <- lapply(feature_abd_list, colnames)
        common_features <- Reduce(intersect, feature_names_list)
        feature_abd_list <- lapply(feature_abd_list,
                                   function(foo) foo[, common_features])
    } else {
        message("Prevalence filtering is ignored; only variance filtering is applied on `metaAbd`...")
        filtered_featuretable <- filter_varfeatures(x = t(metaAbd), topV = topfeatures) ## sample-by-feature data frame
    }
    ##################################################



    #################################################################
    #### Network estimation & pseudovalue calculation per study #####
    #################################################################
    pseudoVal_list <- vector("list", nstudy)
    estimatedNet <- vector("list", nstudy)

    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if(!is.null(fixseed)) {
        # to ensure master uses reproducible parallel RNG
        RNGkind("L'Ecuyer-CMRG")
        set.seed(fixseed)

        # to ensure workers use same RNG engine
        parallel::clusterEvalQ(cl, {
            RNGkind("L'Ecuyer-CMRG")
            NULL
        })

        # generate independent streams for workers
        parallel::clusterSetRNGStream(cl, iseed = fixseed)
    }
    parallel::clusterEvalQ(cl, {
        if(!requireNamespace("SpiecEasi", quietly = TRUE)) {
            stop("Package 'SpiecEasi' is required.")
        }
        NULL
    })
    parallel::clusterExport(cl, varlist = c("networkEst", "comp",
                                            "net.est.method", "net_ctrl"),
                            envir = environment())
    for(i in 1:nstudy) {
        if(filter_perstudy) {
            feature_abd <- feature_abd_list[[i]]
        } else {
            feature_abd <- filtered_featuretable[ni_starts[i]:ni_ends[i], ]
        }

        if(is.character(exposurevar)) {
            if(length(unique(exposurevar)) != 2) {
                stop("If exposure is categorical, it must have two categories!")
            } else {
                uniqueExposure <- unique(exposurevar)
                groupA <- which(exposurevar[batchvar == studies[i]] == uniqueExposure[1])
                groupB <- which(exposurevar[batchvar == studies[i]] == uniqueExposure[2])
                xA <- feature_abd[groupA, , drop = FALSE]
                xB <- feature_abd[groupB, , drop = FALSE]

                if(nrow(xA) < 5 | nrow(xB) < 5) {
                    if(verbose) {
                        message(sprintf("Skipping study %s: too few samples per exposure (nA = %d, nB = %d)", studies[i], nrow(xA), nrow(xB)))
                    }
                    estimatedNet[[i]] <- NULL
                    pseudoVal_list[[i]] <- NULL
                    next
                }

                resA <- extractScore(x = xA, cl = cl, comp = comp,
                                     net.est.method = net.est.method,
                                     net_ctrl = net_ctrl)

                pseudoVal_groupA <- resA$pseudoVal
                resB <- extractScore(x = xB, cl = cl, comp = comp,
                                     net.est.method = net.est.method,
                                     net_ctrl = net_ctrl)

                pseudoVal_groupB <- resB$pseudoVal
                estimatedNet[[i]] <- list(estimatedNet_groupA = resA$estimatedNet,
                                          estimatedNet_groupB = resB$estimatedNet)


                pseudoVal <- matrix(NA, nrow(feature_abd), ncol(feature_abd))
                pseudoVal[groupA, ] <- pseudoVal_groupA
                pseudoVal[groupB, ] <- pseudoVal_groupB
            }
        } else {
            res <- extractScore(x = feature_abd, cl = cl, comp = comp,
                                net.est.method = net.est.method,
                                net_ctrl = net_ctrl)

            pseudoVal <- res$pseudoVal
            estimatedNet[[i]] <- res$estimatedNet
        }

        pseudoVal <- scale(pseudoVal, center = TRUE, scale = TRUE)
        colnames(pseudoVal) <- colnames(feature_abd) # Map the column names (feature names)

        pseudoVal_list[[i]] <- pseudoVal ## standardized pseudovalues
        rm(pseudoVal)
        if(verbose) {
            cat(sprintf("Network estimated for study %d / %d\n", i, nstudy))
        }
    }

    pseudoVal_list <- Filter(Negate(is.null), pseudoVal_list)

    if(length(pseudoVal_list) == 0) {
        stop("No studies remained after filtering/skipping. `pseudoVal_abd` will be `NULL`.")
    }

    pseudoVal_abd <- do.call(rbind, pseudoVal_list) ## sample-by-feature matrix of pseudovalues
    #dimnames(pseudoVal_abd) <- dimnames(filtered_featuretable)
    colnames(pseudoVal_abd) <- colnames(pseudoVal_list[[1]])
    rownames(pseudoVal_abd) <- unlist(lapply(pseudoVal_list, rownames), use.names = FALSE)
    #########################################
    #########################################

    if(nstudy > 1) {
        keep_samples <- rownames(pseudoVal_abd)

        data_meta <- data.frame(sampleID = keep_samples,
                                study = as.factor(batchvar[keep_samples]),
                                exposure = exposurevar[keep_samples])

        if(!is.null(covariates)) {
            data_meta <- cbind(data_meta, covariates[keep_samples, , drop = FALSE])
        }

        if(!is.null(covariates_random)) {
            data_meta <- cbind(data_meta, covariates_random[keep_samples, , drop = FALSE])
        }
        rownames(data_meta) <- data_meta$sampleID
        # data_meta <- data.frame(sampleID = rownames(pseudoVal_abd),
        #                             study = as.factor(batchvar),
        #                             exposure = exposurevar,
        #                             covariates = covariates,
        #                             covariates_random = covariates_random)
        # rownames(data_meta) <- data_meta$sampleID

        #########################################
        ### Batch-correction of pseudo-values ###
        #########################################
        mod_combat <- stats::model.matrix(~ exposure, data = data_meta)
        pseudoVal_abd_batchCtd <- sva::ComBat(dat = t(pseudoVal_abd),
                                              batch = data_meta$study,
                                              mod = mod_combat) ## feature-by-sample matrix batch-corrected pseudovalues


        #########################################
        #### Meta-analysis of pseudo-values  ####
        #########################################
        meta_ctrl_use <- utils::modifyList(meta_ctrl, list(output = getwd(),
                                                           forest_plot = "forest.pdf",
                                                           normalization = "NONE",
                                                           transform = "NONE"))
        fit_lmmeta <- lm_meta_muff(feature_abd = pseudoVal_abd_batchCtd,
                                   exposure = "exposure",
                                   batch = "study",
                                   covariates = covariates,
                                   covariates_random = covariates_random,
                                   data = data_meta, control = meta_ctrl_use)
        if(verbose){
            cat("Meta-analysis done...")
        }
        metafits <- fit_lmmeta$meta_fits[order(abs(fit_lmmeta$meta_fits$qval), decreasing = FALSE), ]
        output <- list(metafits = metafits, pseudoValues = pseudoVal_abd,
                       pseudoValues_batchCtd = pseudoVal_abd_batchCtd,
                       estimatedNet = estimatedNet)
    } else {
        output <- list(pseudoValues = pseudoVal_abd,
                       estimatedNet = estimatedNet)
    }

    output
}


