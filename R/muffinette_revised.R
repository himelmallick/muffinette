#' @title muffinette_revised
#'
#' @description Meta-analysis of feature abundances using a network connectivity approach.
#' Normalized abundances must be provided for per-study network construction and connectivity score calculation.
#' Standardized connectivity scores are corrected for batch (study) effect for subsequent meta-analysis.
#'
#' @param metaAbd a feature-by-sample matrix of abundances across all studies. The abundances must be normalized.
#' @param batchvar a character vector mapping each observation to its originating study.
#' @param exposurevar a character vector mapping each observation to its level of a binary exposure variable, else a vector of values of a continuous exposure variable.
#' @param metaData a data frame of metadata, columns must include exposure, batch as factor, covariates to adjust for, if any.
#' @param prevfilter logical; if TRUE, features are filtered according to specified abundance and prevalence thresholds. Default value: TRUE.
#' @param filter_perstudy logical; if TRUE, prevalence and/or variance filtering is applied on each study instead of \code{metaAbd}. Default value: TRUE.
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
#'   See Details in [MMUPHin::maaslin_meta()]. Default values are taken from internal settings.
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

muffinette_revised <- function(metaAbd, batchvar, exposurevar, metaData,
                       prevfilter = TRUE, filter_perstudy = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
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

    sample_names <- rownames(metaData)

    if(is.null(sample_names)) {
        stop("`metaData` must have sample IDs as row names.")
    }

    if(anyDuplicated(sample_names)) {
        stop("Sample IDs in `rownames(metaData)` must be unique.")
    }

    if(is.null(colnames(metaAbd))) {
        stop("`metaAbd` must have sample IDs as column names.")
    }

    if(anyDuplicated(colnames(metaAbd))) {
        stop("Sample IDs in `colnames(metaAbd)` must be unique.")
    }

    ## Check if the colnames in metaAbd and rownames in metaData
    ## are same or not, identify missing/extra samples in metaAbd considering
    ## metaData as the reference, and if the dimension names are same, the
    ## columns of metaAbd are ordered according to the rows (samples) in metaData
    if(!setequal(colnames(metaAbd), sample_names)) {
        missing_in_abd <- setdiff(sample_names, colnames(metaAbd))
        extra_in_abd <- setdiff(colnames(metaAbd), sample_names)

        stop(paste0(
                "`metaAbd` and `metaData` do not contain the same sample IDs.",
                "\nMissing from `metaAbd`: ",
                paste(missing_in_abd, collapse = ", "),
                "\nExtra in `metaAbd`: ",
                paste(extra_in_abd, collapse = ", ")))
    }

    ## Reorder metaAbd columns to exactly match metaData rows
    metaAbd <- metaAbd[, sample_names, drop = FALSE]
    stopifnot(identical(colnames(metaAbd), rownames(metaData)))

    ## Align batch and exposure variable with metaData
    if(length(batchvar) != length(sample_names)) {
        stop("`batchvar` must have one value per sample.")
    }

    if(length(exposurevar) != length(sample_names)) {
        stop("`exposurevar` must have one value per sample.")
    }

    if(!is.null(names(batchvar))) {
        if(!setequal(names(batchvar), sample_names)) {
            stop("Names of `batchvar` must match `rownames(metaData)`.")
        }
        batchvar <- batchvar[sample_names]
    }

    if(!is.null(names(exposurevar))) {
        if(!setequal(names(exposurevar), sample_names)) {
            stop("Names of `exposurevar` must match `rownames(metaData)`.")
        }
        exposurevar <- exposurevar[sample_names]
    }

    names(batchvar) <- sample_names
    names(exposurevar) <- sample_names

    ## Align fixed-effect covariates with metaData
    if(!is.null(covariates)) {

        covariates <- as.data.frame(covariates)
        if(nrow(covariates) != length(sample_names)) {
            stop("`covariates` must have one row per sample.")
        }
        covariate_names <- rownames(covariates)

        ## Determine whether informative row names were supplied
        has_covariate_names <- !is.null(covariate_names) &&
            !identical(covariate_names, as.character(seq_len(nrow(covariates))))

        if(has_covariate_names) {
            if(anyDuplicated(covariate_names)) {
                stop("Row names of `covariates` must be unique.")
            }
            if(!setequal(covariate_names, sample_names)) {
                stop("`rownames(covariates)` must match ",
                    "`rownames(metaData)`.")
            }

            ## Reorder covariate rows to match metaData
            covariates <- covariates[sample_names, , drop = FALSE]

        } else {
            ## Without sample IDs, the existing row order
            ## must correspond to metaData
            rownames(covariates) <- sample_names
        }

        ## Convert character covariates to factors
        char_covs <- vapply(covariates, is.character, logical(1))

        if(any(char_covs)) {
            covariates[char_covs] <- lapply(covariates[char_covs],
                                            as.factor)
        }

        stopifnot(identical(rownames(covariates), sample_names))
    }


    ## Align random-effect covariates with metaData

    if(!is.null(covariates_random)) {
        covariates_random <- as.data.frame(covariates_random)

        if(nrow(covariates_random) != length(sample_names)) {
            stop("`covariates_random` must have one row per sample.")
        }

        random_covariate_names <- rownames(covariates_random)

        ## Determine whether informative row names were supplied
        has_random_covariate_names <- !is.null(random_covariate_names) &&
            !identical(random_covariate_names,
                       as.character(seq_len(nrow(covariates_random))))

        if(has_random_covariate_names) {
            if(anyDuplicated(random_covariate_names)) {
                stop("Row names of `covariates_random` must be unique.")
            }

            if(!setequal(random_covariate_names, sample_names)) {
                stop("`rownames(covariates_random)` must match ",
                    "`rownames(metaData)`.")
            }

            ## Reorder rows to match metaData
            covariates_random <- covariates_random[sample_names, , drop = FALSE]

        } else {
            ## Without sample IDs, assume the current row order
            ## corresponds to metaData
            rownames(covariates_random) <- sample_names
        }

        ## Convert character random-effect variables to factors
        char_covs_random <- vapply(covariates_random, is.character,
                                   logical(1))

        if(any(char_covs_random)) {covariates_random[char_covs_random] <- lapply(
            covariates_random[char_covs_random],
            as.factor)
        }

        stopifnot(identical(rownames(covariates_random), sample_names))
    }

    # batchvar <- as.factor(batchvar)
    batchvar <- factor(batchvar, levels = unique(as.character(batchvar)))
    names(batchvar) <- sample_names
    studies <- levels(batchvar) ## vector of names of studies

    samples_by_study <- split(sample_names,
                              factor(batchvar, levels = studies))

    ni <- lengths(samples_by_study)
    nstudy <- length(samples_by_study)

    ## Construct study-wise abundance matrix
    feature_abd_list <- lapply(studies, function(study) {
        sample_ids <- samples_by_study[[study]]
        metaAbd[, sample_ids, drop = FALSE]
    })

    names(feature_abd_list) <- studies

    ##################################################
    ########### filtering abundant features ##########
    ##################################################
    if(prevfilter == TRUE & filter_perstudy == TRUE) {
        if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
            stop("Threshold values needed for prevalence filtering!")
        }

        for(study in studies) {

            ## Prevalence filtering by converting feature-by-sample to sample-by-feature matrix
            feature_abd_list[[study]] <- filter_abdfeatures(
                x = t(feature_abd_list[[study]]),
                abd_threshold = abd_threshold,
                prev_threshold = prev_threshold)

            ## Still sample-by-feature matrix
            feature_abd_list[[study]] <- filter_varfeatures(
                x = feature_abd_list[[study]],
                topV = topfeatures)
        }
        feature_names_list <- lapply(feature_abd_list, colnames)
        common_features <- Reduce(intersect, feature_names_list)
        if(length(common_features) == 0) {
            stop("No common features remain after study-wise filtering!")
        }
        feature_abd_list <- lapply(studies, function(study) {
            feature_abd_list[[study]][samples_by_study[[study]],
                                      common_features,
                                      drop = FALSE]
        })
        names(feature_abd_list) <- studies
        #feature_abd_list <- lapply(feature_abd_list,
        #                           function(foo) foo[, common_features])

    } else if(prevfilter == TRUE & filter_perstudy == FALSE) {
        if(is.null(abd_threshold) == TRUE | is.null(prev_threshold) == TRUE) {
            stop("Threshold values needed for prevalence filtering!")
        }
        filtered_featuretable <- filter_abdfeatures(x = t(metaAbd), abd_threshold = abd_threshold,
                                                    prev_threshold = prev_threshold) ## sample-by-feature data frame
        filtered_featuretable <- filter_varfeatures(x = filtered_featuretable, topV = topfeatures) ## sample-by-feature data frame

    } else if(prevfilter == FALSE & filter_perstudy == TRUE) {
        message("Prevalence filtering is ignored; only variance filtering is applied per study...")
        for(study in studies) {
            feature_abd_list[[study]] <- filter_varfeatures(
                x = t(feature_abd_list[[study]]),
                topV = topfeatures) ## sample-by-feature data frame
        }
        feature_names_list <- lapply(feature_abd_list, colnames)
        common_features <- Reduce(intersect, feature_names_list)
        if(length(common_features) == 0) {
            stop("No common features remain after study-wise filtering!")
        }
        feature_abd_list <- lapply(studies, function(study) {
            feature_abd_list[[study]][samples_by_study[[study]],
                                      common_features,
                                      drop = FALSE]
        })
        names(feature_abd_list) <- studies
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

    names(pseudoVal_list) <- studies
    names(estimatedNet) <- studies

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

    if (comp && net.est.method == "SpiecEasi") {
        parallel::clusterEvalQ(cl, {
            if(!requireNamespace("SpiecEasi", quietly = TRUE)) {
                stop("Package 'SpiecEasi' is required.")
            }
            NULL
        })
    }

    parallel::clusterExport(cl, varlist = c("networkEst", "comp",
                                            "net.est.method", "net_ctrl"),
                            envir = environment())

    if(is.character(exposurevar) || is.factor(exposurevar)) {
        if(anyNA(exposurevar)) {
            stop("Missing exposure values are not currently supported.")
        }

        exposure_levels <- unique(as.character(exposurevar))
        if(length(exposure_levels) != 2L) {
            stop("If exposure is categorical, it must have exactly two categories.")
        }
    }

    for(i in seq_along(studies)) {
        study <- studies[i]
        study_samples <- samples_by_study[[study]]

        if(filter_perstudy) {
            feature_abd <- feature_abd_list[[study]]
        } else {
            feature_abd <- filtered_featuretable[
                study_samples, , drop = FALSE]
        }

        stopifnot(identical(rownames(feature_abd), study_samples))

        if(is.character(exposurevar) || is.factor(exposurevar)) {
                exposure_study <- as.character(exposurevar[study_samples])
                groupA <- which(exposure_study == exposure_levels[1])
                groupB <- which(exposure_study == exposure_levels[2])
                xA <- feature_abd[groupA, , drop = FALSE]
                xB <- feature_abd[groupB, , drop = FALSE]

                if(nrow(xA) < 10 | nrow(xB) < 10) {
                    if(verbose) {
                        message(sprintf("Skipping study %s: too few samples per exposure (nA = %d, nB = %d)", studies[i], nrow(xA), nrow(xB)))
                    }
                    # estimatedNet[[study]] <- NULL
                    # pseudoVal_list[[study]] <- NULL
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
                estimatedNet[[study]] <- list(estimatedNet_groupA = resA$estimatedNet,
                                          estimatedNet_groupB = resB$estimatedNet)


                pseudoVal <- matrix(NA, nrow(feature_abd), ncol(feature_abd),
                                    dimnames = dimnames(feature_abd))
                pseudoVal[groupA, ] <- pseudoVal_groupA
                pseudoVal[groupB, ] <- pseudoVal_groupB
            } else {
            res <- extractScore(x = feature_abd, cl = cl, comp = comp,
                                net.est.method = net.est.method,
                                net_ctrl = net_ctrl)

            pseudoVal <- res$pseudoVal
            estimatedNet[[study]] <- res$estimatedNet
        }

        pseudoVal <- scale(pseudoVal, center = TRUE, scale = TRUE)
        dimnames(pseudoVal) <- dimnames(feature_abd) # Map the dim names (feature names)

        pseudoVal_list[[study]] <- pseudoVal ## standardized pseudovalues
        rm(pseudoVal)
        if(verbose) {
            cat(sprintf("Network estimated for study %d / %d\n", i, nstudy))
        }
    }

    pseudoVal_list <- Filter(Negate(is.null), pseudoVal_list)

    if(length(pseudoVal_list) == 0) {
        stop("No studies remained after filtering/skipping. `pseudoVal_abd` will be `NULL`.")
    }

    pseudo_feature_names <- lapply(pseudoVal_list, colnames)

    same_pseudo_features <- all(
        vapply(pseudo_feature_names, identical, logical(1),
               pseudo_feature_names[[1]]))

    if(!same_pseudo_features) {
        stop("Pseudovalue matrices do not have identical feature columns ",
            "in identical order across studies.")
    }

    pseudoVal_abd <- do.call(rbind, pseudoVal_list) ## sample-by-feature matrix of pseudovalues
    #########################################

    #########################################
    n_remain <- length(pseudoVal_list)

    if(n_remain > 1) {
        keep_samples <- rownames(pseudoVal_abd)

        data_meta <- data.frame(sampleID = keep_samples,
                                study = droplevels(factor(batchvar[keep_samples])),
                                exposure = exposurevar[keep_samples])

        if(!is.null(covariates)) {
            data_meta <- cbind(data_meta, covariates[keep_samples, , drop = FALSE])
        }

        if(!is.null(covariates_random)) {
            data_meta <- cbind(data_meta, covariates_random[keep_samples, , drop = FALSE])
        }
        rownames(data_meta) <- data_meta$sampleID
        stopifnot(identical(rownames(data_meta), keep_samples),
                  identical(colnames(pseudoVal_abd), colnames(pseudoVal_list[[1]])),
                  identical(colnames(t(pseudoVal_abd)), keep_samples))


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
        covariate_names <- if(is.null(covariates)) {
            NULL
        } else {
            colnames(covariates)
        }

        random_covariate_names <- if(is.null(covariates_random)) {
            NULL
        } else {
            colnames(covariates_random)
        }
        fit_lmmeta <- lm_meta_muff(feature_abd = pseudoVal_abd_batchCtd,
                                   exposure = "exposure",
                                   batch = "study",
                                   covariates = covariate_names,
                                   covariates_random = random_covariate_names,
                                   data = data_meta, control = meta_ctrl_use)
        if(verbose){
            cat("Meta-analysis done...")
        }
        metafits <- fit_lmmeta$meta_fits[order(abs(fit_lmmeta$meta_fits$qval), decreasing = FALSE), ]
        output <- list(metafits = metafits, pseudoValues = pseudoVal_abd,
                       pseudoValues_batchCtd = pseudoVal_abd_batchCtd,
                       estimatedNet = estimatedNet,
                       metaAbd = metaAbd,
                       metaData = metaData)
    } else {
        output <- list(pseudoValues = pseudoVal_abd,
                       estimatedNet = estimatedNet,
                       metaAbd = metaAbd,
                       metaData = metaData)
    }

    output
}


