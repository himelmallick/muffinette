#' @title muffinette
#'
#' @description Meta-analysis of networks using a pseudo-value approach
#'
#' @param metaAbd a feature-by-sample matrix of abundances across all studies. The abundances must be normalized.
#' @param batchvar a character vector mapping each observation to its originating study.
#' @param exposurevar a character vector mapping each observation to its level of a binary exposure variable, else a vector of values of a continuous exposure variable.
#' @param metaData a data frame of metadata, columns must include exposure, batch as factor, covariates to adjust for, if any.
#' @param filter logical. If TRUE, features are filtered according to specified abundance and prevalence thresholds. Default value: TRUE.
#' @param abd_threshold numeric; threshold for abundance of features, only used when \code{filter = TRUE}. Default value: 0.
#' @param prev_threshold numeric; threshold for prevalence of features, only used when \code{filter = TRUE}. Default value: 0.1.
#' @param topfeatures number of features with highest variability to retain for analysis. Default value is \code{NULL} such that all the remaining features after removing the near zero predictors are kept.
#' @param comp logical. If FALSE, the abundances are non-compositional and network construction step is avoided. Default value: TRUE.
#' @param net.est.method method for network estimation. Default value: "SparCC".
#' @param covariates optional covariates for batch effect correction. Default value: NULL.
#' @param ncores number of cores to use for leave-one-subject-out network estimation in parallel.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param fixseed seed for reproducibility
#' @param control A named list of control parameters. May contain the following
#' components:
#' \describe{
#'   \item{batch}{A named list of control arguments for batch adjustment.}
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
#'   \item \code{control$batch}: Controls batch adjustment with a named list of
#'   components, such as \code{zero_inflation}.
#'   See Details in [MMUPHin::adjust_batch()]. Default values are taken from internal settings.
#'
#'   \item \code{control$meta}: Controls the meta-analysis step for compositional
#'   data, including normalization, transformation and meta-analysis method.
#'   See Details in [MMUPHin::lm_meta()]. Default values are taken from internal settings.
#'
#'   \item \code{control$network}: Additional arguments passed to the selected
#'   network estimation method. The valid arguments depend on
#'   \code{net.est.method}.
#' }
#'
#' Any unspecified arguments fall back to their default values.
#' @return a list with
#' \item{metafits}{data frame with columns as feature names, effect sizes, p-values and q-values, heterogeneity statistics; returned only when there are multiple studies.}
#' \item{pseudoValues}{sample-by-feature matrix of pseudovalues.}
#' \item{pseudoValues_batchCtd}{sample-by-feature matrix of batch-corrected pseudovalues, in presence of multiple studies.}
#' \item{estimatedNet}{the estimated network based on all the subjects; returned only when \code{returm_net} is TRUE.}
#' @export

muffinette <- function(metaAbd, batchvar, exposurevar, metaData,
                       filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
                       comp = TRUE, net.est.method = "SparCC",
                       covariates = NULL, ncores = 4, verbose = TRUE, fixseed = NULL, control = NULL){

    if(is.null(control))
        control <- list()
    if(is.null(control$batch))
        control$batch <- list()
    if(is.null(control$meta))
        control$meta <- list()
    if(is.null(control$network))
        control$network <- list()
    batch_ctrl <- utils::modifyList(control_adjust_batch, control$batch)
    meta_ctrl <- utils::modifyList(control_lm_meta, control$meta)
    net_ctrl <- control$network

    metaAbd <- check_features_abd(metaAbd)

    ni <- as.numeric(table(factor(batchvar, levels = unique(batchvar)))) ## vector of sample sizes for different studies
    ni_ends <- cumsum(ni)
    ni_starts <- c(1, ni_ends[-length(ni)] + 1)
    if(sum(ni) != nrow(metaData)) {
        stop("Total number of samples in feature_abd_list mismatch with metaData!")
    }

    if(comp) {
        if(!all(round(colSums(metaAbd)) == 1)) {
            stop("If comp is TRUE, abundances must be compositional.")
        }
    }

    nstudy <- length(ni) ## number of studies
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
        filtered_featuretable <- t(metaAbd) ## sample-by-feature data frame
    }
    ##################################################

    #print(dim(filtered_featuretable))
    #################################################################
    #### Network estimation & pseudovalue calculation per study #####
    #################################################################
    #feature_abd_list_batchcor <- vector("list", nstudy)
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
        feature_abd <- filtered_featuretable[ni_starts[i]:ni_ends[i], ]
        res <- extractScore(x = feature_abd, cl = cl, comp = comp,
                            net.est.method = net.est.method,
                            net_ctrl = net_ctrl)
        pseudoVal <- scale(res$pseudoVal, center = TRUE, scale = TRUE)
        colnames(pseudoVal) <- colnames(filtered_featuretable) # Map the column names (feature names)

        pseudoVal_list[[i]] <- pseudoVal
        estimatedNet[[i]] <- res$estimatedNet
        rm(pseudoVal)
        if(verbose) {
            cat(sprintf("Network estimated for study %d / %d\n", i, nstudy))
        }
    }

    pseudoVal_abd <- do.call(rbind, pseudoVal_list) ## sample-by-feature matrix of pseudovalues
    dimnames(pseudoVal_abd) <- dimnames(filtered_featuretable)
    #########################################
    #########################################

    if(nstudy > 1) {
        if(!is.null(covariates)) {
            data_meta <- data.frame(sampleID = rownames(pseudoVal_abd),
                                    study = as.factor(batchvar),
                                    exposure = exposurevar,
                                    covariates = covariates)
        } else {
            data_meta <- data.frame(sampleID = rownames(pseudoVal_abd),
                                    study = as.factor(batchvar),
                                    exposure = exposurevar)
        }
        rownames(data_meta) <- data_meta$sampleID

        #########################################
        ### Batch-correction of pseudo-values ###
        #########################################
        mod_combat <- stats::model.matrix(~ exposure, data = data_meta)
        pseudoVal_abd_batchCtd <- sva::ComBat(dat = t(pseudoVal_abd),
                                              batch = batchvar,
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


