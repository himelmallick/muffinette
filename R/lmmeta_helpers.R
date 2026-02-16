

check_exposure <- function (exposure, batch)
{
    ind_exposure <- as.vector(tapply(exposure, batch, function(x) length(setdiff(unique(x),
                                                                                 NA)) > 1))
    names(ind_exposure) <- levels(batch)
    if (is.factor(exposure)) {
        lvl_exposure <- levels(exposure)
        ind_exposure_cat <- as.vector(tapply(exposure, batch,
                                             function(x) all(lvl_exposure %in% x)))
        if (any(ind_exposure & !ind_exposure_cat))
            stop("Exposure is character/factor and does not have common levels ",
                 "in the following batches.\n", paste(names(ind_exposure)[ind_exposure &
                                                                    !ind_exposure_cat], collapse = ", "))
    }
    return(ind_exposure)
}

rename_maaslin <- function (old_names, prefix)
{
    if (is.null(old_names) | length(old_names) == 0)
        return(NULL)
    new_names <- paste0(prefix, seq_along(old_names))
    names(new_names) <- old_names
    return(new_names)
}

catchToList <- function (expr)
{
    val <- NULL
    myWarnings <- NULL
    wHandler <- function(w) {
        myWarnings <<- c(myWarnings, w$message)
        invokeRestart("muffleWarning")
    }
    myError <- NULL
    eHandler <- function(e) {
        myError <<- e$message
        NULL
    }
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler),
                    error = eHandler)
    list(value = val, warnings = myWarnings, error = myError)
}

shorten_name <- function (x, cutoff = 3, replacement = "..")
{
    if (anyDuplicated(x))
        stop("x should be unique character strings!")
    x_sub <- x
    length_x <- nchar(x_sub, type = "c")
    ind_change <- length_x > cutoff * 2 + nchar(replacement,
                                                type = "c")
    x_sub[ind_change] <- paste0(substr(x_sub[ind_change], start = 1,
                                       stop = cutoff), replacement, substr(x_sub[ind_change],
                                                                           start = length_x[ind_change] - cutoff + 1, stop = length_x[ind_change]))
    if (anyDuplicated(x_sub))
        return(x)
    else return(x_sub)
}

check_covariates <- function (data_covariates, batch)
{
    ind_covariates <- matrix(NA, nrow = nlevels(batch), ncol = ncol(data_covariates))
    dimnames(ind_covariates) <- list(levels(batch), names(data_covariates))
    ind_covariates[] <- vapply(data_covariates, function(covariate) as.vector(tapply(covariate,
                                                                                     batch, function(x) {
                                                                                         length(unique(x[!is.na(x)])) > 1
                                                                                     })), rep_len(TRUE, nlevels(batch)))
    return(ind_covariates)
}

check_covariates_random <- function (data_covariates, batch)
{
    ind_covariates_random <- matrix(NA, nrow = nlevels(batch),
                                    ncol = ncol(data_covariates))
    dimnames(ind_covariates_random) <- list(levels(batch), names(data_covariates))
    ind_covariates_random[] <- vapply(data_covariates, function(covariate) as.vector(tapply(covariate,
                                                                                            batch, function(x) {
                                                                                                length(unique(x[!is.na(x)])) > 1 & any(table(x) >
                                                                                                                                           1)
                                                                                            })), rep_len(TRUE, nlevels(batch)))
    if (all(!ind_covariates_random) & ncol(ind_covariates_random) >
        0)
        stop("Random covariates are provided,", " but no batch has clustered observations!")
    return(ind_covariates_random)
}

create_table_maaslin <- function (features, exposure, lvl_exposure)
{
    if (is.null(lvl_exposure))
        values_exposure <- exposure
    else values_exposure <- lvl_exposure[-1]
    names(features) <- NULL
    table_maaslin <- expand.grid(features, exposure, values_exposure,
                                 stringsAsFactors = FALSE)
    names(table_maaslin) <- c("feature", "metadata", "value")
    return(table_maaslin)
}

Maaslin2_wrapper <- function (feature_abd, data, exposure, covariates = NULL, covariates_random = NULL,
          output = tempdir(), normalization = "TSS", transform = "AST",
          analysis_method = "LM")
{
    feature_abd_rename <- feature_abd
    data_rename <- data[, c(exposure, covariates, covariates_random),
                        drop = FALSE]
    features_rename <- rename_maaslin(rownames(feature_abd_rename),
                                      prefix = "T")
    samples_rename <- rename_maaslin(colnames(feature_abd_rename),
                                     prefix = "S")
    exposure_rename <- rename_maaslin(exposure, prefix = "E")
    covariates_rename <- rename_maaslin(covariates, prefix = "X")
    covariates_random_rename <- rename_maaslin(covariates_random,
                                               prefix = "RX")
    dimnames(feature_abd_rename) <- list(features_rename, samples_rename)
    dimnames(data_rename) <- list(samples_rename, c(exposure_rename,
                                                    covariates_rename, covariates_random_rename))
    ind_features <- apply(feature_abd_rename > 0, 1, any)
    message_maaslin <- utils::capture.output(log_maaslin <- catchToList(Maaslin2::Maaslin2(input_data = feature_abd_rename[ind_features,
                                                                                                                    , drop = TRUE], input_metadata = data_rename, output = output,
                                                                                    min_abundance = 0, min_prevalence = 0, normalization = normalization,
                                                                                    transform = transform, analysis_method = analysis_method,
                                                                                    max_significance = 1, random_effects = covariates_random_rename,
                                                                                    fixed_effects = c(exposure_rename, covariates_rename),
                                                                                    standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE)$results))
    if (!is.null(log_maaslin$error)) {
        ch_error <- log_maaslin$error
        variables_rename <- c(exposure_rename, covariates_rename,
                              covariates_random_rename)
        for (i_variable in names(variables_rename)) {
            i_pattern <- paste0("'", variables_rename[i_variable],
                                "'")
            i_pattern_replace <- paste0("'", i_variable, "'")
            ch_error <- gsub(i_pattern, i_pattern_replace, x = ch_error,
                             fixed = TRUE)
        }
        stop("Internal Maaslin run error!\n", ch_error)
    }
    res_rename <- log_maaslin$value
    lvl_exposure <- NULL
    if (is.factor(data[[exposure]]))
        lvl_exposure <- levels(data[[exposure]])
    table_maaslin <- dplyr::left_join(data.frame(feature = names(features_rename),
                                                 feature_rename = features_rename, stringsAsFactors = FALSE),
                                      create_table_maaslin(features_rename, exposure_rename,
                                                           lvl_exposure), by = c(feature_rename = "feature"))
    res <- dplyr::left_join(table_maaslin, res_rename, by = c(feature_rename = "feature",
                                                              "metadata", "value"))
    res <- dplyr::select(res, -feature_rename, -name)
    res$metadata <- exposure
    if (all(res$value == exposure_rename))
        res$value <- exposure
    res$qval <- stats::p.adjust(res$pval, method = "fdr")
    return(res)
}

rma_wrapper <- function (maaslin_fits, method = "REML", output = tempdir(),
          forest_plot = NULL, rma_conv = 1e-06, rma_maxit = 1000,
          pvalAdjust = "fdr", alpha_thresh = 0.05, verbose = TRUE) {
    lvl_batch <- names(maaslin_fits)
    n_batch <- length(lvl_batch)
    exposure <- unique(maaslin_fits[[1]]$metadata)
    #cat("expo: ", exposure, "\n")
    values_exposure <- unique(maaslin_fits[[1]]$value)
    #cat("value_expo: ", values_exposure, "\n")
    if (isTRUE(verbose)) {
        message("expo: ", paste(exposure, collapse = ", "))
        message("value_expo: ", paste(values_exposure, collapse = ", "))
    }

    features <- unique(maaslin_fits[[1]]$feature)
    l_results <- list()
    for (value_exposure in values_exposure) {
        i_result <- data.frame(matrix(NA, nrow = length(features),
                                      ncol = 11 + length(lvl_batch)))
        rma_fits <- vector("list", length(features))
        names(rma_fits) <- features
        colnames(i_result) <- c("feature", "exposure", "coef",
                                "stderr", "pval", "k", "tau2", "stderr.tau2", "pval.tau2",
                                "I2", "H2", paste0("weight_", lvl_batch))
        i_result$feature <- features
        i_result$exposure <- value_exposure
        rownames(i_result) <- i_result$feature
        if (!is.null(forest_plot))
            grDevices::pdf(paste0(output, "/", exposure, "_", value_exposure,
                       "_", forest_plot), width = 6, height = 4 + ifelse(n_batch >
                                                                             4, (n_batch - 4) * 0.5, 0))
        if (any(features != maaslin_fits[[2]][maaslin_fits[[2]]$value ==
                                              value_exposure, "feature"]))
            stop("Feature names don't match between maaslin_fits components!")
        betas <- vapply(maaslin_fits, function(i_maaslin_fit) i_maaslin_fit[i_maaslin_fit$value ==
                                                                                value_exposure, "coef", drop = TRUE], rep_len(0,
                                                                                                                              length(features)))
        sds <- vapply(maaslin_fits, function(i_maaslin_fit) i_maaslin_fit[i_maaslin_fit$value ==
                                                                              value_exposure, "stderr", drop = TRUE], rep_len(0,
                                                                                                                              length(features)))
        pvals <- vapply(maaslin_fits, function(i_maaslin_fit) i_maaslin_fit[i_maaslin_fit$value ==
                                                                                value_exposure, "pval", drop = TRUE], rep_len(0,
                                                                                                                              length(features)))
        rownames(betas) <- rownames(sds) <- rownames(pvals) <- features
        ind_features <- !is.na(betas) & !is.na(sds) & (sds !=
                                                           0)
        count_feature <- apply(ind_features, 1, sum)
        for (feature in features) {
            if (count_feature[feature] >= 2) {
                i_log <- catchToList(metafor::rma.uni(yi = betas[feature,
                                                                 ind_features[feature, ]], sei = sds[feature,
                                                                                                     ind_features[feature, ]], slab = lvl_batch[ind_features[feature,
                                                                                                     ]], method = method, control = list(threshold = rma_conv,
                                                                                                                                         maxiter = rma_maxit)))
                if (!is.null(i_log$error)) {
                    warning("Fitting rma on feature ", feature,
                            ";\n", i_log$error)
                    next
                }
                if (!is.null(i_log$warnings))
                    warning("Fitting rma on feature ", feature,
                            ";\n", i_log$warnings)
                i_rma_fit <- i_log$value
                wts <- metafor::weights.rma.uni(i_rma_fit)
                i_result[feature, c("coef", "stderr", "pval",
                                    "k", "tau2", "stderr.tau2", "pval.tau2", "I2",
                                    "H2", paste0("weight_", names(wts)))] <- c(unlist(i_rma_fit[c("beta",
                                                                                                  "se", "pval", "k", "tau2", "se.tau2", "QEp",
                                                                                                  "I2", "H2")]), wts)
                # if (i_rma_fit$pval < 0.05 & !is.null(forest_plot))
                #     metafor::forest(i_rma_fit, xlab = shorten_name(feature,
                #                                                    cutoff = 15), slab = shorten_name(lvl_batch[ind_features[feature,
                #                                                    ]], cutoff = 5))
                rma_fits[[feature]] <- i_rma_fit
            }
            if (count_feature[feature] == 1) {
                i_ind_features <- ind_features[feature, ]
                tmp_batch <- lvl_batch[i_ind_features]
                i_result[feature, c("coef", "stderr", "pval",
                                    "k", paste0("weight_", tmp_batch))] <- c(betas[feature,
                                                                                   i_ind_features], sds[feature, i_ind_features],
                                                                             pvals[feature, i_ind_features], 1, 100)
            }
        }
        # if (!is.null(forest_plot))
        #     grDevices::dev.off()
        #i_result$pval.bonf <- stats::p.adjust(i_result$pval, method = "bonf")
        #i_result$qval.fdr <- stats::p.adjust(i_result$pval, method = "fdr")
        i_result$qval <- NA_real_
        meta_idx <- which(i_result$k >= 2)
        i_result$qval[meta_idx] <- stats::p.adjust(i_result$pval[meta_idx], method = pvalAdjust)

        if (!is.null(forest_plot)) {
            sig_features <- rownames(i_result)[!is.na(i_result$qval) & i_result$qval < alpha_thresh]
            for (feature in sig_features) {
                fit <- rma_fits[[feature]]
                if (!is.null(fit)) {
                    metafor::forest(
                        fit,
                        xlab = shorten_name(feature, cutoff = 15),
                        slab = shorten_name(lvl_batch[ind_features[feature, ]], cutoff = 5)
                    )
                }
            }
            grDevices::dev.off()
        }
        l_results[[value_exposure]] <- i_result
    }

    results <- Reduce("rbind", l_results)
    return(results)
}


lm_meta_muff <- function (feature_abd, batch, exposure, covariates = NULL, covariates_random = NULL,
          data, control)
{
    control <- match_control(default = control_lm_meta, control = control)
    verbose <- control$verbose
    feature_abd <- as.matrix(feature_abd)
    data <- as.data.frame(data)
    samples <- check_samples(feature_abd = feature_abd, data = data)
    if (length(batch) > 1)
        stop("Only one batch variable is supported!")
    df_batch <- check_metadata(data = data, variables = batch)
    df_meta <- check_metadata(data = data, variables = c(exposure,
                                                         covariates, covariates_random), no_missing = FALSE)
    var_batch <- check_batch(df_batch[[batch]], min_n_batch = 2)
    n_batch <- nlevels(var_batch)
    lvl_batch <- levels(var_batch)
    if (verbose)
        message("Found ", n_batch, " batches")
    if (is.character(df_meta[[exposure]]))
        df_meta[[exposure]] <- factor(df_meta[[exposure]], levels = stringr::str_sort(unique(df_meta[[exposure]])))
    ind_exposure <- check_exposure(df_meta[[exposure]], var_batch)
    if (any(!ind_exposure))
        warning("Exposure variable is missing or has only one non-missing value",
                " in the following batches; Maaslin2 won't be fitted on them\n",
                paste(lvl_batch[!ind_exposure], collapse = ", "))
    ind_covariates <- check_covariates(df_meta[covariates], var_batch)
    for (covariate in covariates) {
        if (any(ind_exposure & !ind_covariates[, covariate]))
            warning("Covariate ", covariate, " is missing or has only one non-missing value",
                    " in the following batches; will be excluded from model for",
                    " these batches:\n", paste(lvl_batch[ind_exposure &
                                                             !ind_covariates[, covariate]], collapse = ", "))
    }
    ind_covariates_random <- check_covariates_random(df_meta[covariates_random],
                                                     var_batch)
    for (covariate in covariates_random) {
        if (!any(ind_exposure & ind_covariates_random[, covariate]))
            warning("Random covariate ", covariate, " has no clustered observations!")
        else if (verbose)
            message("Random covariate ", covariate, "will be fitted for the following batches:\n",
                    paste(lvl_batch[ind_exposure & ind_covariates_random[,
                                                                         covariate]], collapse = ", "))
    }
    dir.create(control$output, recursive = TRUE, showWarnings = FALSE)
    maaslin_fits <- list()
    for (i in seq_len(n_batch)) {
        i_batch <- lvl_batch[i]
        if (!ind_exposure[i_batch])
            next
        if (verbose)
            message("Fitting Maaslin2 on batch ", i_batch, "...")
        i_feature_abd <- feature_abd[, var_batch == i_batch,
                                     drop = FALSE]
        i_data <- df_meta[var_batch == i_batch, , drop = FALSE]
        i_covariates <- covariates[ind_covariates[i_batch, ,
                                                  drop = TRUE]]
        i_covariates_random <- covariates_random[ind_covariates_random[i_batch,
                                                                       , drop = TRUE]]
        i_output <- paste0(control$output, "/", i_batch)
        dir.create(i_output, showWarnings = FALSE)
        i_maaslin <- Maaslin2_wrapper(feature_abd = i_feature_abd,
                                      data = i_data, exposure = exposure, covariates = i_covariates,
                                      covariates_random = i_covariates_random, output = i_output,
                                      normalization = control$normalization, transform = control$transform,
                                      analysis_method = control$analysis_method)
        maaslin_fits[[i_batch]] <- i_maaslin
        maaslin_fits[[i_batch]]$batch <- i_batch
    }
    if (verbose)
        message("Fitting meta-analysis model.")
    meta_fits <- rma_wrapper(maaslin_fits, method = control$rma_method,
                             output = control$output, forest_plot = control$forest_plot,
                             rma_conv = control$rma_conv, rma_maxit = control$rma_maxit,
                             verbose = verbose)
    return(list(meta_fits = meta_fits, maaslin_fits = maaslin_fits,
                control = control))
}


