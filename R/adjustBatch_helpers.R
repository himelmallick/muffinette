relocate_scale <- function (s_data, l_params_shrink, batchmod, n_batch, l_ind)
{
    adj_data <- s_data
    for (i_batch in seq_len(n_batch)) {
        i_ind_feature <- !is.na(l_params_shrink$gamma_star[,
                                                           i_batch]) & !is.na(l_params_shrink$delta_star[, i_batch])
        if (!all(i_ind_feature == l_ind$ind_gamma[, i_batch]))
            stop("Features determined to be eligible for batch estimation do not ",
                 "agree with the ones with valid per-batch shrinked parameters!")
        for (i_feature in seq_len(nrow(adj_data))) {
            if (i_ind_feature[i_feature]) {
                i_ind_sample <- l_ind$ind_data[i_feature, ] &
                    as.logical(batchmod[, i_batch])
                adj_data[i_feature, i_ind_sample] <- (adj_data[i_feature,
                                                               i_ind_sample] - l_params_shrink$gamma_star[i_feature,
                                                                                                          i_batch])/sqrt(l_params_shrink$delta_star[i_feature,
                                                                                                                                                    i_batch])
            }
        }
    }
    return(adj_data)
}

add_back_covariates <- function (adj_data, l_stand_feature, l_ind)
{
    for (i_feature in seq_len(nrow(adj_data))) {
        if (l_ind$ind_feature[i_feature]) {
            i_stand_feature <- l_stand_feature[[i_feature]]
            adj_data[i_feature, l_ind$ind_data[i_feature, ]] <- adj_data[i_feature,
                                                                         l_ind$ind_data[i_feature, ]] * sqrt(i_stand_feature$var_pooled) +
                i_stand_feature$stand_mean
        }
    }
    return(adj_data)
}

adjust_EB <- function (s_data, l_params_shrink, l_stand_feature, batchmod,
          n_batch, l_ind)
{
    if (n_batch != ncol(batchmod))
        stop("n_batch does not agree with batchmod!")
    if (n_batch != ncol(l_params_shrink[[1]]))
        stop("n_batch does not agree with l_params_shrink!")
    adj_data <- relocate_scale(s_data, l_params_shrink, batchmod,
                               n_batch, l_ind)
    adj_data <- add_back_covariates(adj_data, l_stand_feature,
                                    l_ind)
    return(adj_data)
}

check_feature_abd <- function (feature_abd)
{
    if (any(is.na(feature_abd)))
        stop("Found missing values in the feature table!")
    if (any(feature_abd < 0))
        stop("Found negative values in the feature table!")
    if (all(feature_abd <= 1)) {
        return("proportions")
    }
    else if (all(feature_abd == floor(feature_abd))) {
        return("counts")
    }
    else stop("Feature table does not appear to be either proportions or counts!")
}

check_samples <- function (feature_abd, data)
{
    if (ncol(feature_abd) != nrow(data))
        stop("Dimensions of feature table and metadata table do not agree!")
    if (is.null(rownames(data)))
        stop("data should not have empty row names!")
    if (!identical(colnames(feature_abd), rownames(data)))
        stop("Sample names in feature_abd and data don't agree!")
    return(rownames(data))
}

check_metadata <- function (data, variables, no_missing = TRUE)
{
    if (is.null(variables))
        return(NULL)
    variables_absent <- setdiff(variables, colnames(data))
    if (length(variables_absent) > 0) {
        stop("Following variable(s) not present in data:\n",
             paste(variables_absent, collapse = ","))
    }
    variables_missing <- vapply(variables, function(variable) {
        any(is.na(data[[variable]]))
    }, TRUE)
    if (any(variables_missing) & no_missing) {
        stop("Following variable(s) in data have missing values:\n",
             paste(variables[variables_missing], collapse = ","))
    }
    return(data[, variables, drop = FALSE])
}

check_batch <- function (x, min_n_batch = 2)
{
    if (!is.factor(x)) {
        warning("Batch variable is not a factor as provided and will be converted ",
                "to one.")
        x <- as.factor(x)
    }
    if (nlevels(x) < min_n_batch)
        stop("Must have at least ", min_n_batch, " batches!")
    return(x)
}

check_pseudo_count <- function (x)
{
    if (x < 0)
        stop("pseudo_count must be non-negative")
    return(x)
}

construct_design <- function (data, with_intercept = TRUE)
{
    if (is.null(data))
        return(NULL)
    if (with_intercept)
        stats::model.matrix(~., data = data)
    else stats::model.matrix(~. - 1, data = data)
}

set_pseudo <- function (features)
{
    type_features <- check_feature_abd(features)
    if (all(features == 0))
        stop("All feature abundances are zero!")
    min(setdiff(features, 0))/2
}

check_rank <- function (design)
{
    if (is.null(design))
        return(TRUE)
    qr(design)$rank == ncol(design)
}

construct_ind <- function (feature_abd, n_batch, design, zero_inflation)
{
    ind_data <- matrix(TRUE, nrow(feature_abd), ncol(feature_abd))
    ind_gamma <- matrix(TRUE, nrow(feature_abd), n_batch)
    ind_mod <- rep(TRUE, ncol(design) - n_batch)
    if (zero_inflation) {
        ind_data[feature_abd == 0] <- FALSE
        for (i_feature in seq_len(nrow(feature_abd))) {
            i_design <- design[ind_data[i_feature, ], , drop = FALSE]
            i_check_batch <- apply(i_design[, seq_len(n_batch),
                                            drop = FALSE] == 1, 2, any)
            i_design <- i_design[, c(i_check_batch, ind_mod),
                                 drop = FALSE]
            if (sum(i_check_batch) > 1 && qr(i_design)$rank ==
                ncol(i_design) && nrow(i_design) > ncol(i_design)) {
                ind_gamma[i_feature, ] <- i_check_batch
            }
            else ind_gamma[i_feature, ] <- FALSE
        }
    }
    ind_gamma[, apply(ind_gamma, 2, sum) < 2] <- FALSE
    ind_feature <- apply(ind_gamma, 1, any)
    if (all(!ind_feature))
        stop("All features are single-batch-specific; MMUPHin cannot perform correction!")
    return(list(ind_data = ind_data, ind_gamma = ind_gamma, ind_mod = ind_mod,
                ind_feature = ind_feature))
}

transform_features <- function (features, transform = "NONE", pseudo_count = 0)
{
    type_features <- check_feature_abd(features)
    if (type_features != "proportions")
        stop("Transformation should only be applied to normalized features",
             " (proportions)")
    features <- features + pseudo_count
    if (transform == "LOG") {
        if (any(features <= 0))
            stop("LOG transformation not applicable to values smaller than or equal",
                 " to zero!")
        features <- apply(features, 2, LOG)
    }
    if (transform == "AST")
        features <- apply(features, 2, AST)
    if (transform == "NONE")
        features <- features
    return(features)
}

standardize_feature <- function (y, i_design, n_batch)
{
    beta_hat <- solve(crossprod(i_design), crossprod(i_design,
                                                     y))
    grand_mean <- mean(i_design[, seq_len(n_batch)] %*% beta_hat[seq_len(n_batch),
    ])
    var_pooled <- stats::var(y - (i_design %*% beta_hat)[, 1])
    if (isTRUE(all.equal(var_pooled, 0)))
        var_pooled <- 1
    stand_mean <- rep(grand_mean, length(y))
    if (ncol(i_design) > n_batch) {
        stand_mean <- stand_mean + (i_design[, -seq_len(n_batch),
                                             drop = FALSE] %*% beta_hat[-seq_len(n_batch), ])[,
                                                                                              1]
    }
    y_stand <- (y - stand_mean)/sqrt(var_pooled)
    return(list(y_stand = y_stand, stand_mean = stand_mean, var_pooled = var_pooled))
}


fit_stand_feature <- function (s_data, design, l_ind)
{
    l_stand_feature <- list()
    for (i_feature in seq_len(nrow(s_data))) {
        if (l_ind$ind_feature[i_feature]) {
            i_design <- design[l_ind$ind_data[i_feature, ], c(l_ind$ind_gamma[i_feature,
            ], l_ind$ind_mod), drop = FALSE]
            if (nrow(i_design) <= 1 | ncol(i_design) <= 1)
                stop("Something wrong happened!")
            stand_fit <- standardize_feature(y = s_data[i_feature,
                                                        l_ind$ind_data[i_feature, ]], i_design = i_design,
                                             n_batch = sum(l_ind$ind_gamma[i_feature, ]))
            s_data[i_feature, l_ind$ind_data[i_feature, ]] <- stand_fit$y_stand
            l_stand_feature[[i_feature]] <- stand_fit
        }
        else l_stand_feature[[i_feature]] <- NULL
    }
    return(list(s_data = s_data, l_stand_feature = l_stand_feature))
}

fit_EB <- function (s_data, l_stand_feature, batchmod, n_batch, l_ind)
{
    if (n_batch != ncol(batchmod))
        stop("n_batch does not agree with batchmod!")
    gamma_hat <- delta_hat <- matrix(NA, nrow = nrow(s_data),
                                     ncol = n_batch)
    for (i_feature in seq_len(nrow(s_data))) {
        if (l_ind$ind_feature[i_feature]) {
            i_batchmod <- batchmod[l_ind$ind_data[i_feature,
            ], l_ind$ind_gamma[i_feature, ], drop = FALSE]
            i_batchmod[i_batchmod == 0] <- NA
            i_s_data_batch <- s_data[i_feature, l_ind$ind_data[i_feature,
            ]] * i_batchmod
            if (nrow(batchmod[l_ind$ind_data[i_feature, ], l_ind$ind_gamma[i_feature,
            ], drop = FALSE]) <= 1 | ncol(batchmod[l_ind$ind_data[i_feature,
            ], l_ind$ind_gamma[i_feature, ], drop = FALSE]) <=
            1)
                stop("Something wrong happened!")
            i_gamma <- apply(i_s_data_batch, 2, mean, na.rm = TRUE)
            i_delta <- apply(i_s_data_batch, 2, stats::sd, na.rm = TRUE)
            i_delta[is.na(i_delta)] <- 1
            i_delta[i_delta == 0] <- 1
            gamma_hat[i_feature, l_ind$ind_gamma[i_feature, ]] <- i_gamma
            delta_hat[i_feature, l_ind$ind_gamma[i_feature, ]] <- i_delta
        }
    }
    gamma_bar <- apply(gamma_hat, 2, mean, na.rm = TRUE)
    t2 <- apply(gamma_hat, 2, stats::var, na.rm = TRUE)
    a_prior <- apply(delta_hat, 2, aprior, na.rm = TRUE)
    b_prior <- apply(delta_hat, 2, bprior, na.rm = TRUE)
    if (any(apply(!is.na(gamma_hat), 2, sum) < 2) | any(apply(!is.na(delta_hat),
                                                              2, sum) < 2))
        stop("One batch has only one feature with valid parameter estimate!")
    return(list(gamma_hat = gamma_hat, delta_hat = delta_hat,
                gamma_bar = gamma_bar, t2 = t2, a_prior = a_prior, b_prior = b_prior))
}

postmean <- function (g_hat, g_bar, n, d_star, t2)
{
    (t2 * n * g_hat + d_star * g_bar)/(t2 * n + d_star)
}

postvar <- function (sum2, n, a, b)
{
    (0.5 * sum2 + b)/(n/2 + a - 1)
}

it_sol <- function (s_data, g_hat, d_hat, g_bar, t2, a, b, control)
{
    n <- rowSums(!is.na(s_data))
    g.old <- g_hat
    d.old <- d_hat
    change <- 1
    count <- 0
    while (change > control$conv) {
        g.new <- postmean(g_hat, g_bar, n, d.old, t2)
        sum2 <- rowSums((s_data - g.new %*% t(rep(1, ncol(s_data))))^2,
                        na.rm = TRUE)
        sum2[sum2 == 0] <- NA
        d.new <- postvar(sum2, n, a, b)
        change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old,
                      na.rm = TRUE)
        g.old <- g.new
        d.old <- d.new
        count <- count + 1
        if (count > control$maxit)
            stop("Maximum iteration reached!")
    }
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g_star", "d_star")
    adjust
}

fit_shrink <- function (s_data, l_params, batchmod, n_batch, l_ind, control)
{
    if (n_batch != ncol(batchmod))
        stop("n_batch does not agree with batchmod!")
    gamma_star <- delta_star <- matrix(NA, nrow = nrow(s_data),
                                       ncol = n_batch)
    results <- lapply(seq_len(n_batch), function(i_batch) {
        i_s_data <- s_data
        i_s_data[!l_ind$ind_data] <- NA
        i_s_data[, !as.logical(batchmod[, i_batch])] <- NA
        i_s_data[!l_ind$ind_gamma[, i_batch], ] <- NA
        temp <- it_sol(s_data = i_s_data, g_hat = l_params$gamma_hat[,
                                                                     i_batch], d_hat = l_params$delta_hat[, i_batch],
                       g_bar = l_params$gamma_bar[i_batch], t2 = l_params$t2[i_batch],
                       a = l_params$a_prior[i_batch], b = l_params$b_prior[i_batch],
                       control = control)
        gamma_star <- temp[1, ]
        delta_star <- temp[2, ]
        list(gamma_star = gamma_star, delta_star = delta_star)
    })
    for (i_batch in seq_len(n_batch)) {
        gamma_star[, i_batch] <- results[[i_batch]]$gamma_star
        delta_star[, i_batch] <- results[[i_batch]]$delta_star
    }
    return(list(gamma_star = gamma_star, delta_star = delta_star))
}

normalize_features <- function (features, normalization = "NONE", pseudo_count = 0)
{
    if (any(features < 0))
        stop("Feature table must be non-negative for normalization!")
    features <- features + pseudo_count
    if (normalization == "TSS")
        features <- apply(features, 2, TSS)
    if (normalization == "NONE")
        features <- features
    return(features)
}

back_transform_abd <- function (adj_data, feature_abd, type_feature_abd)
{
    adj_data <- 2^adj_data
    adj_data[feature_abd == 0] <- 0
    adj_data <- normalize_features(adj_data, normalization = "TSS")
    adj_data <- t(t(adj_data) * apply(feature_abd, 2, sum))
    dimnames(adj_data) <- dimnames(feature_abd)
    if (type_feature_abd == "counts")
        adj_data <- round(adj_data)
    return(adj_data)
}

fill_dimnames <- function (x, row_prefix, col_prefix)
{
    if (!(is.matrix(x) | is.data.frame(x)))
        stop("x must either be a matrix or a data frame!")
    if (missing(row_prefix) | missing(col_prefix))
        stop("Row/column prefixes must be specified!")
    if (is.null(rownames(x)))
        rownames(x) <- paste0(row_prefix, seq_len(nrow(x)))
    if (is.null(colnames(x)))
        colnames(x) <- paste0(col_prefix, seq_len(ncol(x)))
    return(x)
}

diagnostic_adjust_batch <- function (feature_abd, feature_abd_adj, var_batch, gamma_hat,
          gamma_star, output)
{
    feature_abd <- fill_dimnames(feature_abd, "Feature", "Sample")
    dimnames(feature_abd_adj) <- dimnames(feature_abd)
    if (!is.factor(var_batch))
        stop("var_batch should be a factor!")
    df_plot <- data.frame(gamma_hat = as.vector(gamma_hat), gamma_star = as.vector(gamma_star),
                          var_batch = factor(rep(levels(var_batch), each = nrow(gamma_hat)),
                                             levels = levels(var_batch)))
    df_plot <- subset(df_plot, !is.na(gamma_hat), !is.na(gamma_star))
    p_shrinkage <- ggplot2::ggplot(df_plot, ggplot2::aes(x = gamma_hat, y = gamma_star,
                                       color = var_batch)) + ggplot2::geom_point() + ggplot2::geom_abline(intercept = 0,
                                                                                        slope = 1) + ggplot2::ggtitle("Shrinkage of batch mean parameters") +
        ggplot2::theme(legend.position = c(0, 1), legend.justification = c(0,
                                                                  1), legend.direction = "horizontal", legend.background = ggplot2::element_blank(),
              legend.text = ggplot2::element_blank()) + ggplot2::xlab("Gamma") +
        ggplot2::ylab("Gamma (shrinked)")
    mat_ra <- normalize_features(feature_abd, "TSS")
    mat_ra_adj <- normalize_features(feature_abd_adj, "TSS")
    df_mean_batch <- as.data.frame(apply(mat_ra, 1, function(x) tapply(x,
                                                                       var_batch, mean)))
    df_mean_batch_adj <- as.data.frame(apply(mat_ra_adj, 1, function(x) tapply(x,
                                                                               var_batch, mean)))
    colnames(df_mean_batch) <- colnames(df_mean_batch_adj) <- rownames(feature_abd)
    df_mean_batch$batch <- df_mean_batch_adj$batch <- levels(var_batch)
    df_mean_batch$Adjustment <- "Original"
    df_mean_batch_adj$Adjustment <- "Adjusted"
    df_batch <- rbind(df_mean_batch, df_mean_batch_adj)
    df_batch$Adjustment <- factor(df_batch$Adjustment, levels = c("Original",
                                                                  "Adjusted"))
    df_batch <- tidyr::gather(df_batch, key = "Feature", value = "mean_batch",
                              -Adjustment, -batch)
    df_mean_overall <- data.frame(mean_overall = apply(mat_ra,
                                                       1, mean))
    df_mean_overall_adj <- data.frame(mean_overall = df_mean_overall$mean_overall +
                                          max(df_mean_overall$mean_overall)/100)
    df_mean_overall$Feature <- df_mean_overall_adj$Feature <- rownames(feature_abd)
    df_mean_overall$Adjustment <- "Original"
    df_mean_overall_adj$Adjustment <- "Adjusted"
    df_overall <- rbind(df_mean_overall, df_mean_overall_adj)
    df_overall$Adjustment <- factor(df_overall$Adjustment, levels = c("Original",
                                                                      "Adjusted"))
    df_plot <- merge(df_batch, df_overall, by = c("Feature",
                                                  "Adjustment"))
    p_mean <- ggplot2::ggplot(df_plot, ggplot2::aes(x = mean_overall, y = mean_batch)) +
        ggplot2::geom_point(ggplot2::aes(color = Adjustment)) + ggplot2::geom_line(ggplot2::aes(color = Adjustment,
                                                            group = paste0(Feature, Adjustment))) + ggplot2::geom_abline(intercept = 0,
                                                                                                                slope = 1) + ggplot2::scale_color_manual(values = c(Original = "black",
                                                                                                                                                           Adjusted = "red")) + ggplot2::theme(legend.position = c(0, 1),
                                                                                                                                                                                      legend.justification = c(0, 1), legend.direction = "horizontal",
                                                                                                                                                                                      legend.background = ggplot2::element_blank()) + ggplot2::ggtitle("Original/adjusted mean abundance") +
        ggplot2::xlab("Overal mean") + ggplot2::ylab("Batch mean")
    plot <- cowplot::plot_grid(p_shrinkage, p_mean, nrow = 1)
    ggplot2::ggsave(plot = plot, filename = output, device = "pdf", width = 8,
           height = 4, units = "in")
    invisible(plot)
}

match_control <- function (default, control)
{
    if (missing(control))
        control <- list()
    control_pos <- pmatch(names(control), names(default))
    default[c(stats::na.omit(control_pos))] <- control[!is.na(control_pos)]
    return(default)
}

aprior <- function (delta_hat, na.rm = FALSE)
{
    m <- mean(delta_hat, na.rm = na.rm)
    s2 <- var(delta_hat, na.rm = na.rm)
    if (s2 == 0)
        s2 <- 1
    (2 * s2 + m^2)/s2
}

bprior <- function (delta_hat, na.rm = FALSE)
{
    m <- mean(delta_hat, na.rm = na.rm)
    s2 <- var(delta_hat, na.rm = na.rm)
    if (s2 == 0)
        s2 <- 1
    (m * s2 + m^3)/s2
}

LOG <- function (x)
{
    return(log2(x))
}

AST <- function (x)
{
    return(asin(sqrt(x)))
}

TSS <- function (x)
{
    if (all(x == 0))
        return(x)
    return(x/sum(x))
}

adjust_batch_muff <- function (feature_abd, batch, covariates = NULL, data, control)
{
    control <- match_control(default = control_adjust_batch,
                             control = control)
    verbose <- control$verbose
    feature_abd <- as.matrix(feature_abd)
    type_feature_abd <- check_feature_abd(feature_abd = feature_abd)
    if (verbose)
        message("feature_abd is ", type_feature_abd)
    data <- as.data.frame(data)
    samples <- check_samples(feature_abd = feature_abd, data = data)
    if (length(batch) > 1)
        stop("Only one batch variable is supported!")
    df_batch <- check_metadata(data = data, variables = batch)
    df_covariates <- check_metadata(data = data, variables = covariates)
    var_batch <- check_batch(df_batch[[batch]], min_n_batch = 2)
    n_batch <- nlevels(x = var_batch)
    if (verbose)
        message("Found ", n_batch, " batches")
    batchmod <- construct_design(data = df_batch, with_intercept = FALSE)
    mod <- construct_design(data = df_covariates, with_intercept = TRUE)[,
                                                                         -1, drop = FALSE]
    if (!check_rank(design = mod))
        stop("Covariates are confounded!")
    design <- cbind(batchmod, mod)
    if (!check_rank(design = design))
        stop("Covariates and batch are confounded!")
    if (verbose)
        message("Adjusting for ", ncol(mod), " covariate(s) or covariate(s) level(s)")
    if (is.null(control$pseudo_count)) {
        pseudo_count <- set_pseudo(features = feature_abd)
        if (verbose)
            message("Pseudo count is not specified and set to half of minimal ",
                    "non-zero value: ", format(pseudo_count, digits = 3,
                                               scientific = TRUE))
    }
    else pseudo_count <- check_pseudo_count(control$pseudo_count)
    log_data <- transform_features(features = normalize_features(features = feature_abd,
                                                                 normalization = "TSS", pseudo_count = pseudo_count),
                                   transform = "LOG")
    l_ind <- construct_ind(feature_abd = feature_abd, n_batch = n_batch,
                           design = design, zero_inflation = control$zero_inflation)
    batch_no_correct <- apply(!l_ind$ind_gamma, 2, all)
    if (any(batch_no_correct))
        stop(paste0("The following batch(es) either have no present features,",
                    " or\nare confounded with the covariates in", " features that are present.\nPlease remove them",
                    " from the data before running batch correction:\n",
                    paste0(levels(var_batch)[batch_no_correct], collapse = ", ")))
    if (verbose)
        message("Adjusting for (after filtering) ", sum(l_ind$ind_feature),
                " features")
    if (verbose)
        message("Standardizing data across features")
    stand_fit <- fit_stand_feature(s_data = log_data, design = design,
                                   l_ind = l_ind)
    s_data <- stand_fit$s_data
    l_stand_feature <- stand_fit$l_stand_feature
    if (verbose)
        message("Estimating batch difference parameters and EB priors")
    params_fit <- fit_EB(s_data = s_data, l_stand_feature = l_stand_feature,
                         batchmod = batchmod, n_batch = n_batch, l_ind = l_ind)
    if (verbose)
        message("Performing shrinkage adjustments on batch difference parameters")
    params_shrinked <- fit_shrink(s_data = s_data, l_params = params_fit,
                                  batchmod = batchmod, n_batch = n_batch, l_ind = l_ind,
                                  control = control)
    if (verbose)
        message("Performing batch corrections")
    adj_data <- adjust_EB(s_data = s_data, l_params_shrink = params_shrinked,
                          l_stand_feature = l_stand_feature, batchmod = batchmod,
                          n_batch = n_batch, l_ind = l_ind)
    if (any(is.na(adj_data)))
        stop("There are missing values in the adjusted data!")
    feature_abd_adj <- back_transform_abd(adj_data = adj_data,
                                          feature_abd = feature_abd, type_feature_abd = type_feature_abd)
    if (!is.null(control$diagnostic_plot))
        diagnostic_adjust_batch(feature_abd = feature_abd, feature_abd_adj = feature_abd_adj,
                                var_batch = var_batch, gamma_hat = params_fit$gamma_hat,
                                gamma_star = params_shrinked$gamma_star, output = control$diagnostic_plot)
    return(list(feature_abd_adj = feature_abd_adj, control = control))
}

