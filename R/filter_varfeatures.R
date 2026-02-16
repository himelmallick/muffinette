#' @title filter_varfeatures
#'
#' @description function to remove features with zero variance and keep a fixed number of the most variable features from the rest.
#' @param x a sample-by-feature data frame.
#' @param topV number of most variable features to extract. Default is same as the number of features.
#' @return a sample-by-feature data frame.
#' @importFrom stats var
#' @export


filter_varfeatures <- function(x, topV = NULL) {
    # x is a data frame with features as variables/columns

    nzv_x <- caret::nearZeroVar(x, names = TRUE)
    x <- x[, setdiff(names(x), nzv_x)]

    if(is.null(topV)) {
        topV <- ncol(x)
    }
    topV <- min(topV, ncol(x))
    vars <- apply(x, 2, var, na.rm = TRUE)
    top_idx <- order(vars, decreasing = TRUE)[1:topV]
    x <- x[, top_idx, drop = FALSE]

    x
}



