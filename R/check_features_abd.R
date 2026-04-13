#' @title check_features_abd
#'
#' @description function to ensure that feature abundances are neither missing nor negative. Missing abundances are replaced by zeros.
#' @param feature_abd a feature-by-sample data frame.
#' @param comp logical. If FALSE, the abundances are assumed non-compositional.
#' @return error messages if feature abundances are (1) negative, and converts to proportions if comp is TRUE yet the abundances are not proportions.
#' @export

check_features_abd <- function(feature_abd, comp) {
    # Subset data
    feature_abd_cleaned <- feature_abd
    feature_abd_cleaned[is.na(feature_abd_cleaned) | is.nan(feature_abd_cleaned)] <- 0

    if(any(feature_abd_cleaned < 0)) {
        stop("Feature abundances cannot be negative!")
    }

    # if(!all(feature_abd_cleaned <= 1)) {
    #     stop("Feature abundances must be proportions!")
    # }
    if(comp) {
        cs <- colSums(feature_abd_cleaned)
        if(!all(abs(cs - 1) < 1e-8)) {
            message("comp is TRUE: Converting to proportions")
            feature_abd_cleaned <- feature_abd_cleaned / cs
        }
    }


    feature_abd_cleaned
}
