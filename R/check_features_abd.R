#' @title check_features_abd
#'
#' @description function to ensure that feature abundances are neither missing nor negative.
#' @param feature_abd a sample-by-feature data frame.
#' @return error messages if feature abundances are (1) missing, (2) negative, (3) not proportions.
#' @export

check_features_abd <- function(feature_abd) {
    # Subset data
    feature_abd_cleaned <- feature_abd
    feature_abd_cleaned[is.na(feature_abd_cleaned) | is.nan(feature_abd_cleaned)] <- 0

    if(any(feature_abd_cleaned < 0)) {
        stop("Feature abundances cannot be negative!")
    }

    if(!all(feature_abd_cleaned <= 1)) {
        stop("Feature abundances must be proportions!")
    }

    feature_abd_cleaned
}
