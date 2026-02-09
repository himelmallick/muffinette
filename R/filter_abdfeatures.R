#' @title filter_abdfeatures
#'
#' @description function to remove features according to user-specified abundance and prevalence threshold values.
#' @param x a sample-by-feature data frame.
#' @param abd_threshold threshold for abundance of features.
#' @param prev_threshold threshold for prevalence of features.
#' @return a sample-by-feature data frame.
#' @export

filter_abdfeatures <- function(x, abd_threshold, prev_threshold) {
    # x is a data frame with features as variables/columns
    x <- as.data.frame(x)
    x <- x[, colMeans(x > abd_threshold) > prev_threshold]
    
    x
}

