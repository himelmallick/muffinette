#' @title filter_features
#'
#' @description function to remove features as per provided abundance and prevalence threshold values,
#'  and zero variance features. 
#' @param x a sample-by-feature data frame.
#' @param abd_threshold threshold for abundance of features.
#' @param prev_threshold threshold for prevalence of features.
#' @return sample-by-feature data frame after removing features.
#' @import caret
#' @export


############################################
############# filter features ##############
############################################


filter_features <- function(x, abd_threshold, prev_threshold) {
    # x is a data frame with features as variables/columns
    x <- as.data.frame(x)
    x <- x[, colMeans(x > abd_threshold) > prev_threshold]
    
    nzv_x <- caret::nearZeroVar(x, names = TRUE)
    x <- x[, setdiff(names(x), nzv_x)] 
    return(x)
}
