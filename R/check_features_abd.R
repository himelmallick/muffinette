#' @title check_features_abd
#'
#' @description function to make sure that feature abundances are neither missing nor negative. 
#' @param feature_abd a sample-by-feature data frame.
#' @return error messages if feature abundances are (1) missing, (2) negative, (3) not proportions.
#' @export

############################################
######## check features abundances #########
############################################

check_features_abd <- function(feature_abd) {
    if(any(is.na(feature_abd))) {
        stop("Found missing values in feature abundance table!")
    }
    
    if(any(feature_abd < 0)) {
        stop("Feature abundances cannot be negative!")
    }
    
    if(as.numeric(all(feature_abd <= 1)) == 0) {
        stop("Feature abundances must be proportions!")
    } #all(feature_abd == floor(feature_abd)) |
}
