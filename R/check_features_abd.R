#' @title check_features_abd
#'
#' @description function to make sure that feature abundances are neither missing nor negative. 
#' @param feature_abd a sample-by-feature data frame.
#' @return error messages if feature abundances are (1) missing, (2) negative, (3) not proportions.
#' @export

############################################
######## check features abundances #########
############################################

# check_features_abd <- function(feature_abd) {
#     if(any(is.na(feature_abd))) {
#         stop("Found missing values in feature abundance table!")
#     }
#     
#     if(any(feature_abd < 0)) {
#         stop("Feature abundances cannot be negative!")
#     }
#     
#     if(as.numeric(all(feature_abd <= 1)) == 0) {
#         stop("Feature abundances must be proportions!")
#     } #all(feature_abd == floor(feature_abd)) |
# }


# check_features_abd <- function(feature_abd) {
#     # Remove rows and columns with any NA or NaN
#     rows_to_remove <- apply(feature_abd, 1, function(x) any(is.na(x) | is.nan(x)))
#     cols_to_remove <- apply(feature_abd, 2, function(x) any(is.na(x) | is.nan(x)))
#     
#     feature_abd_cleaned <- feature_abd[!rows_to_remove, !cols_to_remove, drop=FALSE]
#     
#     if(nrow(feature_abd_cleaned) == 0 || ncol(feature_abd_cleaned) == 0) {
#         stop("After removing missing values, no data remains!")
#     }
#     
#     if(any(feature_abd_cleaned < 0)) {
#         stop("Feature abundances cannot be negative!")
#     }
#     
#     if(!all(feature_abd_cleaned <= 1)) {
#         stop("Feature abundances must be proportions!")
#     }
#     
#     # Return cleaned data frame
#     feature_abd_cleaned
# }


check_features_abd <- function(feature_abd) {
    # Identify rows to remove
    rows_to_remove <- apply(feature_abd, 1, function(x) any(is.na(x) | is.nan(x)))
    
    # Capture row names
    removed_rownames <- rownames(feature_abd)[rows_to_remove]
    if(is.null(removed_rownames)) {
        removed_rownames <- which(rows_to_remove)
    }
    
    # Subset data
    feature_abd_cleaned <- feature_abd[!rows_to_remove, , drop=FALSE]
    
    # Check if empty
    if(nrow(feature_abd_cleaned) == 0) {
        stop("After removing missing values, no data remains!")
    }
    
    # Other checks
    if(any(feature_abd_cleaned < 0)) {
        stop("Feature abundances cannot be negative!")
    }
    
    if(!all(feature_abd_cleaned <= 1)) {
        stop("Feature abundances must be proportions!")
    }
    
    # Return both cleaned data and removed row names
    list(feature_abd_cleaned = feature_abd_cleaned,
         removed_rownames = removed_rownames)
}
