#' @title sohpie_jk
#'
#' @description function to compute standardized jackknifes using sparCC. 
#' @param OTUdat sample-by-feature abundance table after batch effect correction.
#' @param groupA Indices of the subjects in the first category of the exposure variable.
#' @param groupB Indices of the subjects in the second category of the exposure variable.
#' @param method method of network estimation.
#' @param seed seed for reproducibility
#' @return sample-by-feature matrix of standardized jackknife values.
#' @import robustbase
#' @import parallel
#' @import dplyr
#' @import gtools
#' @import fdrtool
#' @export


############################################
### Standardized jackknifes using sparCC ###
############################################

sohpie_jk <- function(OTUdat, groupA, groupB, method, seed) {
    
    set.seed(seed)
    
    ##########################################################################
    #STEP 1. Estimate an association matrix via SparCC from the microbiome data.
    ##########################################################################
    if(method == "sparcc") {
        OTUtabA <- OTUdat[groupA, ]
        OTUtabB <- OTUdat[groupB, ]
        n_A <- length(groupA) # Sample size for Group A
        n_B <- length(groupB) # Sample size for Group B
        
        sparcc.matA <- sparcc(OTUtabA)$cor.w
        sparcc.matB <- sparcc(OTUtabB)$cor.w
        
        ##########################################################################
        #STEP 2. Calculate \hat{\theta}_{k} by taking the marginal sum of the
        #association matrix to obtain the network centrality measure for each taxon.
        ##########################################################################
        thetahat_grpA <- SOHPIE::thetahats(sparcc.matA)
        thetahat_grpB <- SOHPIE::thetahats(sparcc.matB)
        
        ##########################################################################
        #STEP 3. Re-estimate association matrix without i-th subject.
        #Then, calculate \hat{\theta}_{k(i)} from the re-estimated association matrix.
        ##########################################################################
        #Re-estimation part
        #Note: You can specify the core size using mc.cores option within mclapply().
        sparcc.mat_drop_grpA <- parallel::mclapply(groupA, function(j) sparcc(OTUtabA[-j, ])$cor.w)
        sparcc.mat_drop_grpB <- parallel::mclapply(groupB, function(j) sparcc(OTUtabB[-j, ])$cor.w)
        # thetahat_{-i} for each taxa
        thetahat_drop_grpA <- sapply(sparcc.mat_drop_grpA, thetahats)
        thetahat_drop_grpB <- sapply(sparcc.mat_drop_grpB, thetahats)
    }
    
    
    ##########################################################################
    #STEP 4. Calculate jackknife pseudovalues (\tilde{\theta}_{ik} using \hat{\theta}_{k} and \hat{\theta}_{k(i)}.
    ########################################################################
    thetatildefun <- function(thetahatinput, thetahatdropinput, sizegroup) {
        thetatildeout <- matrix(NA, ncol=length(thetahatinput), nrow=sizegroup)
        thetatildeout <- sapply(1:nrow(thetahatdropinput), function(k) {
            sizegroup * thetahatinput[k] - (sizegroup - 1) * thetahatdropinput[k, ]
        })
        return(thetatildeout)
    }
    
    thetatilde_grpA <- thetatildefun(thetahat_grpA, thetahat_drop_grpA, n_A)
    thetatilde_grpB <- thetatildefun(thetahat_grpB, thetahat_drop_grpB, n_B)
    thetatilde <- matrix(NA, nrow(OTUdat), ncol(OTUdat))
    thetatilde[groupA, ] <- thetatilde_grpA
    thetatilde[groupB, ] <- thetatilde_grpB
    thetatilde_std <- scale(thetatilde, center = TRUE, scale = TRUE)
    colnames(thetatilde_std) <- colnames(OTUdat) # Map the column names (feature names)
    
    thetatilde_std
}
