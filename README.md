README
================
2026-04-12

# muffinette - Meta-analysis of networks in microbiome studies

The repository houses the **`muffinette`** R package for multi-study
meta-analysis of microbial network properties.

## Dependencies

`muffinette` requires the following `R` package: `devtools` (for
installation only). Please install it before installing `muffinette`,
which can be done as follows (execute from within a fresh R session):

``` r
install.packages("devtools")
library(devtools)
```

## Installation

Once the dependencies are installed, `muffinette` can be loaded using
the following command:

``` r
devtools::install_github("himelmallick/muffinette", quiet = TRUE)
library(muffinette)
```

## Example Implementation

We showcase the effectiveness of our network-connectivity-based
meta-analytic approach (muffinette) by drawing a comparison with
feature-abundance-based meta-analysis (MMUPHin) based on the CRC data
available through the R package `MMUPHin`. This data consists of species
level relative abundance profiles of CRC and control patients in the
five public studies used in Thomas et al. (2019). The dataset is sourced
from the R package `curatedMetagenomicData`.

### Data Pre-processing

``` r
rm(list=ls())
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(MMUPHin)
library(muffinette)
library(vegan)
## Loading required package: permute
## 
## Attaching package: 'permute'
## The following object is masked from 'package:devtools':
## 
##     check
## This is vegan 2.7-3
library(ggplot2)
## Use suppressPackageStartupMessages() to eliminate package startup messages
library(cowplot)
library(ggpubr)
## 
## Attaching package: 'ggpubr'
## The following object is masked from 'package:cowplot':
## 
##     get_legend


data("CRC_abd", "CRC_meta")

set.seed(2310)

data_meta <- data.frame(sampleID = colnames(CRC_abd),
                        age = CRC_meta$age,
                        gender = as.factor(CRC_meta$gender),
                        response = CRC_meta$study_condition,
                        study = as.factor(CRC_meta$studyID))

data_abd <- CRC_abd # feature-by-sample matrix
rownames(data_meta) <- colnames(data_abd)

filtered_featuretable <- muffinette::filter_abdfeatures(x = t(data_abd),
                                                        abd_threshold = 0,
                                                        prev_threshold = 0.1)
filtered_featuretable <- muffinette::filter_varfeatures(x = filtered_featuretable, topV = NULL)

meta_abd_mat <- t(as.matrix(filtered_featuretable)) ## feature-by-sample matrix
```

### Batch (Study) Effect Correction

``` r
batch_corrected_abd <- MMUPHin::adjust_batch(feature_abd = meta_abd_mat,
                                             batch = "study",
                                             covariates = "response",
                                             data = data_meta)$feature_abd_adj
## feature_abd is proportions
## Found 5 batches
## Adjusting for 1 covariate(s) or covariate(s) level(s)
## Pseudo count is not specified and set to half of minimal non-zero value: 5e-08
## Adjusting for (after filtering) 232 features
## Standardizing data across features
## Estimating batch difference parameters and EB priors
## Performing shrinkage adjustments on batch difference parameters
## Performing batch corrections
```

### Visualization of Batch (Study) Effect Correction

``` r
D_before <- vegan::vegdist(t(meta_abd_mat))
D_after <- vegan::vegdist(t(batch_corrected_abd))

fit_adonis_before <- vegan::adonis2(D_before ~ study, data = data_meta)
fit_adonis_after <- vegan::adonis2(D_after ~ study, data = data_meta)

print(fit_adonis_before)
## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999
## 
## vegan::adonis2(formula = D_before ~ study, data = data_meta)
##           Df SumOfSqs      R2      F Pr(>F)    
## Model      4   13.503 0.08202 12.196  0.001 ***
## Residual 546  151.130 0.91798                  
## Total    550  164.633 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
print(fit_adonis_after)
## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999
## 
## vegan::adonis2(formula = D_after ~ study, data = data_meta)
##           Df SumOfSqs      R2      F Pr(>F)    
## Model      4    5.405 0.03385 4.7822  0.001 ***
## Residual 546  154.276 0.96615                  
## Total    550  159.681 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##############
# Ordination #
##############
# Before
R2_before <- round(fit_adonis_before$R2[1]*100, 1)
pcoa_before <- cmdscale(D_before, eig = TRUE)
ord_before <- as.data.frame(pcoa_before$points)
percent_var_before <- round(pcoa_before$eig / sum(pcoa_before$eig) * 100, 1)[1:2]
before_labels <- c(paste('Axis 1 (', percent_var_before[1], '%)', sep = ''), 
                   paste('Axis 2 (', percent_var_before[2], '%)', sep = ''))

before_phrase <- paste('Before (PERMANOVA R2 = ', R2_before, '%)', sep = '')
colnames(ord_before) <- c('PC1', 'PC2')
ord_before$Study <- data_meta$study

# After
R2_after <- round(fit_adonis_after$R2[1]*100, 1)
pcoa_after <- cmdscale(D_after, eig = TRUE)
ord_after <- as.data.frame(pcoa_after$points)
percent_var_after <- round(pcoa_after$eig / sum(pcoa_after$eig) * 100, 1)[1:2]
after_labels <- c(paste('Axis 1 (', percent_var_after[1], '%)', sep = ''), 
                  paste('Axis 2 (', percent_var_after[2], '%)', sep = ''))

after_phrase <- paste('After (PERMANOVA R2 = ', R2_after, '%)', sep = '')
colnames(ord_after) <- c('PC1', 'PC2')
ord_after$Study <- data_meta$study


# Ordination Plot
p_before <- ggplot(ord_before, aes(x = PC1, y = PC2, color = Study)) + 
    geom_point(size = 4) + 
    theme_bw() + 
    xlab(before_labels[1]) + 
    ylab(before_labels[2]) + 
    ggtitle(before_phrase) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position ="none")

p_after <- ggplot(ord_after, aes(x = PC1, y = PC2, color = Study)) + 
    geom_point(size = 4) + 
    theme_bw() + 
    xlab(after_labels[1]) + 
    ylab(after_labels[2]) + 
    ggtitle(after_phrase) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position ="none")

p <- plot_grid(p_before, p_after, ncol = 2)
p
```

![](README_files/figure-gfm/visualize-batch-correction-1.png)<!-- -->

### Meta-analysis of feature abundance data (MMUPHin)

``` r
out_mmuphin <- MMUPHin::lm_meta(feature_abd = batch_corrected_abd,
                            exposure = "response",
                            batch = "study",
                            data = data_meta, control = list(forest_plot = "forest.pdf",
                                                             normalization = 'NONE',
                                                             transform = 'NONE'))
## Found 5 batches
## Fitting Maaslin2 on batch FengQ_2015.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch HanniganGD_2017.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch VogtmannE_2016.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch YuJ_2015.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch ZellerG_2014.metaphlan_bugs_list.stool...
## Fitting meta-analysis model.


metafits_mmuphin <- out_mmuphin$meta_fits
selected_features_mmuphin1 <- metafits_mmuphin[metafits_mmuphin$qval.fdr < 0.05, 1]
selected_features_mmuphin1
## [1] "s__Ruminococcus_torques"                "s__Lachnospiraceae_bacterium_7_1_58FAA"
selected_features_mmuphin2 <- metafits_mmuphin[metafits_mmuphin$qval.fdr < 0.25, 1]
selected_features_mmuphin2
##  [1] "s__Bacteroides_fragilis"                "s__Ruminococcus_torques"               
##  [3] "s__Alistipes_onderdonkii"               "s__Streptococcus_thermophilus"         
##  [5] "s__Lachnospiraceae_bacterium_1_1_57FAA" "s__Clostridium_sp_L2_50"               
##  [7] "s__Bifidobacterium_catenulatum"         "s__Adlercreutzia_equolifaciens"        
##  [9] "s__Solobacterium_moorei"                "s__Fusobacterium_nucleatum"            
## [11] "s__Lachnospiraceae_bacterium_7_1_58FAA"
```

### Meta-network-analysis (muffinette)

``` r
batchvar <- data_meta$study
exposurevar <- data_meta$response

out_muffinette <- muffinette::muffinette(metaAbd = data_abd,
                                batchvar = batchvar, exposurevar = exposurevar,
                                metaData = data_meta,
                                filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
                                comp = TRUE, net.est.method = "SparCC",
                                covariates = NULL, ncores = 4, verbose = TRUE, fixseed = 2310,
                                control = list(network = list(iter = 20, inner_iter = 10, th = 0.1)))
## Registered S3 methods overwritten by 'huge':
##   method    from
##   plot.roc  pROC
##   plot.sim  lava
##   print.roc pROC
##   print.sim lava
## extractScore completed for n = 107
## Network estimated for study 1 / 5
## extractScore completed for n = 55
## Network estimated for study 2 / 5
## extractScore completed for n = 104
## Network estimated for study 3 / 5
## extractScore completed for n = 128
## Network estimated for study 4 / 5
## extractScore completed for n = 157
## Network estimated for study 5 / 5
## Registered S3 methods overwritten by 'GenomeInfoDb':
##   method                from   
##   as.data.frame.Seqinfo Seqinfo
##   merge.Seqinfo         Seqinfo
##   summary.Seqinfo       Seqinfo
##   update.Seqinfo        Seqinfo
## Found5batches
## Adjusting for1covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
## Found 5 batches
## Fitting Maaslin2 on batch FengQ_2015.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch HanniganGD_2017.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch VogtmannE_2016.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch YuJ_2015.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch ZellerG_2014.metaphlan_bugs_list.stool...
## Fitting meta-analysis model.
## expo: exposure
## value_expo: CRC
## Meta-analysis done...

metafits_muffinette <- out_muffinette$metafits
selected_features_muffinette1 <- metafits_muffinette[metafits_muffinette$qval < 0.05, 1]
selected_features_muffinette1
## character(0)
selected_features_muffinette2 <- metafits_muffinette[metafits_muffinette$qval < 0.25, 1]
selected_features_muffinette2
## character(0)
```

``` r
batchvar <- data_meta$study
exposurevar <- data_meta$response

out_muffinette_group <- muffinette::muffinette_group(metaAbd = data_abd,
                                batchvar = batchvar, exposurevar = exposurevar,
                                metaData = data_meta,
                                filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = NULL,
                                comp = TRUE, net.est.method = "SparCC",
                                covariates = NULL, ncores = 4, verbose = TRUE, fixseed = 2310,
                                control = list(network = list(iter = 20, inner_iter = 10, th = 0.1)))
## extractScore completed for n = 46
## extractScore completed for n = 61
## Network estimated for study 1 / 5
## extractScore completed for n = 27
## extractScore completed for n = 28
## Network estimated for study 2 / 5
## extractScore completed for n = 52
## extractScore completed for n = 52
## Network estimated for study 3 / 5
## extractScore completed for n = 75
## extractScore completed for n = 53
## Network estimated for study 4 / 5
## extractScore completed for n = 91
## extractScore completed for n = 66
## Network estimated for study 5 / 5
## Found5batches
## Adjusting for1covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
## Found 5 batches
## Fitting Maaslin2 on batch FengQ_2015.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch HanniganGD_2017.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch VogtmannE_2016.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch YuJ_2015.metaphlan_bugs_list.stool...
## Fitting Maaslin2 on batch ZellerG_2014.metaphlan_bugs_list.stool...
## Fitting meta-analysis model.
## expo: exposure
## value_expo: CRC
## Meta-analysis done...

metafits_muffinette_group <- out_muffinette_group$metafits
selected_features_muffinette_group1 <- metafits_muffinette_group[metafits_muffinette_group$qval < 0.05, 1]
selected_features_muffinette_group1
##  [1] "s__Clostridium_innocuum"                           
##  [2] "s__Flavonifractor_plautii"                         
##  [3] "s__Gemella_sanguinis"                              
##  [4] "s__candidate_division_TM7_single_cell_isolate_TM7b"
##  [5] "s__Subdoligranulum_unclassified"                   
##  [6] "s__Butyrivibrio_unclassified"                      
##  [7] "s__Bacteroides_ovatus"                             
##  [8] "s__Lachnospiraceae_bacterium_1_1_57FAA"            
##  [9] "s__Oribacterium_sinus"                             
## [10] "s__Coprococcus_comes"                              
## [11] "s__Erysipelotrichaceae_bacterium_2_2_44A"          
## [12] "s__Pseudomonas_unclassified"                       
## [13] "s__Parabacteroides_distasonis"                     
## [14] "s__Coprococcus_catus"                              
## [15] "s__Eggerthella_unclassified"                       
## [16] "s__Subdoligranulum_sp_4_3_54A2FAA"                 
## [17] "s__Bifidobacterium_bifidum"                        
## [18] "s__Veillonella_unclassified"                       
## [19] "s__Ruminococcus_bromii"                            
## [20] "s__Rothia_mucilaginosa"                            
## [21] "s__Erysipelotrichaceae_bacterium_5_2_54FAA"
selected_features_muffinette_group2 <- metafits_muffinette_group[metafits_muffinette_group$qval < 0.25, 1]
selected_features_muffinette_group2
##  [1] "s__Clostridium_innocuum"                           
##  [2] "s__Flavonifractor_plautii"                         
##  [3] "s__Gemella_sanguinis"                              
##  [4] "s__candidate_division_TM7_single_cell_isolate_TM7b"
##  [5] "s__Subdoligranulum_unclassified"                   
##  [6] "s__Butyrivibrio_unclassified"                      
##  [7] "s__Bacteroides_ovatus"                             
##  [8] "s__Lachnospiraceae_bacterium_1_1_57FAA"            
##  [9] "s__Oribacterium_sinus"                             
## [10] "s__Coprococcus_comes"                              
## [11] "s__Erysipelotrichaceae_bacterium_2_2_44A"          
## [12] "s__Pseudomonas_unclassified"                       
## [13] "s__Parabacteroides_distasonis"                     
## [14] "s__Coprococcus_catus"                              
## [15] "s__Eggerthella_unclassified"                       
## [16] "s__Subdoligranulum_sp_4_3_54A2FAA"                 
## [17] "s__Bifidobacterium_bifidum"                        
## [18] "s__Veillonella_unclassified"                       
## [19] "s__Ruminococcus_bromii"                            
## [20] "s__Rothia_mucilaginosa"                            
## [21] "s__Erysipelotrichaceae_bacterium_5_2_54FAA"        
## [22] "s__Dorea_unclassified"                             
## [23] "s__Lachnospiraceae_bacterium_3_1_57FAA_CT1"        
## [24] "s__Bacteroides_vulgatus"                           
## [25] "s__Clostridium_sp_L2_50"                           
## [26] "s__Paraprevotella_clara"                           
## [27] "s__Lachnospiraceae_bacterium_5_1_63FAA"            
## [28] "s__Bifidobacterium_catenulatum"                    
## [29] "s__Blautia_hydrogenotrophica"                      
## [30] "s__Erysipelotrichaceae_bacterium_3_1_53"           
## [31] "s__Coriobacteriaceae_bacterium_phI"                
## [32] "s__Subdoligranulum_variabile"                      
## [33] "s__Bifidobacterium_adolescentis"                   
## [34] "s__Bacteroides_coprocola"                          
## [35] "s__Lactobacillus_delbrueckii"                      
## [36] "s__Paraprevotella_xylaniphila"                     
## [37] "s__Lactococcus_lactis"                             
## [38] "s__C2likevirus_unclassified"                       
## [39] "s__Alistipes_onderdonkii"                          
## [40] "s__Oscillibacter_unclassified"                     
## [41] "s__Clostridium_bartlettii"                         
## [42] "s__Roseburia_hominis"                              
## [43] "s__Lachnospiraceae_bacterium_2_1_58FAA"
```

### Comparison of the Three Meta-Analysis Approaches (threshold as 0.05)

![](README_files/figure-gfm/compare0.05-1.png)<!-- -->

### Comparison of the Three Meta-Analysis Approaches (threshold as 0.25)

![](README_files/figure-gfm/compare0.25-1.png)<!-- -->
