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
library(MMUPHin)
library(muffinette)
library(vegan)
library(ggplot2)
library(cowplot)
library(ggpubr)


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
filtered_featuretable <- muffinette::filter_varfeatures(x = filtered_featuretable, topV = 100)

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
## Pseudo count is not specified and set to half of minimal non-zero value: 1e-07
## Adjusting for (after filtering) 100 features
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
## Model      4   13.454 0.08281 12.324  0.001 ***
## Residual 546  149.013 0.91719                  
## Total    550  162.467 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
print(fit_adonis_after)
## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999
## 
## vegan::adonis2(formula = D_after ~ study, data = data_meta)
##           Df SumOfSqs      R2      F Pr(>F)    
## Model      4    5.306 0.03366 4.7541  0.001 ***
## Residual 546  152.344 0.96634                  
## Total    550  157.650 1.00000                  
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
selected_features_mmuphin <- metafits_mmuphin[metafits_mmuphin$qval.fdr < 0.05, 1]
selected_features_mmuphin
## [1] "s__Ruminococcus_torques"
```

### Meta-network-analysis (muffinette)

``` r
batchvar <- data_meta$study
exposurevar <- data_meta$response

out_muffinette <- muffinette::muffinette(metaAbd = data_abd,
                                batchvar = batchvar, exposurevar = exposurevar,
                                metaData = data_meta,
                                filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = 100,
                                comp = TRUE, net.est.method = "SparCC",
                                covariates = NULL, ncores = 4, verbose = TRUE, fixseed = 2310,
                                control = list(network = list(iter = 20, inner_iter = 10, th = 0.1)))
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
selected_features_muffinette <- metafits_muffinette[metafits_muffinette$qval < 0.05, 1]
selected_features_muffinette
## character(0)
```

``` r
batchvar <- data_meta$study
exposurevar <- data_meta$response

out_muffinette_group <- muffinette::muffinette_group(metaAbd = data_abd,
                                batchvar = batchvar, exposurevar = exposurevar,
                                metaData = data_meta,
                                filter = TRUE, abd_threshold = 0, prev_threshold = 0.1, topfeatures = 100,
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
selected_features_muffinette_group <- metafits_muffinette_group[metafits_muffinette_group$qval < 0.05, 1]
selected_features_muffinette_group
## [1] "s__Bacteroides_intestinalis" "s__Bacteroides_plebeius"     "s__Megasphaera_unclassified"
## [4] "s__Bacteroides_eggerthii"    "s__Prevotella_copri"         "s__Butyrivibrio_crossotus"  
## [7] "s__Megamonas_hypermegale"    "s__Clostridium_sp_L2_50"
```

### Comparison of the Three Meta-Analysis Approaches

``` r
fig_mmuphin <- metafits_mmuphin %>%
  filter(qval.fdr < 0.05) %>%
  mutate(feature = sub("^species:", "", feature)) %>%
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>%
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity", fill = "lightpink") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(y = NULL, x = NULL,
       title = "MMUPHin")

fig_muffinette <- metafits_muffinette %>%
  filter(qval < 0.05) %>%
  mutate(feature = sub("^species:", "", feature)) %>%
  ggplot(aes(y = coef, x = reorder(feature, coef))) +
  geom_bar(stat = "identity", fill = "lightpink") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 5)) +
  labs(y = NULL, x = NULL, title = "muffinette")

fig <- ggarrange(fig_mmuphin, fig_muffinette, nrow = 2, heights = c(0.75, 1))
fig <- annotate_figure(fig,
                       left = text_grob("Differentially Abundant Species",
                                        size = 8, face = "bold", rot = 90),
                       bottom = text_grob("Estimated Coefficients",
                                          size = 8, hjust = 0.5, face = "bold"))
fig
```
