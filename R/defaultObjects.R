control_adjust_batch <- list(zero_inflation = TRUE, pseudo_count = NULL, diagnostic_plot = NULL,
                             conv = 1e-04, maxit = 1000, verbose = TRUE)
#"adjust_batch_diagnostic.pdf"

control_lm_meta <- list(normalization = "TSS", transform = "AST", analysis_method = "LM",
     rma_method = "REML", rma_conv = 1e-04, rma_maxit = 1000,
     output = tempdir(), forest_plot = "forest.pdf",
     verbose = TRUE)
