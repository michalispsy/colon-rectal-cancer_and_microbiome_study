#!/usr/bin/env Rscript

# Usage: mmuphin_wrapper.R <feature_matrix_csv> <metadata_csv> <output_csv> <covariates_comma_separated>

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: mmuphin_wrapper.R <feature_matrix_csv> <metadata_csv> <output_csv> <covariates>", call.=FALSE)
}

# Suppress warnings for clean output
suppressMessages(library(MMUPHin))
suppressMessages(library(dplyr))

feat_file = args[1]
meta_file = args[2]
out_file = args[3]
covariates_str = args[4]

covariates_list = unlist(strsplit(covariates_str, ","))

# Read data
# feat: rows are samples, cols are features
feat_df = read.csv(feat_file, row.names=1)
# meta: rows are samples
meta_df = read.csv(meta_file, row.names=1)

# MMUPHin expects features as rows, samples as columns for adjust_batch
feat_mat = t(as.matrix(feat_df))

# adjust_batch using Study as batch
fit <- adjust_batch(feature_abd = feat_mat,
                    batch = "Study",
                    covariates = covariates_list,
                    data = meta_df,
                    control = list(verbose = FALSE))

# Transpose back so rows are samples, cols are features
corrected_mat = t(fit$feature_abd_adj)

# Write output
write.csv(corrected_mat, file=out_file, row.names=TRUE)
