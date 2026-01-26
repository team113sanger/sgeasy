#' sgeasy: Downstream Analysis for Saturation Genome Editing Screens
#'
#' Provides tools for analyzing saturation genome editing (SGE) screen data,
#' including normalization using control oligos, differential abundance
#' analysis with DESeq2, and visualization utilities.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{run_differential_analysis}}}{Run complete differential
#'     abundance analysis pipeline}
#'   \item{\code{\link{create_normalization_matrix}}}{Create matrix of neutral
#'     variants for normalization}
#'   \item{\code{\link{create_count_matrix}}}{Create count matrix from
#'     annotated data}
#'   \item{\code{\link{filter_by_counts}}}{Filter variants by count thresholds}
#'   \item{\code{\link{estimate_size_factors}}}{Estimate DESeq2 size factors}
#'   \item{\code{\link{plot_sample_scatter}}}{Create sample comparison plots}
#'   \item{\code{\link{remove_artefacts}}}{Remove artefact variants}
#'   \item{\code{\link{reweight_replicated_variants}}}{Combine replicated
#'     variants using inverse variance weighting}
#'   \item{\code{\link{sample_distance_matrix}}}{Create sample distance heatmap}
#'   \item{\code{\link{build_sample_pca}}}{Build PCA plot from rlog data}
#' }
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors sizeFactors
#' @importFrom DESeq2 DESeq rlog results resultsNames counts plotPCA
#' @importFrom SummarizedExperiment assay
#' @importFrom DEGreport degComps deg
#' @importFrom dplyr filter select mutate rename rename_with left_join
#' @importFrom dplyr group_by ungroup slice_head across rowwise between n
#' @importFrom dplyr starts_with contains all_of summarise case_when
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom purrr map map2 imap reduce iwalk list_flatten
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme element_text
#' @importFrom ggplot2 xlab ylab coord_fixed theme_classic
#' @importFrom GGally ggpairs wrap
#' @importFrom rlang .data :=
#' @importFrom stats as.formula relevel pnorm p.adjust dist
#' @importFrom grDevices png dev.off colorRampPalette
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cli cli_h1 cli_ul cli_text
#' @import S7
NULL

# Suppress R CMD check notes for NSE variables
utils::globalVariables(c(
  "Consequence", "SEQUENCE", "total_counts", "condition",
  "baseMean", "supplier_name", ".",
  # New variables for reweighting functions
  "HGVSc", "HGVSp", "pam_mut_sgrna_id", "lfcSE_continuous",
  "adj_lfc_continuous", "sgRNA_id", "weight", "weighted_lfc",
  "has_non_pam", "n_observations", "total_weight", "sum_weighted_lfc",
  "combined_lfc", "combined_SE", "combined_Z", "pval",
  # PCA variables
  "PC1", "PC2"
))
