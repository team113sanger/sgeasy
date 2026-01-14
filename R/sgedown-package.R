#' sgedown: Downstream Analysis for Saturation Genome Editing Screens
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
#' }
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors sizeFactors
#' @importFrom DESeq2 DESeq rlog results resultsNames counts
#' @importFrom SummarizedExperiment assay
#' @importFrom DEGreport degComps deg
#' @importFrom dplyr filter select mutate rename rename_with left_join
#' @importFrom dplyr group_by ungroup slice_head across rowwise between
#' @importFrom dplyr starts_with contains all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom purrr map map2 imap reduce
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme element_text
#' @importFrom GGally ggpairs wrap
#' @importFrom rlang .data
#' @importFrom stats as.formula relevel
NULL

# Suppress R CMD check notes for NSE variables
utils::globalVariables(c(

  "Consequence", "SEQUENCE", "total_counts", "condition",
  "baseMean", "supplier_name", "."
))
