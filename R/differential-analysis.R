#' Run differential abundance analysis
#'
#' Performs DESeq2-based differential abundance analysis on SGE screen data,
#' using control oligos for size factor normalization.
#'
#' @param count_matrix Numeric matrix of variant counts with sequences as
#'   row names and samples as columns.
#' @param normalization_matrix Matrix of neutral variants for size factor
#'   estimation (e.g., synonymous/intronic variants).
#' @param sample_metadata Data frame with sample information. Must contain
#'   `condition` and `supplier_name` columns.
#' @param condition_levels Character vector specifying condition factor levels
#'   in order (default: c("Day4", "Day7", "Day10", "Day15")).
#' @param shrinkage_type Type of LFC shrinkage: "normal", "apeglm", or "ashr"
#'   (default: "normal").
#' @param alpha Significance threshold for adjusted p-values (default: 0.05).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{results}{DESeq2 results table as a tibble with SEQUENCE column}
#'     \item{rlog}{Regularized log-transformed DESeqDataSet object}
#'     \item{contrast_summary}{Combined summary across all contrasts including
#'       log2 fold changes, p-values, and z-scores for each sample}
#'   }
#'
#' @details
#' This function implements the full SGE differential abundance pipeline:
#' \enumerate{
#'   \item Extracts condition information from metadata for samples in matrix
#'   \item Estimates size factors from neutral variants (synonymous/intronic)
#'   \item Creates DESeqDataSet with proper experimental design
#'   \item Applies control-derived size factors to all variants
#'   \item Runs DESeq2 Wald tests for all condition contrasts
#'   \item Computes regularized log transformation for visualization
#'   \item Generates z-scores across samples for each variant
#' }
#'
#' The control-based normalization strategy accounts for differences in
#' library size while avoiding bias from variants under selection.
#'
#' @examples
#' \dontrun{
#' results <- run_differential_analysis(
#'   count_matrix = counts,
#'   normalization_matrix = norm_counts,
#'   sample_metadata = metadata
#' )
#'
#' # Access components
#' deseq_results <- results$results
#' rlog_data <- results$rlog
#' summary_table <- results$contrast_summary
#' }
#'
#' @seealso [estimate_size_factors()], [create_normalization_matrix()]
#' @export
run_differential_analysis <- function(count_matrix,
                                      normalization_matrix,
                                      sample_metadata,
                                      condition_levels = c("Day4", "Day7",
                                                           "Day10", "Day15"),
                                      shrinkage_type = "normal",
                                      alpha = 0.05) {

  # Build condition data for samples in the count matrix
 condition_data <- data.frame(
    condition = sample_metadata$condition[
      match(colnames(count_matrix), sample_metadata$supplier_name)
    ]
  )
  condition_data$condition <- factor(
    condition_data$condition,
    levels = condition_levels
  )
  rownames(condition_data) <- colnames(count_matrix)

  # Estimate size factors from control variants
  size_factors <- compute_control_size_factors(
    count_matrix = normalization_matrix,
    sample_metadata = condition_data
  )

  # Create DESeq dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = condition_data,
    design = ~ condition
  )

  # Apply control-derived size factors
  dds <- DESeq2::estimateSizeFactors(dds)
  DESeq2::sizeFactors(dds) <- size_factors

  # Run DESeq2
  dds <- DESeq2::DESeq(dds)

  # Regularized log transform for visualization
  rld <- DESeq2::rlog(dds)

  # Base results
  res <- DESeq2::results(dds) |>
    tibble::as_tibble(rownames = "SEQUENCE")

  # Calculate z-scores
  z_score <- SummarizedExperiment::assay(rld) |>
    as.matrix() |>
    t() |>
    scale() |>
    t() |>
    tibble::as_tibble(rownames = "SEQUENCE")

  colnames(z_score) <- paste0(colnames(z_score), "_z_score")
  colnames(z_score)[1] <- "SEQUENCE"

  # Get all contrasts
  available_contrasts <- DESeq2::resultsNames(dds)
  available_contrasts <- available_contrasts[2:length(available_contrasts)]

  # Build multi-contrast summary using DEGreport
  table_wald <- DEGreport::degComps(
    dds,
    combs = "condition",
    contrast = available_contrasts,
    alpha = alpha,
    skip = FALSE,
    type = shrinkage_type,
    pairs = FALSE,
    fdr = "default"
  )

  first_contrast <- names(table_wald)[1]

  all_contrast_summary <- purrr::imap(table_wald, .create_deg_table) |>
    purrr::reduce(dplyr::left_join, by = "SEQUENCE") |>
    dplyr::rename(baseMean = paste0("baseMean_", first_contrast)) |>
    dplyr::select(-dplyr::starts_with("baseMean_"))

  # Add z-scores
  all_contrast_summary <- all_contrast_summary |>
    dplyr::left_join(z_score, by = "SEQUENCE")

  # Return named list
  list(
    results = res,
    rlog = rld,
    contrast_summary = all_contrast_summary
  )
}


#' Create DEG summary table for a single contrast
#'
#' Internal function to extract differential expression results for one
#' contrast and format with suffixed column names.
#'
#' @param table DEGComps result object.
#' @param name Contrast name to use as column suffix.
#'
#' @return A tibble with contrast-specific column names.
#'
#' @keywords internal
.create_deg_table <- function(table, name) {
  DEGreport::deg(table, "raw") |>
    tibble::as_tibble(rownames = "SEQUENCE") |>
    dplyr::rename_with(~ paste0(.x, "_", name), .cols = -"SEQUENCE")
}
