#' Estimate size factors from control data
#'
#' Estimates DESeq2 size factors using a subset of data (typically control
#' oligos or neutral variants) for normalization across samples.
#'
#' @param count_data A matrix of raw counts with genes/variants as rows and
#'   samples as columns.
#' @param col_data A data frame with sample metadata. Must contain a
#'   `condition` column.
#' @param design A formula string specifying the experimental design
#'   (default: "~ condition").
#' @param min_row_sum Minimum total count across samples to retain a variant
#'   (default: 10).
#' @param reference_level Reference level for condition factor (default: NULL).
#'
#' @return A named numeric vector of size factors, one per sample.
#'
#' @details
#' This function creates a DESeqDataSet from the provided count matrix,
#' filters low-count variants, and estimates size factors using DESeq2's
#' median-of-ratios method. The size factors can then be applied to a
#' separate count matrix for normalization.
#'
#' @examples
#' \dontrun{
#' # Using synonymous variants for normalization
#' size_factors <- estimate_size_factors(
#'   count_data = norm_matrix,
#'   col_data = sample_info,
#'   reference_level = "Day4"
#' )
#' }
#'
#' @seealso \code{\link[DESeq2]{estimateSizeFactors}}
#' @export
estimate_size_factors <- function(count_data,
                                  col_data,
                                  design = "~ condition",
                                  min_row_sum = 10,
                                  reference_level = NULL) {
  message("Estimating size factors...")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = stats::as.formula(design)
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) > min_row_sum, ]

  if (!is.null(reference_level)) {
    dds$condition <- relevel(dds$condition, ref = reference_level)
  }

  control_size_factors <- DESeq2::sizeFactors(
    DESeq2::estimateSizeFactors(dds)
  )

  return(control_size_factors)
}


#' Compute size factors from control oligos
#'
#' A convenience wrapper around [estimate_size_factors()] with
#' sensible defaults for SGE control oligo normalization.
#'
#' @param count_matrix A matrix of control oligo counts.
#' @param sample_metadata A data frame with sample information. Must contain
#'   a `condition` column matching the sample columns.
#' @param reference_level Reference condition level (default: "Day4").
#'
#' @return A named numeric vector of size factors.
#'
#' @examples
#' \dontrun{
#' size_factors <- compute_control_size_factors(
#'   count_matrix = norm_matrix,
#'   sample_metadata = metadata,
#'   reference_level = "Day4"
#' )
#' }
#'
#' @export
compute_control_size_factors <- function(count_matrix,
                                         sample_metadata,
                                         reference_level = "Day4") {
  control_size_factors <- estimate_size_factors(
    count_data = count_matrix,
    col_data = sample_metadata,
    design = "~ condition",
    min_row_sum = 10,
    reference_level = reference_level
  )
  return(control_size_factors)
}
