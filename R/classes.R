#' S7 Class Definitions for sgeasy
#'
#' @description
#' S7 class definitions for SGE analysis results. The `SGEResults` class
#' provides a formal container for differential analysis output with
#' validation and a clean print method.
#'
#' @name sgeasy-classes
NULL

#' SGEResults Class
#'
#' An S7 class that holds results from SGE differential abundance analysis.
#'
#' @param results A data frame (tibble) containing DESeq2 results with a
#'   SEQUENCE column identifying each variant.
#' @param rlog A DESeqTransform object containing regularized log-transformed
#'   counts for visualization.
#' @param contrast_summary A data frame (tibble) containing combined summary
#'   statistics across all contrasts, including log2 fold changes, p-values,
#'   z-scores, and optional rate estimates.
#' @param metadata A list of analysis metadata including condition levels,
#'   shrinkage type, alpha threshold, and analysis timestamp.
#'
#' @details
#' Access properties using the `@` operator:
#' \itemize{
#'   \item `results@results` - DESeq2 results table
#'   \item `results@rlog` - Regularized log-transformed data
#'   \item `results@contrast_summary` - Combined contrast statistics
#'   \item `results@metadata` - Analysis parameters and timestamp
#' }
#'
#' @return An S7 object of class SGEResults.
#'
#' @examples
#' \dontrun{
#' # Run analysis (returns SGEResults object)
#' results <- run_differential_analysis(counts, norm_counts, metadata)
#'
#' # Print shows summary
#' results
#'
#' # Access properties
#' results@results
#' results@contrast_summary
#' results@rlog
#'
#' # Get rlog matrix
#' SummarizedExperiment::assay(results@rlog)
#'
#' # Check analysis metadata
#' results@metadata$analysis_date
#' }
#'
#' @export
SGEResults <- S7::new_class(
  "SGEResults",
  properties = list(
    results = S7::class_data.frame,
    rlog = S7::new_property(
      class = S7::class_any,
      validator = function(value) {
        if (!inherits(value, "DESeqTransform")) {
          "rlog must be a DESeqTransform object"
        }
      }
    ),
    contrast_summary = S7::class_data.frame,
    metadata = S7::new_property(
      class = S7::class_list,
      default = list()
    )
  ),
  validator = function(self) {
    if (!"SEQUENCE" %in% names(self@results)) {
      return("results must contain a SEQUENCE column")
    }
    if (!"SEQUENCE" %in% names(self@contrast_summary)) {
      return("contrast_summary must contain a SEQUENCE column")
    }
    if (nrow(self@results) != nrow(self@contrast_summary)) {
      return("results and contrast_summary must have same number of rows")
    }
    NULL
  }
)

# Print method for SGEResults - displays summary of analysis results
S7::method(print, SGEResults) <- function(x, ...) {
  n_variants <- nrow(x@results)
  n_samples <- ncol(SummarizedExperiment::assay(x@rlog))
  contrasts <- names(x@contrast_summary) |>
    stringr::str_subset("^log2FoldChange_") |>
    stringr::str_remove("^log2FoldChange_")

  cli_h1("SGE Analysis Results")
  cli_ul(c(
    paste0("Variants: ", format(n_variants, big.mark = ",")),
    paste0("Samples: ", n_samples),
    paste0("Contrasts: ", length(contrasts))
  ))
  cli_text("")
  cli_text("Access data with {.code @results}, {.code @rlog}, {.code @contrast_summary}")
  invisible(x)
}
