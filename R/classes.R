#' S7 Class Definitions for sgeasy
#'
#' @description
#' S7 class definitions for SGE analysis data and results. These classes
#' provide formal containers with validation and clean print methods for
#' an improved REPL experience.
#'
#' @name sgeasy-classes
NULL


#' SGEData Class
#'
#' An S7 class that holds loaded SGE screen data before analysis.
#'
#' @param counts A data frame of raw count data with SEQUENCE and SAMPLE columns.
#' @param metadata A data frame of sample metadata with condition information.
#' @param annotation A data frame of variant annotations (optional, can be NULL).
#' @param complete_dataset A data frame of counts joined with metadata.
#'
#' @details
#' Access properties using the `@` operator:
#' \itemize{
#'   \item `data@counts` - Raw count data
#'   \item `data@metadata` - Sample metadata
#'   \item `data@annotation` - Variant annotations (may be NULL)
#'   \item `data@complete_dataset` - Counts joined with metadata
#' }
#'
#' @return An S7 object of class SGEData.
#'
#' @examples
#' \dontrun{
#' # Load data (returns SGEData object)
#' data <- load_sge_data(
#'   count_files = "inputs/",
#'   metadata_file = "metadata/screen_metadata.tsv"
#' )
#'
#' # Print shows summary
#' data
#'
#' # Access properties
#' data@counts
#' data@metadata
#' data@annotation
#' data@complete_dataset
#' }
#'
#' @export
SGEData <- S7::new_class(
 "SGEData",
  properties = list(
    counts = S7::class_data.frame,
    metadata = S7::class_data.frame,
    annotation = S7::new_property(
      class = S7::class_any,
      default = NULL,
      validator = function(value) {
        if (!is.null(value) && !is.data.frame(value)) {
          return("annotation must be a data frame or NULL")
        }
        NULL
      }
    ),
    complete_dataset = S7::class_data.frame
  ),
  validator = function(self) {
    if (!"SEQUENCE" %in% names(self@counts)) {
      return("counts must contain a SEQUENCE column")
    }
    if (!"SAMPLE" %in% names(self@counts)) {
      return("counts must contain a SAMPLE column")
    }
    if (!"condition" %in% names(self@metadata)) {
      return("metadata must contain a condition column")
    }
    if (!is.null(self@annotation) && !"Seq" %in% names(self@annotation)) {
      return("annotation must contain a Seq column")
    }
    NULL
  }
)

# Print method for SGEData - displays summary of loaded data
S7::method(print, SGEData) <- function(x, ...) {
  n_variants <- length(unique(x@counts$SEQUENCE))
  n_samples <- length(unique(x@counts$SAMPLE))
  n_conditions <- length(unique(x@metadata$condition))
  has_annotation <- !is.null(x@annotation)

  cli_h1("SGE Data")
  cli_ul(c(
    paste0("Variants: ", format(n_variants, big.mark = ",")),
    paste0("Samples: ", n_samples),
    paste0("Conditions: ", n_conditions),
    paste0("Annotation: ", if (has_annotation) "loaded" else "not loaded")
  ))
  cli_text("")
  cli_text("Access data with {.code @counts}, {.code @metadata}, {.code @annotation}, {.code @complete_dataset}")
  invisible(x)
}


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
