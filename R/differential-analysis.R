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
#'   `condition`, `duration`, and `supplier_name` columns.
#' @param condition_levels Character vector specifying condition factor levels
#'   in order (default: c("Day4", "Day7", "Day10", "Day15")).
#' @param shrinkage_type Type of LFC shrinkage: "normal", "apeglm", or "ashr"
#'   (default: "normal").
#' @param alpha Significance threshold for adjusted p-values (default: 0.05).
#' @param include_rate Logical; whether to include continuous rate analysis
#'   (default: TRUE).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{results}{DESeq2 results table as a tibble with SEQUENCE column}
#'     \item{rlog}{Regularized log-transformed DESeqDataSet object}
#'     \item{contrast_summary}{Combined summary across all contrasts including
#'       log2 fold changes, p-values, z-scores, and rate estimates}
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
#'   \item Optionally runs continuous analysis for rate estimation
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
                                      alpha = 0.05,
                                      include_rate = TRUE) {
  message("Running DESeq2 differential abundance analysis on timepoint contrasts")

  # Build condition data for samples in the count matrix
  condition_data <- data.frame(
    condition = sample_metadata$condition[
      match(colnames(count_matrix), sample_metadata$supplier_name)
    ],
    duration = sample_metadata$duration[
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

  # Continuous analysis for rate estimation
 if (include_rate) {
    message("Running DESeq2 differential abundance analysis on ",
            "continuous data for rate estimation")

    # Estimate size factors for continuous design
    continuous_size_factors <- estimate_size_factors(
      count_data = normalization_matrix,
      col_data = condition_data,
      design = "~ duration",
      min_row_sum = 10,
      reference_level = NULL
    )

    continuous_dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = condition_data,
      design = ~ duration
    )
    continuous_dds <- DESeq2::estimateSizeFactors(continuous_dds)
    DESeq2::sizeFactors(continuous_dds) <- continuous_size_factors

    continuous_dds <- DESeq2::DESeq(continuous_dds)

    wald_continuous <- DEGreport::degComps(
      continuous_dds,
      combs = "duration",
      alpha = alpha,
      skip = FALSE,
      type = shrinkage_type,
      pairs = FALSE,
      fdr = "default"
    )

    rate <- wald_continuous[[1]] |>
      tibble::as_tibble(rownames = "SEQUENCE") |>
      dplyr::rename_with(~ paste0(.x, "_continuous"), .cols = -"SEQUENCE") |>
      dplyr::select(-dplyr::starts_with("baseMean_"))

    all_contrast_summary <- dplyr::left_join(
      all_contrast_summary,
      rate,
      by = "SEQUENCE"
    )
  }

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


#' Calculate median scores for log2 fold changes
#'
#' Computes the median of each log2FoldChange column from neutral variants
#' (synonymous and intronic) and adds new columns with the median values.
#'
#' @param data A data frame containing log2FoldChange columns and a
#'   Consequence column for filtering neutral variants.
#' @param neutral_variants Character vector of consequence types to use for
#'   median calculation (default: c("intron_variant", "synonymous_variant")).
#'
#' @return Data frame with additional median_log2FoldChange_* columns.
#'
#' @details
#' For proper statistical centering, the medians are calculated from
#' neutral variants (synonymous and intronic by default) that are expected
#' to have no functional effect. The function filters the data internally
#' to identify these variants, computes the median for each log2FoldChange
#' column, then broadcasts these values to all rows in the original data.
#'
#' @examples
#' \dontrun{
#' # Using default neutral variants (recommended)
#' data_with_medians <- calculate_median_scores(contrast_summary)
#'
#' # Using custom neutral variant types
#' data_with_medians <- calculate_median_scores(
#'   contrast_summary,
#'   neutral_variants = c("synonymous_variant")
#' )
#' }
#'
#' @export
calculate_median_scores <- function(data,
                                    neutral_variants = c("intron_variant",
                                                         "synonymous_variant")) {
  # Calculate medians from neutral variants using across()
  medians <- data |>
    dplyr::filter(.data$Consequence %in% neutral_variants) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("log2FoldChange"),
        ~ median(.x, na.rm = TRUE),
        .names = "median_{.col}"
      )
    ) |>
    dplyr::select(dplyr::starts_with("median_")) |>
    dplyr::slice_head(n = 1)

  # Bind median columns to original data
 dplyr::bind_cols(data, medians)
}


#' Add adjusted statistics for a single contrast
#'
#' Computes adjusted log fold change, z-score, p-value and FDR for a
#' specified contrast suffix.
#'
#' @param data A data frame with log2FoldChange, lfcSE, and median columns.
#' @param suffix The contrast suffix (e.g., "condition_Day7_vs_Day4").
#'
#' @return Data frame with additional adj_lfc_*, adj_score_*, pval_*, and
#'   FDR_* columns for the specified suffix.
#'
#' @details
#' The adjusted log fold change subtracts the median LFC (from neutral variants)
#' to center the distribution. The adjusted score is then computed as
#' adj_lfc / lfcSE, which follows approximately a standard normal distribution
#' under the null hypothesis. P-values are computed from a two-sided normal test.
#'
#' @keywords internal
add_adj_columns <- function(data, suffix) {
  lfc_col <- paste0("log2FoldChange_", suffix)
  se_col <- paste0("lfcSE_", suffix)
  median_col <- paste0("median_", lfc_col)

  data |>
    dplyr::mutate(
      "adj_lfc_{suffix}" := .data[[lfc_col]] - .data[[median_col]],
      "adj_score_{suffix}" := .data[[paste0("adj_lfc_", suffix)]] / .data[[se_col]],
      "pval_{suffix}" := stats::pnorm(abs(.data[[paste0("adj_score_", suffix)]]),
                                       lower.tail = FALSE) * 2,
      "FDR_{suffix}" := stats::p.adjust(.data[[paste0("pval_", suffix)]],
                                         method = "BH")
    )
}


#' Adjust statistics for all contrasts
#'
#' Applies [add_adj_columns()] to all contrasts found in the data.
#'
#' @param data A data frame with log2FoldChange columns for multiple contrasts.
#'
#' @return Data frame with adjusted statistics for all contrasts.
#'
#' @examples
#' \dontrun{
#' adjusted_data <- adjust_all_contrasts(contrast_summary)
#' }
#'
#' @export
adjust_all_contrasts <- function(data) {
  suffixes <- names(data) |>
    stringr::str_subset("^log2FoldChange_") |>
    stringr::str_remove("^log2FoldChange_")

  purrr::reduce(suffixes, \(d, s) add_adj_columns(d, s), .init = data)
}


#' Classify variants by consequence type
#'
#' Adds a simplified variant_class column based on the Consequence annotation.
#'
#' @param data A data frame with a Consequence column.
#'
#' @return Data frame with an additional variant_class column.
#'
#' @examples
#' \dontrun{
#' classified_data <- classify_variants(annotated_data)
#' }
#'
#' @export
classify_variants <- function(data) {
  data |>
    dplyr::mutate(
      variant_class = dplyr::case_when(
        .data$Consequence %in% "synonymous_variant" ~ "synonymous",
        .data$Consequence %in% "missense_variant" ~ "missense",
        .data$Consequence %in% "stop_gained" ~ "nonsense",
        .data$Consequence %in% "intron_variant" ~ "intronic",
        TRUE ~ "other"
      )
    )
}


#' Add functional classification for a single contrast
#'
#' Classifies variants as depleted, enriched, or unchanged based on
#' FDR threshold and log fold change direction.
#'
#' @param data A data frame with FDR and adj_lfc columns.
#' @param suffix The contrast suffix.
#' @param fdr_threshold FDR threshold for significance (default: 0.01).
#'
#' @return Data frame with an additional functional_classification_* column.
#'
#' @keywords internal
add_classification <- function(data, suffix, fdr_threshold = 0.01) {
  fdr_col <- paste0("FDR_", suffix)
  lfc_col <- paste0("adj_lfc_", suffix)

  data |>
    dplyr::mutate(
      "functional_classification_{suffix}" := dplyr::case_when(
        .data[[fdr_col]] < fdr_threshold & .data[[lfc_col]] < 0 ~ "depleted",
        .data[[fdr_col]] < fdr_threshold & .data[[lfc_col]] > 0 ~ "enriched",
        .data[[fdr_col]] >= fdr_threshold ~ "unchanged",
        TRUE ~ NA_character_
      )
    )
}


#' Classify variants for all contrasts
#'
#' Applies [add_classification()] to all contrasts found in the data.
#'
#' @param data A data frame with FDR columns for multiple contrasts.
#' @param fdr_threshold FDR threshold for significance (default: 0.01).
#'
#' @return Data frame with functional classification columns for all contrasts.
#'
#' @examples
#' \dontrun{
#' classified_data <- classify_all_contrasts(adjusted_data, fdr_threshold = 0.01)
#' }
#'
#' @export
classify_all_contrasts <- function(data, fdr_threshold = 0.01) {
  suffixes <- names(data) |>
    stringr::str_subset("^FDR_") |>
    stringr::str_remove("^FDR_")

  purrr::reduce(suffixes, \(d, s) add_classification(d, s, fdr_threshold),
                .init = data)
}


#' Recalculate screen statistics
#'
#' Applies the full post-processing pipeline: calculates medians from neutral
#' variants, adjusts contrasts, and classifies variants.
#'
#' @param data A data frame with contrast results from differential analysis.
#'   Must contain a Consequence column for identifying neutral variants.
#' @param fdr_threshold FDR threshold for functional classification
#'   (default: 0.01).
#'
#' @return Data frame with all computed statistics and classifications.
#'
#' @details
#' This function applies the standard SGE post-processing pipeline:
#' \enumerate{
#'   \item Calculates median log2FoldChange from neutral variants (synonymous
#'     and intronic) using [calculate_median_scores()]
#'   \item Adjusts all contrasts by subtracting medians and computing z-scores
#'     using [adjust_all_contrasts()]
#'   \item Classifies variants as depleted/enriched/unchanged based on FDR
#'     using [classify_all_contrasts()]
#' }
#'
#' @examples
#' \dontrun{
#' processed_data <- recalculate_screen_statistics(contrast_summary)
#'
#' # With custom FDR threshold
#' processed_data <- recalculate_screen_statistics(
#'   contrast_summary,
#'   fdr_threshold = 0.05
#' )
#' }
#'
#' @export
recalculate_screen_statistics <- function(data, fdr_threshold = 0.01) {
  data |>
    calculate_median_scores() |>
    adjust_all_contrasts() |>
    classify_all_contrasts(fdr_threshold = fdr_threshold)
}


#' Post-process screen results with annotation
#'
#' Combines differential analysis results with variant annotations, applies
#' statistical adjustments, and simplifies consequence types.
#'
#' @param data A data frame with contrast results from differential analysis.
#' @param targeton_id Character string identifying the targeton.
#' @param annotation A data frame with variant annotations. Must contain
#'   `Seq` and `Targeton_ID` columns for joining.
#' @param fdr_threshold FDR threshold for functional classification
#'   (default: 0.01).
#'
#' @return Data frame with annotations, adjusted statistics, functional
#'   classifications, and simplified consequence types.
#'
#' @details
#' This function provides a complete post-processing pipeline:
#' \enumerate{
#'   \item Adds the targeton identifier to the data
#'   \item Joins with the annotation file by SEQUENCE and Targeton_ID
#'   \item Applies [recalculate_screen_statistics()] for statistical adjustment
#'   \item Applies [slim_consequence()] for consequence simplification
#' }
#'
#' @examples
#' \dontrun{
#' processed <- post_process(
#'   data = contrast_summary,
#'   targeton_id = "GENE1_exon2",
#'   annotation = vep_annotations
#' )
#' }
#'
#' @seealso [recalculate_screen_statistics()], [slim_consequence()]
#' @export
post_process <- function(data,
                         targeton_id,
                         annotation,
                         fdr_threshold = 0.01) {
  data |>
    dplyr::mutate(Targeton_ID = targeton_id) |>
    dplyr::left_join(
      annotation,
      by = c("SEQUENCE" = "Seq", "Targeton_ID" = "Targeton_ID")
    ) |>
    recalculate_screen_statistics(fdr_threshold = fdr_threshold) |>
    slim_consequence()
}
