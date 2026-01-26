#' Stage Wrapper Functions for SGE Analysis
#'
#' These facade functions provide simplified interfaces for each major stage
#' of the SGE analysis workflow, designed for integration with UIs and
#' streamlined scripting.
#'
#' @name facades
#' @keywords internal
NULL


#' Load SGE data from files
#'
#' Loads count files, sample metadata, and variant annotations. This is the
#' first stage of the SGE analysis workflow.
#'
#' @param count_files Character vector of paths to count TSV files, or a
#'   directory path to search for count files.
#' @param metadata_file Path to the sample metadata TSV file.
#' @param annotation_file Path to the VEP annotation file. If NULL, will
#'   attempt to extract from metadata `vep_anno` column.
#' @param count_file_pattern Glob pattern for finding count files if
#'   `count_files` is a directory (default: "*lib_counts.tsv.gz").
#' @param exclude_patterns Character vector of patterns to exclude from
#'   count file discovery (default: NULL).
#' @param metadata_join_by Column name in metadata to join with SAMPLE
#'   column in counts (default: "sequencing_id").
#' @param calculate_duration Logical; whether to calculate duration from
#'   condition column (default: TRUE).
#'
#' @return A list containing:
#'   \describe{
#'     \item{counts}{The raw count data frame}
#'     \item{metadata}{The sample metadata data frame}
#'     \item{annotation}{The variant annotation data frame}
#'     \item{complete_dataset}{Counts joined with metadata}
#'   }
#'
#' @examples
#' \dontrun{
#' data <- load_sge_data(
#'   count_files = "inputs/",
#'   metadata_file = "metadata/screen_metadata.tsv",
#'   annotation_file = "annotations/vep_annotations.tsv"
#' )
#' }
#'
#' @export
load_sge_data <- function(count_files,
                          metadata_file,
                          annotation_file = NULL,
                          count_file_pattern = "*lib_counts.tsv.gz",
                          exclude_patterns = NULL,
                          metadata_join_by = "sequencing_id",
                          calculate_duration = TRUE) {
  message("Loading SGE data...")

  # Find count files if directory provided
  if (length(count_files) == 1 && dir.exists(count_files)) {
    count_files <- list.files(
      count_files,
      pattern = utils::glob2rx(count_file_pattern),
      recursive = TRUE,
      full.names = TRUE
    )
  }

  # Apply exclusion patterns
  if (!is.null(exclude_patterns)) {
    for (pattern in exclude_patterns) {
      count_files <- count_files[!grepl(pattern, count_files)]
    }
  }

  message(sprintf("  Found %d count files", length(count_files)))

  # Read counts
  counts <- readr::read_tsv(count_files, comment = "##", show_col_types = FALSE)
  message(sprintf("  Loaded %d rows of count data", nrow(counts)))

  # Read metadata
  metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE)
  message(sprintf("  Loaded %d rows of metadata", nrow(metadata)))

  # Calculate duration if requested
  if (calculate_duration && "condition" %in% names(metadata)) {
    metadata <- metadata |>
      dplyr::mutate(
        duration = as.numeric(
          stringr::str_remove(.data$condition, pattern = "Day")
        ) - 4
      )
  }

  # Get annotation file from metadata if not provided
  if (is.null(annotation_file) && "vep_anno" %in% names(metadata)) {
    annotation_file <- unique(metadata$vep_anno)[1]
    message(sprintf("  Using annotation file from metadata: %s", annotation_file))
  }

  # Read annotation
  annotation <- NULL
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    annotation <- readr::read_tsv(annotation_file, show_col_types = FALSE)
    message(sprintf("  Loaded %d rows of annotation data", nrow(annotation)))
  }

  # Join counts with metadata
  join_spec <- stats::setNames("SAMPLE", metadata_join_by)
  complete_dataset <- dplyr::left_join(counts, metadata, by = join_spec)
  message(sprintf("  Created complete dataset with %d rows", nrow(complete_dataset)))

  list(
    counts = counts,
    metadata = metadata,
    annotation = annotation,
    complete_dataset = complete_dataset
  )
}


#' Prepare count matrices for analysis
#'
#' Transforms the loaded data into analysis-ready count matrices. This includes
#' filtering, pivoting, deduplication, annotation joining, and artefact removal.
#'
#' @param complete_dataset Data frame from `load_sge_data()$complete_dataset`.
#' @param annotation Data frame of variant annotations (from `load_sge_data()`).
#' @param exclude_samples Character vector of sample names to exclude
#'   (default: NULL).
#' @param exclude_targetons Character vector of targeton IDs to exclude
#'   (default: NULL).
#' @param min_counts Minimum total count threshold (default: 10).
#' @param sample_col Column name for sample identifier (default: "supplier_name").
#' @param count_col Column name for count values (default: "COUNT").
#' @param targeton_col Column name for targeton ID (default: "targeton_id").
#'
#' @return A list containing:
#'   \describe{
#'     \item{annotated_counts}{List of annotated count data frames by targeton}
#'     \item{normalisation_matrices}{List of normalization matrices by targeton}
#'     \item{complete_matrices}{List of full count matrices by targeton}
#'     \item{targeton_ids}{Character vector of targeton IDs}
#'   }
#'
#' @examples
#' \dontrun{
#' matrices <- prepare_count_matrices(
#'   complete_dataset = data$complete_dataset,
#'   annotation = data$annotation,
#'   exclude_samples = c("bad_sample_1"),
#'   min_counts = 10
#' )
#' }
#'
#' @export
prepare_count_matrices <- function(complete_dataset,
                                   annotation,
                                   exclude_samples = NULL,
                                   exclude_targetons = NULL,
                                   min_counts = 10,
                                   sample_col = "supplier_name",
                                   count_col = "COUNT",
                                   targeton_col = "targeton_id") {
  message("Preparing count matrices...")

  # Filter samples
  if (!is.null(exclude_samples)) {
    complete_dataset <- complete_dataset |>
      dplyr::filter(!.data[[sample_col]] %in% exclude_samples) |>
      dplyr::filter(!.data$SAMPLE %in% exclude_samples)
    message(sprintf("  Excluded %d sample(s)", length(exclude_samples)))
  }

  # Filter targetons
  if (!is.null(exclude_targetons)) {
    complete_dataset <- complete_dataset |>
      dplyr::filter(!.data[[targeton_col]] %in% exclude_targetons)
    message(sprintf("  Excluded %d targeton(s)", length(exclude_targetons)))
  }

  # Split by targeton
  message("  Splitting by targeton...")
  targeton_matrices <- complete_dataset |>
    dplyr::group_by(.data[[targeton_col]]) |>
    dplyr::group_split()

  targeton_ids <- purrr::map_chr(
    targeton_matrices,
    ~ unique(.x[[targeton_col]])
  )
  targeton_matrices <- purrr::set_names(targeton_matrices, targeton_ids)
  message(sprintf("  Found %d targeton(s)", length(targeton_ids)))

  # Pivot to wide format
  message("  Pivoting to wide format...")
  count_matrices <- purrr::map(
    targeton_matrices,
    ~ tidyr::pivot_wider(
      .x,
      names_from = dplyr::all_of(sample_col),
      values_from = dplyr::all_of(count_col),
      names_prefix = "count_",
      id_cols = c("NAME", "SEQUENCE", "LENGTH", targeton_col)
    )
  )

  # Deduplicate
  message("  Deduplicating sequences...")
  unique_matrices <- purrr::map(
    count_matrices,
    ~ dplyr::group_by(.x, .data$SEQUENCE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
  )

  # Filter by counts
  message(sprintf("  Filtering by minimum counts >= %d...", min_counts))
  filtered_matrices <- purrr::map(
    unique_matrices,
    filter_by_counts,
    min_counts = min_counts
  )

  # Join annotations
  if (!is.null(annotation)) {
    message("  Joining annotations...")
    annotated_counts <- purrr::map(
      filtered_matrices,
      ~ dplyr::left_join(
        .x,
        annotation,
        by = c("SEQUENCE" = "Seq", stats::setNames("Targeton_ID", targeton_col))
      )
    )

    # Remove artefacts
    message("  Removing artefacts...")
    annotated_counts <- purrr::map(annotated_counts, remove_artefacts)
  } else {
    annotated_counts <- filtered_matrices
  }

  # Create matrices
  message("  Creating normalization matrices...")
  normalisation_matrices <- purrr::map(
    annotated_counts,
    create_normalization_matrix
  )

  message("  Creating complete matrices...")
  complete_matrices <- purrr::map(
    annotated_counts,
    create_count_matrix
  )

  message("  Done preparing count matrices")

  list(
    annotated_counts = annotated_counts,
    normalisation_matrices = normalisation_matrices,
    complete_matrices = complete_matrices,
    targeton_ids = targeton_ids
  )
}


#' Run differential analysis on all screens
#'
#' Runs DESeq2-based differential analysis on all targetons in the prepared
#' count matrices.
#'
#' @param complete_matrices List of count matrices from `prepare_count_matrices()`.
#' @param normalisation_matrices List of normalization matrices from
#'   `prepare_count_matrices()`.
#' @param metadata Sample metadata data frame.
#' @param condition_levels Character vector of condition factor levels
#'   (default: c("Day4", "Day7", "Day10", "Day15")).
#' @param ... Additional arguments passed to `run_differential_analysis()`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{deseq_results}{List of [SGEResults] objects by targeton}
#'     \item{contrast_tables}{Combined contrast table for all targetons}
#'     \item{rlog_results}{List of rlog-transformed data by targeton}
#'   }
#'
#' @examples
#' \dontrun{
#' results <- analyze_screens(
#'   complete_matrices = matrices$complete_matrices,
#'   normalisation_matrices = matrices$normalisation_matrices,
#'   metadata = data$metadata
#' )
#' }
#'
#' @export
analyze_screens <- function(complete_matrices,
                            normalisation_matrices,
                            metadata,
                            condition_levels = c("Day4", "Day7", "Day10", "Day15"),
                            ...) {
  message("Running differential analysis...")

  deseq_results <- purrr::map2(
    complete_matrices,
    normalisation_matrices,
    run_differential_analysis,
    sample_metadata = metadata,
    condition_levels = condition_levels,
    ...
  )

  message("  Combining contrast tables...")
  contrast_tables <- purrr::map_dfr(
    deseq_results,
    \(x) x@contrast_summary,
    .id = "Targeton_ID"
  )

  rlog_results <- purrr::map(deseq_results, \(x) x@rlog)

  message(sprintf("  Analysis complete for %d targeton(s)", length(deseq_results)))

  list(
    deseq_results = deseq_results,
    contrast_tables = contrast_tables,
    rlog_results = rlog_results
  )
}


#' Generate QC plots
#'
#' Creates quality control visualizations including sample distance matrices
#' and PCA plots.
#'
#' @param rlog_results List of rlog-transformed data by targeton (from
#'   `analyze_screens()$rlog_results`).
#' @param output_dir Directory to save plots (default: "results/qc").
#' @param create_distance_matrices Logical; create distance matrix plots
#'   (default: TRUE).
#' @param create_pca_plots Logical; create PCA plots (default: TRUE).
#'
#' @return Invisibly returns a list of created plot objects.
#'
#' @examples
#' \dontrun{
#' generate_qc_plots(
#'   rlog_results = results$rlog_results,
#'   output_dir = "results/qc"
#' )
#' }
#'
#' @export
generate_qc_plots <- function(rlog_results,
                              output_dir = "results/qc",
                              create_distance_matrices = TRUE,
                              create_pca_plots = TRUE) {
  message("Generating QC plots...")

  plots <- list()

  if (create_distance_matrices) {
    message("  Creating distance matrices...")
    plots$distance <- purrr::map(rlog_results, sample_distance_matrix)
    plot_distance_matrices(rlog_results, output_dir = output_dir)
  }

  if (create_pca_plots) {
    message("  Creating PCA plots...")
    plots$pca <- purrr::map(rlog_results, build_sample_pca)
    plot_screen_pcas(rlog_results, output_dir = output_dir)
  }

  message(sprintf("  QC plots saved to: %s", output_dir))
  invisible(plots)
}


#' Generate analysis results
#'
#' Post-processes contrast tables, creates result visualizations, and
#' optionally performs replicate reweighting.
#'
#' @param contrast_tables Combined contrast table from `analyze_screens()`.
#' @param annotation Variant annotation data frame.
#' @param output_dir Directory to save results (default: "results").
#' @param create_waterfall_plots Logical; create waterfall plots (default: TRUE).
#' @param reweight_replicates Logical; perform replicate reweighting
#'   (default: TRUE).
#' @param fdr_threshold FDR threshold for classification (default: 0.01).
#'
#' @return A list containing:
#'   \describe{
#'     \item{processed_data}{Post-processed contrast table}
#'     \item{reweighted_data}{Reweighted variant data (if reweight_replicates = TRUE)}
#'     \item{plots}{List of created plot objects}
#'   }
#'
#' @examples
#' \dontrun{
#' final_results <- generate_results(
#'   contrast_tables = results$contrast_tables,
#'   annotation = data$annotation,
#'   output_dir = "results"
#' )
#' }
#'
#' @export
generate_results <- function(contrast_tables,
                             annotation,
                             output_dir = "results",
                             create_waterfall_plots = TRUE,
                             reweight_replicates = TRUE,
                             fdr_threshold = 0.01) {
  message("Generating results...")

  # Post-process
  message("  Post-processing contrast tables...")
  processed_data <- post_process(
    data = contrast_tables,
    annotation = annotation,
    fdr_threshold = fdr_threshold
  )

  plots <- list()

  # Create waterfall plots
  if (create_waterfall_plots) {
    message("  Creating waterfall plots...")
    individual_screens <- processed_data |>
      dplyr::group_by(.data$Targeton_ID) |>
      dplyr::group_split() |>
      purrr::set_names(unique(processed_data$Targeton_ID))

    plots$waterfall <- purrr::imap(
      individual_screens,
      ~ plot_all_conditions(.x, .y)
    )

    write_waterfall_plots(
      plots$waterfall,
      output_dir = file.path(output_dir, "screen")
    )
  }

  # Reweight replicates
  reweighted_data <- NULL
  if (reweight_replicates &&
      all(c("HGVSc", "HGVSp", "pam_mut_sgrna_id") %in% names(processed_data))) {
    message("  Reweighting replicated variants...")
    reweight_input <- processed_data |>
      dplyr::select(
        "Targeton_ID", "SEQUENCE", "HGVSc", "HGVSp",
        "pam_mut_sgrna_id", dplyr::contains("continuous")
      )
    reweighted_data <- reweight_replicated_variants(
      reweight_input,
      fdr_threshold = fdr_threshold
    )
    message(sprintf("    Reweighted %d unique variants", nrow(reweighted_data)))
  }

  # Save results
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  readr::write_tsv(
    processed_data,
    file.path(output_dir, "processed_contrast_tables.tsv")
  )
  message(sprintf("  Results saved to: %s", output_dir))

  if (!is.null(reweighted_data)) {
    readr::write_tsv(
      reweighted_data,
      file.path(output_dir, "reweighted_variants.tsv")
    )
  }

  list(
    processed_data = processed_data,
    reweighted_data = reweighted_data,
    plots = plots
  )
}
