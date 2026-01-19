#' Create normalization matrix from neutral variants
#'
#' Extracts synonymous and intronic variants from annotated data to create
#' a matrix suitable for size factor estimation.
#'
#' @param data A data frame with variant annotations including a
#'   `Consequence` column and sample count columns.
#' @param consequence_types Character vector of consequence types to include
#'   (default: c("synonymous_variant", "intron_variant")).
#' @param sample_pattern Regex pattern to identify sample columns
#'   (default: "pel").
#'
#' @return A numeric matrix with variants as rows and samples as columns.
#'   Row names are set to the SEQUENCE column values.
#'
#' @details
#' This function filters the input data to retain only neutral variants
#' (synonymous and intronic by default) which serve as controls for
#' normalization. These variants are expected to have no functional effect
#' and thus provide a baseline for estimating library size factors.
#'
#' @examples
#' \dontrun{
#' norm_matrix <- create_normalization_matrix(
#'   data = annotated_counts,
#'   consequence_types = c("synonymous_variant", "intron_variant")
#' )
#' }
#'
#' @export
create_normalization_matrix <- function(data,
                                        consequence_types = c("synonymous_variant",
                                                              "intron_variant"),
                                        sample_pattern = "pel") {
  data |>
    dplyr::filter(.data$Consequence %in% consequence_types) |>
    dplyr::select("SEQUENCE", dplyr::contains(sample_pattern)) |>
    tibble::column_to_rownames("SEQUENCE") |>
    as.matrix()
}


#' Create count matrix from annotated data
#'
#' Extracts sequence identifiers and sample counts into a numeric matrix
#' suitable for differential analysis.
#'
#' @param data A data frame with SEQUENCE column and sample count columns.
#' @param sample_pattern Regex pattern to identify sample columns
#'   (default: "pel").
#'
#' @return A numeric matrix with sequences as row names and samples as columns.
#'
#' @examples
#' \dontrun{
#' count_matrix <- create_count_matrix(annotated_data)
#' }
#'
#' @export
create_count_matrix <- function(data, sample_pattern = "pel") {
  data |>
    dplyr::select("SEQUENCE", dplyr::contains(sample_pattern)) |>
    tibble::column_to_rownames("SEQUENCE") |>
    as.matrix()
}


#' Filter variants by count thresholds
#'
#' Removes variants with total counts outside specified bounds and
#' optionally excludes specific sample patterns.
#'
#' @param data A data frame with count data.
#' @param min_counts Minimum total count to retain (default: 10).
#' @param max_counts Maximum total count to retain (default: Inf).
#' @param exclude_pattern Pattern for columns to exclude from output
#'   (default: "_0_"). Set to NULL to keep all columns.
#' @param sample_pattern Regex pattern to identify sample columns for
#'   count summation (default: "pel").
#'
#' @return Filtered data frame with variants meeting count criteria.
#'
#' @details
#' This function calculates the total count for each variant across all
#' sample columns matching `sample_pattern`, then filters to retain only
#' variants within the specified count range. Optionally removes columns
#' matching `exclude_pattern` (e.g., day 0 samples).
#'
#' @examples
#' \dontrun{
#' filtered_data <- filter_by_counts(
#'   data = count_data,
#'   min_counts = 10,
#'   exclude_pattern = "_0_"
#' )
#' }
#'
#' @export
filter_by_counts <- function(data,
                             min_counts = 10,
                             max_counts = Inf,
                             exclude_pattern = "_0_",
                             sample_pattern = "pel") {
  result <- data |>
    dplyr::rowwise() |>
    dplyr::mutate(
      total_counts = sum(dplyr::across(dplyr::contains(sample_pattern)),
                         na.rm = TRUE)
    ) |>
    dplyr::filter(dplyr::between(.data$total_counts, min_counts, max_counts)) |>
    dplyr::ungroup() |>
    dplyr::select(-"total_counts")

  if (!is.null(exclude_pattern)) {
    result <- result |>
      dplyr::select(-dplyr::contains(exclude_pattern))
  }

  result
}


#' Simplify variant consequence annotations
#'
#' Creates a simplified consequence annotation column by mapping VEP consequences
#' to a reduced set of categories, considering mutator type where relevant.
#'
#' @param data A data frame with Consequence and mutator columns.
#'
#' @return Data frame with an additional slim_consequence column.
#'
#' @details
#' The function maps VEP consequences to simplified categories:
#' \itemize{
#'   \item Clinical inframe deletions/insertions (custom mutator)
#'   \item Synonymous (including stop_retained with appropriate mutators)
#'   \item Start/stop lost variants
#'   \item Codon deletions (inframe mutator)
#'   \item Stop gained, frameshift, splice variants
#'   \item Missense, intron, UTR variants
#' }
#'
#' @examples
#' \dontrun{
#' data_with_slim <- slim_consequence(annotated_data)
#' }
#'
#' @export
slim_consequence <- function(data) {
  data |>
    dplyr::mutate(
      slim_consequence = dplyr::case_when(
        # More specific conditions first (with mutator requirements)
        stringr::str_detect(.data$Consequence, "inframe_deletion") &
          .data$mutator == "custom" ~ "clinical_inframe_deletion",
        stringr::str_detect(.data$Consequence, "inframe_insertion") &
          .data$mutator == "custom" ~ "clinical_inframe_insertion",
        stringr::str_detect(.data$Consequence, "stop_retained_variant") &
          .data$mutator %in% c("snvre", "snv", "custom") ~ "synonymous",
        stringr::str_detect(.data$Consequence, "start_lost") &
          .data$mutator %in% c("snvre", "snv", "custom", "ala") ~ "start_lost",
        stringr::str_detect(.data$Consequence, "stop_lost") &
          .data$mutator %in% c("snvre", "snv", "custom", "ala") ~ "stop_lost",
        stringr::str_detect(.data$mutator, "inframe") ~ "codon_deletion",

        # Simple pattern matches
        stringr::str_detect(.data$Consequence, "stop_gained") ~ "stop_gained",
        stringr::str_detect(.data$Consequence, "frameshift_variant") ~ "frameshift",
        stringr::str_detect(.data$Consequence, "splice_acceptor_variant") ~ "splice_acceptor",
        stringr::str_detect(.data$Consequence, "splice_donor_variant") ~ "splice_donor",
        stringr::str_detect(.data$Consequence, "missense_variant") ~ "missense",
        stringr::str_detect(.data$Consequence, "synonymous_variant") ~ "synonymous",
        stringr::str_detect(.data$Consequence, "intron_variant") ~ "intron",
        stringr::str_detect(.data$Consequence, "UTR_variant") ~ "UTR",

        TRUE ~ NA_character_
      )
    )
}
