#' Create normalization matrix from neutral variants
#'
#' Extracts synonymous and intronic variants from annotated data to create
#' a matrix suitable for size factor estimation.
#'
#' @param data A data frame with variant annotations including a
#'   `Consequence` column and sample count columns.
#' @param consequence_types Character vector of consequence types to include
#'   (default: c("synonymous_variant", "intron_variant")).
#' @param sample_pattern Pattern to identify sample columns
#'   (default: "count").
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
#' The sample_pattern is used with `starts_with()` to select count columns.
#' When using pivot_wider with `names_prefix = "count_"`, columns will be
#' named like "count_SAMPLE_A_10_pel0000", and "count" will match these.
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
                                        sample_pattern = "count") {
  create_count_matrix(data, sample_pattern, consequence_types)
}


#' Create count matrix from annotated data
#'
#' Extracts sequence identifiers and sample counts into a numeric matrix
#' suitable for differential analysis. Optionally filters by consequence type.
#'
#' @param data A data frame with SEQUENCE column and sample count columns.
#' @param sample_pattern Pattern to identify sample columns
#'   (default: "count").
#' @param consequence_types Optional character vector of consequence types to
#'   filter on (requires a `Consequence` column in data). If NULL (default),
#'   no filtering is applied.
#'
#' @return A numeric matrix with sequences as row names and samples as columns.
#'
#' @details
#' The sample_pattern is used with `starts_with()` to select count columns.
#' When using pivot_wider with `names_prefix = "count_"`, columns will be
#' named like "count_SAMPLE_A_10_pel0000", and "count" will match these.
#'
#' @examples
#' \dontrun{
#' count_matrix <- create_count_matrix(annotated_data)
#'
#' # Filter to specific consequence types
#' norm_matrix <- create_count_matrix(
#'   annotated_data,
#'   consequence_types = c("synonymous_variant", "intron_variant")
#' )
#' }
#'
#' @export
create_count_matrix <- function(data,
                                sample_pattern = "count",
                                consequence_types = NULL) {
  result <- data
  if (!is.null(consequence_types)) {
    result <- dplyr::filter(result, .data$Consequence %in% consequence_types)
  }
  result |>
    dplyr::select("SEQUENCE", dplyr::starts_with(sample_pattern)) |>
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
#' @param sample_pattern Pattern to identify sample columns for
#'   count summation (default: "count").
#'
#' @return Filtered data frame with variants meeting count criteria.
#'
#' @details
#' This function calculates the total count for each variant across all
#' sample columns matching `sample_pattern` using `starts_with()`, then
#' filters to retain only variants within the specified count range.
#' Optionally removes columns matching `exclude_pattern` (e.g., day 0 samples).
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
                             sample_pattern = "count") {
  result <- data |>
    dplyr::rowwise() |>
    dplyr::mutate(
      total_counts = sum(dplyr::across(dplyr::starts_with(sample_pattern)),
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
          .data$mutator %in% c("snvre", "snv", "custom", "ala", "aa") ~ "start_lost",
        stringr::str_detect(.data$Consequence, "stop_lost") &
          .data$mutator %in% c("snvre", "snv", "custom", "ala", "aa") ~ "stop_lost",
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


#' Calculate position effect ratios
#'
#' Calculates the ratio of mean counts between a specified timepoint and a
#' reference condition for each sequence. Counts are averaged across
#' replicates within each condition before computing the ratio.
#'
#' @param dataframe A data frame with COUNT, condition, NAME, SEQUENCE, and
#'   targeton_id columns.
#' @param annotation A data frame with sequence annotations containing Seq
#'   and Targeton_ID columns.
#' @param timepoint Character string specifying the numerator timepoint
#'   (e.g., "Day4", "Day10").
#' @param reference Character string specifying the denominator condition
#'   (default: "Day0").
#'
#' @return A data frame with mean counts per condition, ratio calculations,
#'   and joined annotation data including Targeton_ID and vcf_pos.
#'
#' @examples
#' \dontrun{
#' pos_effect <- calculate_position_effect(
#'   dataframe = count_data,
#'   annotation = annotation_df,
#'   timepoint = "Day4"
#' )
#'
#' # Using Day4 as reference instead of plasmid
#' pos_effect <- calculate_position_effect(
#'   dataframe = count_data,
#'   annotation = annotation_df,
#'   timepoint = "Day10",
#'   reference = "Day4"
#' )
#' }
#'
#' @export
calculate_position_effect <- function(dataframe, annotation, timepoint,
                                      reference = "Day0") {
  pos_ratio <- dataframe |>
    dplyr::filter(.data$condition %in% c(reference, timepoint)) |>
    dplyr::group_by(.data$targeton_id, .data$NAME, .data$condition,
                    .data$SEQUENCE) |>
    dplyr::summarise(mean_counts = mean(.data$COUNT), .groups = "drop") |>
    tidyr::pivot_wider(names_from = "condition", values_from = "mean_counts") |>
    dplyr::mutate(ratio = .data[[timepoint]] / .data[[reference]]) |>
    dplyr::ungroup()

  pos_ratio <- pos_ratio |>
    dplyr::left_join(
      annotation,
      by = c("SEQUENCE" = "Seq", "targeton_id" = "Targeton_ID")
    )

  pos_ratio
}


#' Extend a loess fit to cover edge positions
#'
#' Fills NA values at the start and end of a fitted vector with the nearest
#' non-NA value. Interior NAs are left unchanged. This handles cases where
#' the loess prediction cannot extrapolate beyond the range of the training
#' data.
#'
#' @param fit_vector A numeric vector of fitted values, possibly with NAs at
#'   the edges.
#'
#' @return A numeric vector of the same length with edge NAs filled.
#'
#' @keywords internal
extend_loess_fit <- function(fit_vector) {
  if (!any(is.na(fit_vector))) return(fit_vector)

  valid_idx <- which(!is.na(fit_vector))
  if (length(valid_idx) == 0) return(fit_vector)

  first_valid <- min(valid_idx)
  last_valid <- max(valid_idx)

  out <- fit_vector
  if (first_valid > 1) {
    out[1:(first_valid - 1)] <- fit_vector[first_valid]
  }
  if (last_valid < length(fit_vector)) {
    out[(last_valid + 1):length(fit_vector)] <- fit_vector[last_valid]
  }
  out
}


#' Correct for positional effects using neutral variants
#'
#' Fits a loess regression on neutral (e.g. synonymous) variants' log2 ratios
#' as a function of genomic position, then subtracts the fitted positional
#' bias from all variants. This is performed independently per targeton.
#'
#' This is analogous to the median-based LFC adjustment performed by
#' [calculate_median_scores()], but operates in the spatial dimension
#' rather than as a single global correction.
#'
#' @param pos_ratio A data frame from [calculate_position_effect()] containing
#'   at minimum: \code{ratio}, \code{vcf_pos}, \code{Consequence}, and
#'   \code{targeton_id} columns.
#' @param neutral_variants Character vector of VEP consequence types to use
#'   as the neutral reference (default: \code{"synonymous_variant"}).
#' @param span Loess span parameter controlling smoothness of the fit.
#'   Lower values produce more local fits that capture cut-site effects;
#'   higher values produce smoother global trends (default: 0.3).
#' @param min_neutral Minimum number of neutral variants with valid ratios
#'   required to fit the model. Targetons with fewer neutral variants will
#'   not be corrected and a warning is issued (default: 20).
#'
#' @return The input data frame with additional columns:
#'   \describe{
#'     \item{pos_effect}{The fitted positional bias (log2 scale) at each
#'       variant's position}
#'     \item{corrected_log2_ratio}{log2(ratio) minus the positional effect}
#'   }
#'
#' @examples
#' \dontrun{
#' pos_data <- calculate_position_effect(count_data, annotation, "Day4")
#' corrected <- correct_position_effect(pos_data)
#'
#' # With custom parameters
#' corrected <- correct_position_effect(
#'   pos_data,
#'   neutral_variants = c("synonymous_variant", "intron_variant"),
#'   span = 0.15
#' )
#' }
#'
#' @seealso [calculate_position_effect()], [calculate_median_scores()]
#' @export
correct_position_effect <- function(pos_ratio,
                                    neutral_variants = c("synonymous_variant"),
                                    span = 0.3,
                                    min_neutral = 20) {
  if (!"targeton_id" %in% names(pos_ratio)) {
    stop("pos_ratio must contain a targeton_id column")
  }

  targetons <- split(pos_ratio, pos_ratio$targeton_id)

  purrr::map_dfr(targetons, function(tgt_data) {
    .fit_targeton_position_effect(
      tgt_data, neutral_variants, span, min_neutral
    )
  })
}


#' Fit positional effect model for a single targeton
#'
#' @param data Data frame for a single targeton with ratio, vcf_pos, and
#'   Consequence columns.
#' @param neutral_variants Character vector of neutral consequence types.
#' @param span Loess span parameter.
#' @param min_neutral Minimum number of neutral variants required.
#'
#' @return Data frame with \code{pos_effect} and \code{corrected_log2_ratio}
#'   columns added.
#'
#' @keywords internal
.fit_targeton_position_effect <- function(data, neutral_variants,
                                          span, min_neutral) {
  targeton_id <- unique(data$targeton_id)


  # Filter to neutral variants with valid ratios, excluding PAM mutations
  neutral_data <- data |>
    dplyr::filter(
      .data$Consequence %in% neutral_variants,
      is.na(.data$pam_mut_sgrna_id) | .data$pam_mut_sgrna_id == "",
      is.finite(.data$ratio),
      .data$ratio > 0
    )

  if (nrow(neutral_data) < min_neutral) {
    warning(
      "Targeton '", targeton_id, "' has only ", nrow(neutral_data),
      " neutral variants (minimum: ", min_neutral, "). ",
      "Skipping positional correction.",
      call. = FALSE
    )
    data$pos_effect <- NA_real_
    data$corrected_log2_ratio <- dplyr::if_else(
      is.finite(data$ratio) & data$ratio > 0,
      log2(data$ratio),
      NA_real_
    )
    return(data)
  }

  # Build position sequence spanning all variants
  all_positions <- data$vcf_pos[!is.na(data$vcf_pos)]
  pos_range <- range(all_positions)
  pos_seq <- seq(from = pos_range[1], to = pos_range[2])

  # Fit loess on neutral variants only
  lo <- stats::loess(
    log2(ratio) ~ vcf_pos,
    data = neutral_data,
    span = span
  )

  # Predict across full position range and extend edges
  pred <- stats::predict(lo, newdata = data.frame(vcf_pos = pos_seq))
  pred <- extend_loess_fit(pred)

  # Create lookup table
  pos_effect_lookup <- data.frame(
    vcf_pos = pos_seq,
    pos_effect = pred
  )

  # Join positional effect and compute corrected ratio
  data |>
    dplyr::left_join(pos_effect_lookup, by = "vcf_pos") |>
    dplyr::mutate(
      corrected_log2_ratio = dplyr::if_else(
        is.finite(.data$ratio) & .data$ratio > 0,
        log2(.data$ratio) - .data$pos_effect,
        NA_real_
      )
    )
}


#' Calculate position effects for multiple timepoints
#'
#' Convenience wrapper that calls [calculate_position_effect()] for each
#' timepoint and returns a named list of data frames.
#'
#' @inheritParams calculate_position_effect
#' @param timepoints Character vector of timepoint conditions to compute
#'   position effects for (e.g., \code{c("Day4", "Day7", "Day10", "Day15")}).
#'
#' @return A named list of data frames, one per timepoint, each containing
#'   the output of [calculate_position_effect()].
#'
#' @examples
#' \dontrun{
#' all_pos <- calculate_all_position_effects(
#'   dataframe = count_data,
#'   annotation = annotation_df,
#'   timepoints = c("Day4", "Day7", "Day10", "Day15")
#' )
#' }
#'
#' @seealso [calculate_position_effect()], [correct_all_position_effects()]
#' @export
calculate_all_position_effects <- function(dataframe, annotation,
                                           timepoints, reference = "Day0") {
  purrr::map(timepoints, ~calculate_position_effect(dataframe, annotation, .x, reference)) |>
    purrr::set_names(timepoints)
}


#' Correct position effects for multiple timepoints
#'
#' Applies [slim_consequence()], [remove_artefacts()], and
#' [correct_position_effect()] to each element of a position effect list
#' (as returned by [calculate_all_position_effects()]).
#'
#' @param pos_effect_list A named list of data frames, as returned by
#'   [calculate_all_position_effects()].
#' @param span Loess span parameter passed to [correct_position_effect()]
#'   (default: 0.3).
#' @param ... Additional arguments passed to [correct_position_effect()].
#'
#' @return A named list of data frames, each with positional correction applied.
#'
#' @examples
#' \dontrun{
#' all_pos <- calculate_all_position_effects(
#'   dataframe = count_data,
#'   annotation = annotation_df,
#'   timepoints = c("Day4", "Day7", "Day10", "Day15")
#' )
#' corrected_list <- correct_all_position_effects(all_pos)
#' }
#'
#' @seealso [calculate_all_position_effects()], [correct_position_effect()]
#' @export
correct_all_position_effects <- function(pos_effect_list, span = 0.3, ...) {
  purrr::map(pos_effect_list, ~{
    slim_consequence(.x) |> remove_artefacts() |> correct_position_effect(span = span, ...)
  })
}
