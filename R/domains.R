# Protein Domain Visualization Functions
# =======================================
# Functions for fetching protein domains from Ensembl and creating
# combined domain track + amino acid heatmap visualizations.

# Exported Data Objects --------------------------------------------------------

#' Classification colours for amino acid heatmap
#'
#' A named character vector of colours for functional classification categories.
#'
#' @export
classification_heatmap_colours <- c(
  "unchanged" = "#CCCCCC",
  "enriched"  = "#FF5733",
  "depleted"  = "#3C5488FF"
)

#' Standard amino acid ordering for heatmap y-axis
#'
#' A character vector defining the order of amino acids and special variant
#' types for the heatmap y-axis.
#'
#' @export
amino_acid_order <- c(
  "codon_deletion", "frameshift", "stop_scan",
  "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
  "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
)

# Domain Fetching Functions ----------------------------------------------------

#' Fetch protein domains from Ensembl REST API
#'
#' Retrieves protein domain annotations for a given transcript from the
#' Ensembl REST API.
#'
#' @param transcript_id Character. Ensembl transcript ID (e.g., "ENST00000355451").
#' @param species Character. Species name (default: "human"). Currently unused
#'   but reserved for future multi-species support.
#' @param domain_sources Character vector. Filter to specific domain sources
#'   (e.g., c("Pfam", "SMART")). If NULL (default), returns all sources.
#' @param timeout Numeric. Request timeout in seconds (default: 30).
#'
#' @return A tibble with columns: domain, source, start, end.
#'   Returns NULL if no domains are found.
#'
#' @details
#' The function first looks up the protein ID associated with the transcript,

#' then fetches all protein features annotated as domains.
#'
#' @examples
#' \dontrun{
#' domains <- fetch_domains_ensembl("ENST00000355451")
#' domains <- fetch_domains_ensembl("ENST00000355451",
#'                                   domain_sources = c("Pfam", "SMART"))
#' }
#'
#' @export
fetch_domains_ensembl <- function(transcript_id,
                                  species = "human",
                                  domain_sources = NULL,
                                  timeout = 30) {
  # Validate input
  if (!is.character(transcript_id) || length(transcript_id) != 1) {
    stop("transcript_id must be a single character string")
  }

  # Get protein ID from transcript

url <- paste0(
    "https://rest.ensembl.org/lookup/id/", transcript_id,
    "?content-type=application/json;expand=1"
  )

  response <- httr::GET(url, httr::timeout(timeout))

  if (httr::status_code(response) != 200) {
    stop("Failed to fetch transcript info. Check transcript ID.")
  }

  transcript_info <- jsonlite::fromJSON(
    httr::content(response, "text", encoding = "UTF-8")
  )
  protein_id <- transcript_info$Translation$id

  if (is.null(protein_id)) {
    stop("No protein associated with this transcript.")
  }

  # Fetch protein features (domains)
  url_features <- paste0(
    "https://rest.ensembl.org/overlap/translation/",
    protein_id,
    "?content-type=application/json;feature=protein_feature"
  )

  response_features <- httr::GET(url_features, httr::timeout(timeout))

  if (httr::status_code(response_features) != 200) {
    stop("Failed to fetch protein features.")
  }

  features <- jsonlite::fromJSON(
    httr::content(response_features, "text", encoding = "UTF-8")
  )

  if (length(features) == 0) {
    message("No domains found for this protein.")
    return(NULL)
  }

  # Parse domain information
  domain_df <- features |>
    dplyr::select(
      domain = "description",
      source = "type",
      "start",
      "end"
    ) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$start)

  # Filter by domain sources if specified
  if (!is.null(domain_sources)) {
    domain_df <- domain_df |>
      dplyr::filter(.data$source %in% domain_sources)
  }

  return(tibble::as_tibble(domain_df))
}

# Range Merging Functions ------------------------------------------------------

#' Merge overlapping domain ranges
#'
#' Merges overlapping ranges for domains with the same name and source
#' using IRanges::reduce() for efficient interval operations.
#'
#' @param domain_df A data frame with columns: domain, source, start, end.
#'
#' @return A tibble with merged ranges (same columns as input).
#'
#' @examples
#' \dontrun{
#' domains <- fetch_domains_ensembl("ENST00000355451")
#' merged <- merge_overlapping_domains(domains)
#' }
#'
#' @export
merge_overlapping_domains <- function(domain_df) {
  if (is.null(domain_df) || nrow(domain_df) == 0) {
    return(domain_df)
  }

  merged <- domain_df |>
    dplyr::group_by(.data$domain, .data$source) |>
    dplyr::summarise(
      ranges = list(IRanges::reduce(IRanges::IRanges(.data$start, .data$end))),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      start = list(IRanges::start(.data$ranges)),
      end = list(IRanges::end(.data$ranges))
    ) |>
    dplyr::select(-"ranges") |>
    tidyr::unnest(cols = c("start", "end")) |>
    dplyr::ungroup()

  tibble::as_tibble(merged)
}

# Data Preparation Functions ---------------------------------------------------

#' Filter dataset to protein_altering variants for heatmap visualization
#'
#' Filters a dataset to keep only protein_altering variants suitable for amino acid
#' heatmap visualization. Removes custom mutator variants, synonymous variants,
#' and variants without valid protein positions.
#'
#' @param data A data frame containing variant information.
#' @param mutator_col Character. Name of the mutator column (default: "mutator").
#' @param consequence_col Character. Name of the consequence column
#'   (default: "slim_consequence").
#' @param amino_acids_col Character. Name of the amino acids column
#'   (default: "Amino_acids").
#' @param position_col Character. Name of the protein position column
#'   (default: "Protein_position").
#' @param exclude_mutators Character vector. Mutator values to exclude
#'   (default: "custom").
#' @param exclude_consequences Character vector. Consequence values to exclude
#'   (default: "synonymous").
#'
#' @return A tibble with filtered variants and AA_change column added.
#'
#' @examples
#' \dontrun{
#' protein_altering_df <- filter_protein_altering_variants(combined_dataset)
#' }
#'
#' @export
filter_protein_altering_variants <- function(data,
                                     mutator_col = "mutator",
                                     consequence_col = "slim_consequence",
                                     amino_acids_col = "Amino_acids",
                                     position_col = "Protein_position",
                                     exclude_mutators = "custom",
                                     exclude_consequences = "synonymous") {
  result <- data

  # Filter out excluded mutators if column exists

  if (mutator_col %in% names(data)) {
    result <- result |>
      dplyr::filter(!.data[[mutator_col]] %in% exclude_mutators)
  }

  # Filter out excluded consequences if column exists
  if (consequence_col %in% names(data)) {
    result <- result |>
      dplyr::filter(!.data[[consequence_col]] %in% exclude_consequences)
  }

  # Create AA_change and ensure numeric position
  result <- result |>
    dplyr::mutate(
      AA_change = paste0(
        substr(.data[[amino_acids_col]], 1, 1),
        .data[[position_col]],
        substr(.data[[amino_acids_col]], 3, 3)
      ),
      Protein_position = as.numeric(.data[[position_col]])
    ) |>
    dplyr::filter(!is.na(.data$Protein_position))

  tibble::as_tibble(result)
}

#' Prepare data for amino acid heatmap visualization
#'
#' Filters and aggregates variant data for creating an amino acid substitution
#' heatmap. Groups by amino acid change and computes mean scores and FDR,
#' then classifies variants as enriched, depleted, or unchanged.
#'
#' @param data A data frame containing variant information.
#' @param score_col Character. Name of the score column (default: "adj_score_continuous").
#' @param fdr_col Character. Name of the FDR column (default: "FDR_continuous").
#' @param fdr_threshold Numeric. FDR threshold for significance (default: 0.01).
#' @param amino_acids_col Character. Name of the column containing amino acid
#'   change in format "X/Y" (default: "Amino_acids").
#' @param position_col Character. Name of the protein position column
#'
#' @return A tibble with columns: AA_change, previous, new_aa, Protein_position,
#'   mean_adj_score, mean_FDR, functional_classification.
#'
#' @details
#' The function:
#' 1. Creates AA_change identifier (e.g., "A123V")
#' 2. Groups by unique amino acid changes and computes mean score/FDR
#' 3. Classifies based on FDR threshold and score direction
#' 4. Converts special characters: "*" -> "stop_scan", "-" -> "codon_deletion", "X" -> "frameshift"
#' 5. Orders new_aa as factor using amino_acid_order
#'
#' @examples
#' \dontrun{
#' plot_data <- prepare_amino_acid_heatmap_data(variant_results)
#' }
#'
#' @export
prepare_amino_acid_heatmap_data <- function(data,
                                            score_col = "adj_score_continuous",
                                            fdr_col = "FDR_continuous",
                                            fdr_threshold = 0.01,
                                            amino_acids_col = "Amino_acids",
                                            position_col = "Protein_position") {
  # Create AA_change identifier
  plot_df <- data |>
    dplyr::mutate(
      AA_change = paste0(
        substr(.data[[amino_acids_col]], 1, 1),
        .data[[position_col]],
        substr(.data[[amino_acids_col]], 3, 3)
      ),
      Protein_position = as.numeric(.data[[position_col]])
    ) |>
    dplyr::filter(!is.na(.data$Protein_position))

  # Group and summarise
  plot_df <- plot_df |>
    dplyr::group_by(.data$AA_change, .data[[amino_acids_col]], .data$Protein_position) |>
    dplyr::summarise(
      mean_adj_score = if (all(is.na(.data[[score_col]]))) {
        NA_real_
      } else {
        mean(.data[[score_col]], na.rm = TRUE)
      },
      mean_FDR = if (all(is.na(.data[[fdr_col]]))) {
        NA_real_
      } else {
        mean(.data[[fdr_col]], na.rm = TRUE)
      },
      .groups = "drop"
    ) |>
    dplyr::filter(!is.na(.data$mean_adj_score) & !is.na(.data$mean_FDR))

  # Rename the amino acids column for processing
  plot_df <- plot_df |>
    dplyr::rename(Amino_acids = dplyr::all_of(amino_acids_col))

  # Classify variants
  plot_df <- plot_df |>
    dplyr::mutate(
      functional_classification = dplyr::case_when(
        .data$mean_FDR < fdr_threshold & .data$mean_adj_score < 0 ~ "depleted",
        .data$mean_FDR < fdr_threshold & .data$mean_adj_score > 0 ~ "enriched",
        TRUE ~ "unchanged"
      )
    )

  # Separate amino acids and convert special characters
  plot_df <- plot_df |>
    tidyr::separate(
      .data$Amino_acids,
      into = c("previous", "new_aa"),
      sep = "/"
    ) |>
    dplyr::mutate(
      new_aa = dplyr::case_when(
        .data$new_aa == "*" ~ "stop_scan",
        .data$new_aa == "-" ~ "codon_deletion",
        .data$new_aa == "X" ~ "frameshift",
        TRUE ~ .data$new_aa
      )
    ) |>
    dplyr::mutate(new_aa = factor(.data$new_aa, levels = amino_acid_order))

  return(tibble::as_tibble(plot_df))
}

# Plotting Functions -----------------------------------------------------------

#' Plot protein domain track
#'
#' Creates a ggplot2 visualization of protein domains as colored rectangles.
#' Designed to be combined with other plots using patchwork.
#'
#' @param domain_df A data frame with columns: domain, source, start, end.
#' @param x_range Numeric vector of length 2 specifying x-axis limits.
#'   If NULL, determined from data.
#' @param show_legend Logical. Whether to show the domain legend (default: TRUE).
#' @param domain_colors Named character vector of colors for domains.
#'   If NULL, uses default ggplot2 colors.
#' @param show_x_axis Logical. Whether to show the x-axis (default: FALSE).
#'   Set to TRUE if using as standalone plot.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' domains <- fetch_domains_ensembl("ENST00000355451")
#' merged <- merge_overlapping_domains(domains)
#' plot_domain_track(merged)
#' }
#'
#' @export
plot_domain_track <- function(domain_df,
                              x_range = NULL,
                              show_legend = TRUE,
                              domain_colors = NULL,
                              show_x_axis = FALSE) {
  if (is.null(domain_df) || nrow(domain_df) == 0) {
    # Return empty plot if no domains
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No domains found")
    return(p)
  }

  # Determine x_range if not provided
  if (is.null(x_range)) {
    x_range <- c(min(domain_df$start), max(domain_df$end))
  }

  p <- ggplot2::ggplot(domain_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = 0,
        ymax = 1,
        fill = .data$domain
      ),
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::ylab("")

  # Apply custom colors if provided
  if (!is.null(domain_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = domain_colors)
  }

  # Configure axis and legend
  legend_position <- if (show_legend) "right" else "none"

  axis_elements <- if (show_x_axis) {
    list(
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        legend.position = legend_position,
        legend.title = ggplot2::element_blank()
      ),
      ggplot2::xlab("Protein position")
    )
  } else {
    list(
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        legend.position = legend_position,
        legend.title = ggplot2::element_blank()
      ),
      ggplot2::xlab("")
    )
  }

  p <- p + axis_elements

  return(p)
}

#' Plot amino acid substitution heatmap
#'
#' Creates a ggplot2 heatmap showing amino acid substitutions colored by
#' functional classification (enriched, depleted, unchanged).
#'
#' @param data A data frame prepared by \code{\link{prepare_amino_acid_heatmap_data}}.
#'   Must have columns: Protein_position, new_aa, functional_classification.
#' @param x_range Numeric vector of length 2 specifying x-axis limits.
#'   If NULL, determined from data.
#' @param classification_colors Named character vector of colors for
#'   functional classifications. If NULL, uses \code{classification_heatmap_colours}.
#' @param show_x_axis Logical. Whether to show the x-axis (default: TRUE).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_data <- prepare_amino_acid_heatmap_data(variant_results)
#' plot_amino_acid_heatmap(plot_data)
#' }
#'
#' @export
plot_amino_acid_heatmap <- function(data,
                                    x_range = NULL,
                                    classification_colors = NULL,
                                    show_x_axis = TRUE) {
  if (is.null(data) || nrow(data) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No data to plot")
    return(p)
  }

  # Use default colors if not provided
  if (is.null(classification_colors)) {
    classification_colors <- classification_heatmap_colours
  }

  # Determine x_range if not provided
  if (is.null(x_range)) {
    x_range <- range(data$Protein_position, na.rm = TRUE)
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$Protein_position, y = .data$new_aa)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$functional_classification),
      color = "black",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(values = classification_colors) +
    ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::ylab("") +
    ggplot2::labs(fill = "Classification")

  # Configure x-axis
  if (show_x_axis) {
    p <- p +
      ggplot2::xlab("Protein position") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8, colour = "black"),
        axis.title.y = ggplot2::element_text(size = 12, colour = "black")
      )
  } else {
    p <- p +
      ggplot2::xlab("") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 8, colour = "black"),
        axis.title.y = ggplot2::element_text(size = 12, colour = "black")
      )
  }

  return(p)
}

#' Plot combined domain track and amino acid heatmap
#'
#' High-level function that creates a combined visualization with protein
#' domains shown above an amino acid substitution heatmap. Uses patchwork
#' to align x-axes.
#'
#' @param data A data frame containing variant data with columns for amino acid
#'   changes, positions, scores, and FDR values. Will be processed by
#'   \code{\link{prepare_amino_acid_heatmap_data}} if not already prepared.
#' @param transcript_id Character. Ensembl transcript ID for fetching domains.
#'   Required if \code{domain_df} is NULL.
#' @param domain_df A data frame of domains (output from
#'   \code{\link{fetch_domains_ensembl}}). If NULL and transcript_id provided,
#'   domains will be fetched automatically.
#' @param merge_domains Logical. Whether to merge overlapping domains
#'   (default: TRUE).
#' @param prepared Logical. Whether \code{data} has already been processed by
#'   \code{\link{prepare_amino_acid_heatmap_data}} (default: FALSE).
#' @param height_ratio Numeric vector of length 2. Relative heights of domain
#'   track and heatmap (default: c(1, 15)).
#' @param domain_sources Character vector. Filter domains to specific sources.
#' @param classification_colors Named character vector for heatmap colors.
#' @param domain_colors Named character vector for domain colors.
#' @param ... Additional arguments passed to
#'   \code{\link{prepare_amino_acid_heatmap_data}}.
#'
#' @return A patchwork object combining domain track and heatmap.
#'
#' @examples
#' \dontrun{
#' # With automatic domain fetching
#' p <- plot_domain_heatmap(
#'   data = variant_results,
#'   transcript_id = "ENST00000355451"
#' )
#'
#' # With pre-fetched domains
#' domains <- fetch_domains_ensembl("ENST00000355451")
#' p <- plot_domain_heatmap(
#'   data = variant_results,
#'   domain_df = domains
#' )
#'
#' # Save to file
#' ggplot2::ggsave("domain_heatmap.png", p, width = 8, height = 7, dpi = 600)
#' }
#'
#' @export
plot_domain_heatmap <- function(data,
                                transcript_id = NULL,
                                domain_df = NULL,
                                merge_domains = TRUE,
                                prepared = FALSE,
                                height_ratio = c(1, 15),
                                domain_sources = NULL,
                                classification_colors = NULL,
                                domain_colors = NULL,
                                ...) {
  # Fetch domains if needed
  if (is.null(domain_df)) {
    if (is.null(transcript_id)) {
      stop("Either 'domain_df' or 'transcript_id' must be provided.")
    }
    domain_df <- fetch_domains_ensembl(
      transcript_id,
      domain_sources = domain_sources
    )
  }

  # Merge domains if requested
  if (merge_domains && !is.null(domain_df)) {
    domain_df <- merge_overlapping_domains(domain_df)
  }

  # Prepare heatmap data if not already done
  if (!prepared) {
    plot_data <- prepare_amino_acid_heatmap_data(data, ...)
  } else {
    plot_data <- data
  }

  # Determine shared x_range
  x_range <- range(plot_data$Protein_position, na.rm = TRUE)

  # Create individual plots
  domain_track <- plot_domain_track(
    domain_df,
    x_range = x_range,
    show_legend = TRUE,
    domain_colors = domain_colors,
    show_x_axis = FALSE
  )

  heatmap <- plot_amino_acid_heatmap(
    plot_data,
    x_range = x_range,
    classification_colors = classification_colors,
    show_x_axis = TRUE
  )

  # Combine with patchwork
  combined <- domain_track / heatmap +
    patchwork::plot_layout(heights = height_ratio, axes = "collect_x")

  return(combined)
}
