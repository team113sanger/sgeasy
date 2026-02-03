# Amino Acid Heatmap Functions
# ============================
# Data objects, data preparation, and plotting functions for amino acid
# substitution heatmaps with protein domain and pLDDT tracks.

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

#' AlphaFold pLDDT confidence colour scheme
#'
#' A named character vector of colours for pLDDT confidence categories,
#' matching the standard AlphaFold colour scheme.
#'
#' @export
plddt_colours <- c(
  "H" = "#0053D6",
  "M" = "#65CBF3",
  "L" = "#FFDB13",
  "D" = "#FF7D45"
)

# Data Preparation Functions ---------------------------------------------------

#' Add AA_change column and convert protein position to numeric
#'
#' Internal helper function that creates the AA_change identifier (e.g., "A123V")
#' from amino acid and position columns, converts position to numeric, and
#' filters out rows with invalid positions.
#'
#' @param data A data frame containing variant information.
#' @param amino_acids_col Character. Name of the amino acids column.
#' @param position_col Character. Name of the protein position column.
#'
#' @return A data frame with AA_change column added, Protein_position as numeric,
#'   and rows with NA positions removed.
#'
#' @keywords internal
add_aa_change <- function(data, amino_acids_col, position_col) {
  data |>
    dplyr::mutate(
      AA_change = paste0(
        substr(.data[[amino_acids_col]], 1, 1),
        .data[[position_col]],
        substr(.data[[amino_acids_col]], 3, 3)
      ),
      Protein_position = as.numeric(.data[[position_col]])
    ) |>
    dplyr::filter(!is.na(.data$Protein_position))
}

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
  result <- add_aa_change(result, amino_acids_col, position_col)

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
  plot_df <- add_aa_change(data, amino_acids_col, position_col)

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

  tibble::as_tibble(plot_df)
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

  p + axis_elements
}

#' Plot pLDDT confidence track
#'
#' Creates a ggplot2 track showing per-residue AlphaFold pLDDT confidence
#' scores, coloured by confidence category using the standard AlphaFold
#' colour scheme.
#'
#' @param plddt_df A data frame as returned by
#'   \code{\link{fetch_plddt_alphafold}}, with columns: position, plddt,
#'   confidence_category.
#' @param x_range Numeric vector of length 2 specifying x-axis limits.
#'   If NULL, determined from data.
#' @param show_legend Logical. Whether to show the colour legend (default: TRUE).
#' @param plddt_colors Named character vector of colours for confidence
#'   categories (H, M, L, D). If NULL, uses \code{plddt_colours}.
#' @param show_x_axis Logical. Whether to show the x-axis (default: FALSE).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plddt <- fetch_plddt_alphafold("Q8IVD9")
#' plot_plddt_track(plddt)
#' }
#'
#' @export
plot_plddt_track <- function(plddt_df,
                             x_range = NULL,
                             show_legend = TRUE,
                             plddt_colors = NULL,
                             show_x_axis = FALSE) {
  if (is.null(plddt_df) || nrow(plddt_df) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No pLDDT data found")
    return(p)
  }

  if (is.null(plddt_colors)) {
    plddt_colors <- plddt_colours
  }

  if (is.null(x_range)) {
    x_range <- range(plddt_df$position, na.rm = TRUE)
  }

  # Order categories for legend
  plddt_df$confidence_category <- factor(
    plddt_df$confidence_category,
    levels = c("H", "M", "L", "D")
  )

  p <- ggplot2::ggplot(plddt_df) +
    ggplot2::geom_tile(
      ggplot2::aes(
        x = .data$position,
        y = 0.5,
        fill = .data$confidence_category
      ),
      height = 1
    ) +
    ggplot2::scale_fill_manual(
      values = plddt_colors,
      labels = c(
        "H" = "Very high (>90)",
        "M" = "Confident (70-90)",
        "L" = "Low (50-70)",
        "D" = "Very low (<50)"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::ylab("") +
    ggplot2::labs(fill = "pLDDT")

  legend_position <- if (show_legend) "right" else "none"

  if (show_x_axis) {
    p <- p +
      ggplot2::xlab("Protein position") +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        legend.position = legend_position
      )
  } else {
    p <- p +
      ggplot2::xlab("") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        legend.position = legend_position
      )
  }

  p
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

  p
}

#' Plot combined domain track and amino acid heatmap
#'
#' High-level function that creates a combined visualization with
#' protein domains and optionally pLDDT confidence scores shown
#' above an amino acid substitution heatmap. Uses patchwork to
#' align x-axes.
#'
#' @param data A data frame containing variant data with columns
#'   for amino acid changes, positions, scores, and FDR values.
#'   Will be processed by
#'   \code{\link{prepare_amino_acid_heatmap_data}} if not already
#'   prepared.
#' @param uniprot_acc Character. UniProt accession for fetching
#'   domains (e.g., "P04637"). Also used to fetch pLDDT scores
#'   when \code{plddt_df = "auto"}.
#' @param transcript_id Character. Ensembl transcript ID for
#'   fetching domains (e.g., "ENST00000355451"). Used with
#'   domain_api = "ensembl".
#' @param domain_api Character. Which API to use for fetching
#'   domains: "ted" (default), "interpro", or "ensembl".
#' @param domain_df A data frame of domains. If provided, skips
#'   API fetch.
#' @param plddt_df A data frame of pLDDT scores as returned by
#'   \code{\link{fetch_plddt_alphafold}}, or the string "auto"
#'   to fetch from AlphaFold using \code{uniprot_acc}. If NULL
#'   (default), no pLDDT track is shown.
#' @param merge_domains Logical. Whether to merge overlapping
#'   domains (default: TRUE).
#' @param prepared Logical. Whether \code{data} has already been
#'   processed by
#'   \code{\link{prepare_amino_acid_heatmap_data}}
#'   (default: FALSE).
#' @param height_ratio Numeric vector. Relative heights of tracks
#'   and heatmap. Length depends on number of tracks: 2 elements
#'   for domain + heatmap, 3 elements for pLDDT + domain +
#'   heatmap. If NULL, uses sensible defaults.
#' @param domain_sources Character vector. Filter domains to
#'   specific sources (e.g., c("pfam", "smart") for InterPro).
#' @param consensus_levels Character vector. For TED API only:
#'   filter by consensus level ("high", "medium", "low").
#' @param classification_colors Named character vector for
#'   heatmap colors.
#' @param domain_colors Named character vector for domain colors.
#' @param plddt_colors Named character vector for pLDDT
#'   confidence category colors (H, M, L, D).
#' @param ... Additional arguments passed to
#'   \code{\link{prepare_amino_acid_heatmap_data}}.
#'
#' @return A patchwork object combining tracks and heatmap.
#'
#' @examples
#' \dontrun{
#' # With domain track only
#' p <- plot_domain_heatmap(
#'   data = variant_results,
#'   uniprot_acc = "P04637",
#'   domain_api = "ted"
#' )
#'
#' # With pLDDT track (auto-fetch from AlphaFold)
#' p <- plot_domain_heatmap(
#'   data = variant_results,
#'   uniprot_acc = "P04637",
#'   plddt_df = "auto"
#' )
#'
#' # With pre-fetched pLDDT data
#' plddt <- fetch_plddt_alphafold("P04637")
#' p <- plot_domain_heatmap(
#'   data = variant_results,
#'   uniprot_acc = "P04637",
#'   plddt_df = plddt
#' )
#'
#' # Save to file
#' ggplot2::ggsave(
#'   "domain_heatmap.png", p,
#'   width = 8, height = 7, dpi = 600
#' )
#' }
#'
#' @export
plot_domain_heatmap <- function(data,
                                uniprot_acc = NULL,
                                transcript_id = NULL,
                                domain_api = "ted",
                                domain_df = NULL,
                                plddt_df = NULL,
                                merge_domains = TRUE,
                                prepared = FALSE,
                                height_ratio = NULL,
                                domain_sources = NULL,
                                consensus_levels = NULL,
                                classification_colors = NULL,
                                domain_colors = NULL,
                                plddt_colors = NULL,
                                ...) {
  # Validate domain_api
  domain_api <- match.arg(domain_api, c("ted", "interpro", "ensembl"))

  # Fetch domains if needed
  if (is.null(domain_df)) {
    if (domain_api == "ted") {
      if (is.null(uniprot_acc)) {
        stop("'uniprot_acc' required for TED API.")
      }
      domain_df <- fetch_domains_ted(
        uniprot_acc,
        consensus_levels = consensus_levels
      )
    } else if (domain_api == "interpro") {
      if (is.null(uniprot_acc)) {
        stop("'uniprot_acc' required for InterPro API.")
      }
      domain_df <- fetch_domains_interpro(
        uniprot_acc,
        domain_sources = domain_sources
      )
    } else if (domain_api == "ensembl") {
      if (is.null(transcript_id)) {
        stop("'transcript_id' required for Ensembl API.")
      }
      domain_df <- fetch_domains_ensembl(
        transcript_id,
        domain_sources = domain_sources
      )
    }
  }

  # Merge domains if requested
  if (merge_domains && !is.null(domain_df)) {
    domain_df <- merge_overlapping_domains(domain_df)
  }

  # Fetch pLDDT if requested
  if (identical(plddt_df, "auto")) {
    if (is.null(uniprot_acc)) {
      stop("'uniprot_acc' required to auto-fetch pLDDT.")
    }
    plddt_df <- fetch_plddt_alphafold(uniprot_acc)
  }

  # Prepare heatmap data if not already done
  if (!prepared) {
    plot_data <- prepare_amino_acid_heatmap_data(data, ...)
  } else {
    plot_data <- data
  }

  # Determine shared x_range
  x_range <- range(
    plot_data$Protein_position, na.rm = TRUE
  )

  # Create individual plots
  has_plddt <- !is.null(plddt_df) && nrow(plddt_df) > 0

  if (has_plddt) {
    plddt_track <- plot_plddt_track(
      plddt_df,
      x_range = x_range,
      show_legend = TRUE,
      plddt_colors = plddt_colors,
      show_x_axis = FALSE
    )
  }

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
  if (has_plddt) {
    if (is.null(height_ratio)) {
      height_ratio <- c(1, 1, 15)
    }
    combined <- plddt_track / domain_track / heatmap +
      patchwork::plot_layout(
        heights = height_ratio,
        axes = "collect_x"
      )
  } else {
    if (is.null(height_ratio)) {
      height_ratio <- c(1, 15)
    }
    combined <- domain_track / heatmap +
      patchwork::plot_layout(
        heights = height_ratio,
        axes = "collect_x"
      )
  }

  combined
}
