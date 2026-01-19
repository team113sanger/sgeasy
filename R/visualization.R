#' Create sample comparison scatter plots
#'
#' Generates a matrix of pairwise scatter plots for samples, useful for
#' visualizing replicate correlation and batch effects.
#'
#' @param data A matrix or DESeqTransform object of transformed counts
#'   (e.g., rlog). If a DESeqTransform, the assay will be extracted.
#' @param sample_names Character vector of sample names to include in the plot.
#' @param point_size Size of scatter plot points (default: 0.4).
#' @param line_color Color for regression line (default: "blue").
#' @param title Plot title (default: "Regularized-log Transformed Read Count").
#'
#' @return Invisibly returns the ggmatrix plot object. The plot is also
#'   printed to the current graphics device.
#'
#' @details
#' This function creates a pairs plot using [GGally::ggpairs()] showing
#' pairwise scatter plots between selected samples. Each panel includes:
#' \itemize{
#'   \item Scatter points showing variant abundance in each sample pair
#'   \item Linear regression line for correlation visualization
#'   \item Correlation coefficients in the upper triangle
#' }
#'
#' This is useful for assessing replicate consistency and identifying
#' potential batch effects or outlier samples.
#'
#' @examples
#' \dontrun{
#' # Plot Day4 replicates
#' plot_sample_scatter(
#'   data = assay(results$rlog),
#'   sample_names = c("SAMPLE_A_4_pel", "SAMPLE_B_4_pel")
#' )
#'
#' # Or pass the rlog object directly
#' plot_sample_scatter(
#'   data = results$rlog,
#'   sample_names = c("SAMPLE_A_4_pel", "SAMPLE_B_4_pel")
#' )
#' }
#'
#' @seealso [GGally::ggpairs()]
#' @export
plot_sample_scatter <- function(data,
                                sample_names,
                                point_size = 0.4,
                                line_color = "blue",
                                title = "Regularized-log Transformed Read Count") {

  # Extract assay if DESeqTransform object
 if (inherits(data, "DESeqTransform")) {
    data <- SummarizedExperiment::assay(data)
  }

  subset_data <- data |>
    tibble::as_tibble(rownames = "SEQUENCE") |>
    dplyr::select(dplyr::all_of(sample_names))

  p <- GGally::ggpairs(
    subset_data,
    lower = list(
      continuous = GGally::wrap(
        .ggpairs_scatter,
        pts = list(size = point_size, colour = "black"),
        smt = list(method = "lm", se = FALSE, size = 0.2, colour = line_color)
      )
    ),
    title = title
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 13, face = "bold")
    )

  print(p)
  invisible(p)
}


#' Custom ggpairs scatter panel function
#'
#' Internal function providing customized scatter plot panels for ggpairs.
#'
#' @param data Data passed by ggpairs.
#' @param mapping Aesthetic mapping passed by ggpairs.
#' @param pts List of arguments for [ggplot2::geom_point()].
#' @param smt List of arguments for [ggplot2::geom_smooth()].
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object for the scatter panel.
#'
#' @keywords internal
.ggpairs_scatter <- function(data, mapping, pts = list(), smt = list(), ...) {
  ggplot2::ggplot(data = data, mapping = mapping) +
    do.call(ggplot2::geom_point, pts) +
    do.call(ggplot2::geom_smooth, smt)
}


#' Default color palette for consequence types
#'
#' A named character vector mapping slim_consequence categories to colors.
#'
#' @format A named character vector.
#'
#' @examples
#' consequence_colours
#'
#' @export
consequence_colours <- c(
  "synonymous" = "steelblue3",
  "stop_gained" = "firebrick2",

  "missense" = "seagreen",
  "codon_deletion" = "orchid4",
  "frameshift" = "gold1",
  "intron" = "lightsteelblue4",
  "splice_donor" = "hotpink2",
  "splice_acceptor" = "darkorange1",
  "UTR" = "royalblue4"
)


#' Default shape palette for functional classifications
#'
#' A named numeric vector mapping functional classifications to point shapes.
#'
#' @format A named numeric vector.
#'
#' @examples
#' classification_shapes
#'
#' @export
classification_shapes <- c(
  "unchanged" = 0,
  "depleted" = 6,
  "enriched" = 2
)


#' Plot z-scores for a single condition contrast
#'
#' Creates a waterfall-style scatter plot showing variant z-scores along
#' genomic coordinates, colored by consequence type and shaped by
#' functional classification.
#'
#' @param data A data frame with adj_score, functional_classification,
#'   slim_consequence, vcf_pos, and FDR columns.
#' @param suffix The contrast suffix (e.g., "condition_Day7_vs_Day4").
#' @param targeton_id Identifier for the targeton (used in plot title).
#' @param colours Named vector of colors for slim_consequence values
#'   (default: [consequence_colours]).
#' @param shapes Named vector of shapes for functional_classification values
#'   (default: [classification_shapes]).
#' @param fdr_threshold FDR threshold for alpha transparency (default: 0.01).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' p <- plot_condition_scores(
#'   data = processed_data,
#'   suffix = "condition_Day7_vs_Day4",
#'   targeton_id = "GENE1"
#' )
#' }
#'
#' @export
plot_condition_scores <- function(data,
                                  suffix,
                                  targeton_id,
                                  colours = consequence_colours,
                                  shapes = classification_shapes,
                                  fdr_threshold = 0.01) {
  score_col <- paste0("adj_score_", suffix)
  class_col <- paste0("functional_classification_", suffix)
  fdr_col <- paste0("FDR_", suffix)

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$vcf_pos,
      y = .data[[score_col]],
      shape = .data[[class_col]],
      colour = .data$slim_consequence,
      fill = .data$slim_consequence
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(alpha = .data[[fdr_col]] <= fdr_threshold),
      size = 1.5
    ) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::scale_alpha_discrete(
      labels = c(
        paste0("FDR > ", fdr_threshold),
        paste0("FDR <= ", fdr_threshold)
      )
    ) +
    ggplot2::xlab("GRCh38 coordinate") +
    ggplot2::ylab(paste0("z-score (", suffix, ")")) +
    ggplot2::ggtitle(paste("z-score Plot:", targeton_id, "-", suffix)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1,
        colour = "black", size = 10
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = ggplot2::margin()
    )
}


#' Generate plots for all condition contrasts
#'
#' Creates z-score plots for all contrasts found in the data.
#'
#' @param data A data frame with processed screen data.
#' @param targeton_id Identifier for the targeton.
#' @param colours Named vector of colors (default: [consequence_colours]).
#' @param shapes Named vector of shapes (default: [classification_shapes]).
#'
#' @return A named list of ggplot objects, one per contrast.
#'
#' @examples
#' \dontrun{
#' plots <- plot_all_conditions(processed_data, "GENE1")
#' }
#'
#' @export
plot_all_conditions <- function(data,
                                targeton_id,
                                colours = consequence_colours,
                                shapes = classification_shapes) {
  suffixes <- names(data) |>
    stringr::str_subset("^adj_score_") |>
    stringr::str_remove("^adj_score_")

  plots <- purrr::map(
    suffixes,
    ~ plot_condition_scores(data, .x, targeton_id, colours, shapes)
  )
  names(plots) <- suffixes
  plots
}


#' Write waterfall plots to files
#'
#' Saves a nested list of plots to PNG files.
#'
#' @param plot_list A named list of plot lists (typically from
#'   [plot_all_conditions()] applied to multiple targetons).
#' @param output_dir Directory to write plots to (default: "results").
#' @param width Plot width in pixels (default: 900).
#' @param height Plot height in pixels (default: 600).
#'
#' @return Invisibly returns NULL. Side effect: writes PNG files.
#'
#' @examples
#' \dontrun{
#' write_waterfall_plots(plot_list, output_dir = "figures")
#' }
#'
#' @export
write_waterfall_plots <- function(plot_list,
                                  output_dir = "results",
                                  width = 900,
                                  height = 600) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  plot_list |>
    purrr::list_flatten() |>
    purrr::iwalk(~ {
      grDevices::png(
        file.path(output_dir, paste0(.y, ".png")),
        width = width,
        height = height
      )
      print(.x)
      grDevices::dev.off()
    })

  invisible(NULL)
}
