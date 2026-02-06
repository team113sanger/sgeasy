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
  "stop_lost" = "purple3",
  "start_lost" = "darkred",
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
  "enriched" = 3
)


#' Shared scale and theme layers for SGE scatter plots
#'
#' Internal helper that returns a list of ggplot2 layers shared by
#' \code{\link{plot_condition_scores}} and \code{\link{plot_gene_level_waterfall}}.
#'
#' @param colours Named colour vector for consequence types.
#' @param shapes Named shape vector for functional classifications.
#' @param fdr_threshold FDR threshold for alpha legend labels.
#'
#' @return A list of ggplot2 scale and theme objects.
#'
#' @keywords internal
.sge_scatter_scales <- function(colours, shapes, fdr_threshold) {
  list(
    ggplot2::scale_colour_manual(values = colours),
    ggplot2::scale_fill_manual(values = colours),
    ggplot2::scale_shape_manual(values = shapes),
    ggplot2::scale_alpha_discrete(
      labels = c(paste0("FDR > ", fdr_threshold),
                 paste0("FDR <= ", fdr_threshold))
    ),
    ggplot2::theme_classic(),
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
  )
}


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
    ggplot2::xlab("GRCh38 coordinate") +
    ggplot2::ylab(paste0("z-score (", suffix, ")")) +
    ggplot2::ggtitle(paste("z-score Plot:", targeton_id, "-", suffix)) +
    .sge_scatter_scales(colours, shapes, fdr_threshold)
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


#' Save a list of plots to PNG files in a directory
#'
#' Internal helper that creates the output directory (if needed) and
#' iterates over a named list, writing each element as a PNG.
#'
#' @param plot_list A named list of plot objects.
#' @param output_dir Directory to write files to.
#' @param suffix Character string appended to each name before
#'   the \code{.png} extension.
#' @param width Plot width in pixels.
#' @param height Plot height in pixels.
#' @param plot_fn A function that renders a single element from
#'   \code{plot_list} to the current graphics device.
#'
#' @return \code{invisible(NULL)}.
#'
#' @keywords internal
.save_plots_to_dir <- function(plot_list, output_dir,
                               suffix, width, height,
                               plot_fn) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  purrr::iwalk(plot_list, ~ {
    grDevices::png(
      file.path(output_dir, paste0(.y, suffix, ".png")),
      width = width, height = height
    )
    plot_fn(.x)
    grDevices::dev.off()
  })
  invisible(NULL)
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
  .save_plots_to_dir(
    purrr::list_flatten(plot_list),
    output_dir, "", width, height,
    function(p) print(p)
  )
}


#' Create sample distance matrix heatmap
#'
#' Computes a sample-to-sample distance matrix from regularized log-transformed
#' data and displays it as a heatmap.
#'
#' @param rlog_data A DESeqTransform object (output from `DESeq2::rlog()`)
#'   or a matrix of transformed counts.
#' @param color_palette Color palette for the heatmap. Default uses reversed
#'   Blues palette.
#'
#' @return A pheatmap object (invisibly). The heatmap is also printed.
#'
#' @examples
#' \dontrun{
#' sample_distance_matrix(results$rlog)
#' }
#'
#' @export
sample_distance_matrix <- function(rlog_data,
                                   color_palette = NULL) {
  if (inherits(rlog_data, "DESeqTransform")) {
    mat <- SummarizedExperiment::assay(rlog_data)
  } else {
    mat <- rlog_data
  }

  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "Blues"))
    )(255)
  }

  dist_matrix <- as.matrix(stats::dist(t(mat)))

  p <- pheatmap::pheatmap(
    dist_matrix,
    trace = "none",
    col = color_palette,
    silent = TRUE
  )

  print(p)
  invisible(p)
}


#' Build PCA plot from rlog data
#'
#' Creates a PCA plot from regularized log-transformed data, colored by
#' condition and shaped by replicate.
#'
#' @param rlog_data A DESeqTransform object (output from `DESeq2::rlog()`).
#' @param intgroup Character vector of variables to use for grouping
#'   (default: c("condition", "replicate")).
#' @param point_size Size of points in the plot (default: 3).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' pca_plot <- build_sample_pca(results$rlog)
#' }
#'
#' @export
build_sample_pca <- function(rlog_data,
                             intgroup = c("condition", "replicate"),
                             point_size = 3) {
  pca_data <- DESeq2::plotPCA(rlog_data, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  # Build aesthetic mapping based on available intgroup variables
  if (length(intgroup) >= 2 && all(intgroup %in% names(pca_data))) {
    p <- ggplot2::ggplot(
      pca_data,
      ggplot2::aes(
        x = .data$PC1,
        y = .data$PC2,
        color = .data[[intgroup[1]]],
        shape = .data[[intgroup[2]]]
      )
    )
  } else if (length(intgroup) >= 1 && intgroup[1] %in% names(pca_data)) {
    p <- ggplot2::ggplot(
      pca_data,
      ggplot2::aes(
        x = .data$PC1,
        y = .data$PC2,
        color = .data[[intgroup[1]]]
      )
    )
  } else {
    p <- ggplot2::ggplot(
      pca_data,
      ggplot2::aes(x = .data$PC1, y = .data$PC2)
    )
  }

  p +
    ggplot2::geom_point(size = point_size) +
    ggplot2::xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic()
}


#' Save distance matrix plots to files
#'
#' Saves multiple sample distance matrix heatmaps to PNG files.
#'
#' @param rlog_list A named list of DESeqTransform objects or matrices.
#' @param output_dir Directory to save plots (default: "results/qc").
#' @param width Plot width in pixels (default: 900).
#' @param height Plot height in pixels (default: 600).
#'
#' @return Invisibly returns NULL. Side effect: writes PNG files.
#'
#' @examples
#' \dontrun{
#' rlog_results <- map(deseq_results, pluck, "rlog")
#' dist_matrices <- map(rlog_results, sample_distance_matrix)
#' plot_distance_matrices(dist_matrices, output_dir = "results/qc")
#' }
#'
#' @export
plot_distance_matrices <- function(rlog_list,
                                   output_dir = "results/qc",
                                   width = 900,
                                   height = 600) {
  .save_plots_to_dir(
    rlog_list, output_dir,
    "_sample_distance_matrix", width, height,
    function(x) sample_distance_matrix(x)
  )
}


#' Save PCA plots to files
#'
#' Saves multiple sample PCA plots to PNG files.
#'
#' @param rlog_list A named list of DESeqTransform objects.
#' @param output_dir Directory to save plots (default: "results/qc").
#' @param width Plot width in pixels (default: 900).
#' @param height Plot height in pixels (default: 600).
#' @param intgroup Character vector of variables to use for grouping
#'   (default: c("condition", "replicate")).
#'
#' @return Invisibly returns NULL. Side effect: writes PNG files.
#'
#' @examples
#' \dontrun{
#' rlog_results <- map(deseq_results, pluck, "rlog")
#' plot_screen_pcas(rlog_results, output_dir = "results/qc")
#' }
#'
#' @export
plot_screen_pcas <- function(rlog_list,
                             output_dir = "results/qc",
                             width = 900,
                             height = 600,
                             intgroup = c("condition", "replicate")) {
  .save_plots_to_dir(
    rlog_list, output_dir,
    "_pca", width, height,
    function(x) print(build_sample_pca(x, intgroup = intgroup))
  )
}




#' Plot position effect ratios
#'
#' Creates a scatter plot showing the ratio of counts at a specified timepoint
#' versus Day0 across genomic coordinates, colored by consequence type with
#' a loess smoothing line.
#'
#' @param pos_ratio A data frame from [calculate_position_effect()] containing
#'   vcf_pos, ratio, and slim_consequence columns.
#' @param timepoint Character string specifying the timepoint (used in y-axis
#'   label, e.g., "Day10").
#' @param colours Named vector of colors for slim_consequence values
#'
#' @return A ggplot object (invisibly). The plot is also printed.
#'
#' @examples
#' \dontrun{
#' pos_effect <- calculate_position_effect(count_data, annotation, "Day10")
#' position_effect_plot(pos_effect, "Day10")
#' }
#'
#' @seealso [calculate_position_effect()]
#' @export
position_effect_plot <- function(pos_ratio,
                          timepoint,
                          colours = consequence_colours) {
  p <- ggplot2::ggplot(
    pos_ratio,
    ggplot2::aes(
      x = .data$vcf_pos,
      y = .data$ratio,
      color = .data$slim_consequence
    )
  ) +
    ggplot2::scale_y_continuous(trans = "log2") +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.5) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::ggtitle("Replicate Mean") +
    ggplot2::xlab("GRCh38 Genomic Coordinate") +
    ggplot2::ylab(
      glue::glue("Average {timepoint} SGE Counts / Plasmid Counts")
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = FALSE,
      span = 0.5,
      color = "blue"
    ) +
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

  print(p)
  invisible(p)
}


# Ridgeline and Jitter Plot Functions -------------------------------------------

#' Default consequence categories for distribution plots
#'
#' A character vector of VEP slim consequence types shown by default in
#' ridgeline and jitter plots. The order determines the display order.
#'
#' @export
default_consequences <- c(
  "synonymous",
  "UTR",
  "missense",
  "codon_deletion",
  "stop_gained",
  "frameshift",
  "stop_lost",
  "intron",
  "splice_donor",
  "splice_acceptor",
  "start_lost"
)

#' Plot ridgeline density distributions by group
#'
#' Creates a ridgeline (joy) plot showing the density distribution of a
#' numeric score across groups, with jittered data points overlaid.
#' Useful for comparing score distributions across variant consequence
#' categories.
#'
#' @param data A data frame containing the score and grouping columns.
#' @param score_col Character. Name of the column containing the numeric
#'   scores to plot on the x-axis.
#' @param group_col Character. Name of the column containing the grouping
#'   variable (e.g., consequence type) for the y-axis ridges.
#' @param consequences Character vector. Only include rows where the
#'   \code{group_col} value is in this vector. The factor levels of the
#'   grouping variable will be set to this order. Set to \code{NULL} to
#'   include all groups without filtering
#'   (default: [default_consequences]).
#' @param colours Named character vector of fill colours keyed by group
#'   values (default: [consequence_colours]).
#' @param xlab Character. X-axis label (default: value of \code{score_col}).
#' @param ylab Character. Y-axis label (default: "Density").
#' @param point_shape Shape for jittered points (default: "|").
#' @param point_size Size for jittered points (default: 3).
#' @param alpha Ridge fill transparency (default: 0.7).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_ridgeline(
#'   data = results,
#'   score_col = "adj_score_condition_Day7_vs_Day4",
#'   group_col = "slim_consequence"
#' )
#'
#' # With specific consequences in a custom order
#' plot_ridgeline(
#'   data = results,
#'   score_col = "adj_score_condition_Day7_vs_Day4",
#'   group_col = "slim_consequence",
#'   consequences = c("synonymous", "missense", "stop_gained"),
#'   xlab = "z-score (Day7 vs Day4)",
#'   ylab = "Consequence"
#' )
#' }
#'
#' @export
plot_ridgeline <- function(data,
                           score_col,
                           group_col,
                           consequences = default_consequences,
                           colours = consequence_colours,
                           xlab = score_col,
                           ylab = "Density",
                           point_shape = "|",
                           point_size = 3,
                           alpha = 0.7) {
  if (!is.null(consequences)) {
    data <- data[data[[group_col]] %in% consequences, , drop = FALSE]
    data[[group_col]] <- factor(data[[group_col]], levels = consequences)
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[score_col]],
      y = .data[[group_col]]
    )
  ) +
    ggridges::geom_density_ridges(
      ggplot2::aes(fill = .data[[group_col]]),
      jittered_points = TRUE,
      position = ggridges::position_points_jitter(
        width = 0.05, height = 0
      ),
      point_shape = point_shape,
      point_size = point_size,
      point_alpha = 1,
      alpha = alpha
    ) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")
}

#' Plot jittered strip chart by group
#'
#' Creates a jitter (strip) plot showing individual data points grouped by
#' a categorical variable, with optional median lines and FDR-based
#' transparency. Inspired by Figure 1E of Olvera-Leon et al. (2024) Cell.
#'
#' @param data A data frame containing the score, grouping, and
#'   optionally FDR columns.
#' @param score_col Character. Name of the column containing the numeric
#'   scores to plot on the y-axis.
#' @param group_col Character. Name of the column containing the grouping
#'   variable (e.g., consequence type) for the x-axis.
#' @param consequences Character vector. Only include rows where the
#'   \code{group_col} value is in this vector. The factor levels of the
#'   grouping variable will be set to this order. Set to \code{NULL} to
#'   include all groups without filtering
#'   (default: [default_consequences]).
#' @param fdr_col Character or NULL. Name of the column containing FDR
#'   values. When provided, points with FDR >= \code{fdr_threshold} are
#'   rendered semi-transparent.
#' @param fdr_threshold Numeric. FDR cutoff for full opacity
#'   (default: 0.01).
#' @param colours Named character vector of point colours keyed by group
#'   values (default: [consequence_colours]).
#' @param xlab Character. X-axis label (default: value of \code{group_col}).
#' @param ylab Character. Y-axis label (default: value of \code{score_col}).
#' @param show_median Logical. Whether to draw a horizontal median line
#'   per group (default: TRUE).
#' @param point_size Size of jittered points (default: 1).
#' @param jitter_width Width of horizontal jitter (default: 0.3).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_consequence_jitter(
#'   data = results,
#'   score_col = "adj_score_condition_Day14_vs_Day4",
#'   group_col = "slim_consequence"
#' )
#'
#' # With FDR-based transparency and selected consequences
#' plot_consequence_jitter(
#'   data = results,
#'   score_col = "adj_score_condition_Day14_vs_Day4",
#'   group_col = "slim_consequence",
#'   fdr_col = "FDR_condition_Day14_vs_Day4",
#'   consequences = c("synonymous", "missense", "stop_gained", "frameshift"),
#'   ylab = "z-score (D4-D14)"
#' )
#' }
#'
#' @export
plot_consequence_jitter <- function(data,
                                    score_col,
                                    group_col,
                                    consequences = default_consequences,
                                    fdr_col = NULL,
                                    fdr_threshold = 0.01,
                                    colours = consequence_colours,
                                    xlab = group_col,
                                    ylab = score_col,
                                    show_median = TRUE,
                                    point_size = 1,
                                    jitter_width = 0.3) {
  if (!is.null(consequences)) {
    data <- data[data[[group_col]] %in% consequences, , drop = FALSE]
    data[[group_col]] <- factor(data[[group_col]], levels = consequences)
  }

  # Build alpha mapping based on FDR if provided
  if (!is.null(fdr_col)) {
    data$.alpha <- ifelse(data[[fdr_col]] < fdr_threshold, 1, 0.2)
  } else {
    data$.alpha <- 1
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[group_col]],
      y = .data[[score_col]],
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_jitter(
      ggplot2::aes(alpha = .data$.alpha),
      width = jitter_width,
      size = point_size,
      shape = 21,
      colour = "black",
      stroke = 0.3
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (show_median) {
    p <- p +
      ggplot2::stat_summary(
        fun = stats::median,
        geom = "crossbar",
        width = 0.5,
        linewidth = 0.4,
        colour = "black"
      )
  }

  p
}

#' Plot functional scores along genomic coordinates
#'
#' Creates a waterfall-style scatter plot showing variant scores along
#' genomic coordinates, colored by consequence type and shaped by
#' functional classification. Supports optional faceting by a grouping
#' variable such as exon. Column names are fully configurable.
#'
#' @param data A data frame containing the columns referenced by the
#'   other parameters.
#' @param score_col Character. Name of the column containing the numeric
#'   score to plot on the y-axis.
#' @param pos_col Character. Name of the column containing the genomic
#'   coordinate for the x-axis (default: "vcf_pos").
#' @param consequence_col Character. Name of the column containing
#'   consequence annotations for colour mapping
#'   (default: "slim_consequence").
#' @param class_col Character or NULL. Name of the column containing
#'   functional classifications for point shapes. If NULL, no shape
#'   mapping is applied (default: "functional_classification").
#' @param fdr_col Character or NULL. Name of the column containing FDR
#'   values for alpha transparency. If NULL, all points are fully opaque
#'   (default: "FDR").
#' @param fdr_threshold Numeric. FDR cutoff for alpha transparency
#'   (default: 0.01).
#' @param facet_col Character or NULL. Name of the column to facet by
#'   (e.g., "Exons"). If NULL, no faceting is applied (default: NULL).
#' @param facet_nrow Integer. Number of rows in facet layout
#'   (default: 1).
#' @param title Character or NULL. Plot title. If NULL, no title is
#'   shown (default: NULL).
#' @param xlab Character. X-axis label (default: "GRCh38 coordinate").
#' @param ylab Character. Y-axis label (default: value of
#'   \code{score_col}).
#' @param colours Named character vector of colours for consequence
#'   values (default: [consequence_colours]).
#' @param shapes Named numeric vector of shapes for classification
#'   values (default: [classification_shapes]).
#' @param point_size Numeric. Size of points (default: 1.5).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_gene_level_waterfall(
#'   data = combined_dataset,
#'   score_col = "combined_lfc",
#'   class_col = "functional_classification",
#'   fdr_col = "FDR",
#'   facet_col = "Exons",
#'   title = "Functional score Plot: NUDCD3 - Weighted LFC",
#'   ylab = "Functional score (Weighted LFC)"
#' )
#' }
#'
#' @export
plot_gene_level_waterfall <- function(data,
                                      score_col,
                                      pos_col = "vcf_pos",
                                      consequence_col = "slim_consequence",
                                      class_col = "functional_classification",
                                      fdr_col = "FDR",
                                      fdr_threshold = 0.01,
                                      facet_col = NULL,
                                      facet_nrow = 1,
                                      title = NULL,
                                      xlab = "GRCh38 coordinate",
                                      ylab = score_col,
                                      colours = consequence_colours,
                                      shapes = classification_shapes,
                                      point_size = 1.5) {
  mapping <- ggplot2::aes(
    x = .data[[pos_col]],
    y = .data[[score_col]],
    colour = .data[[consequence_col]],
    fill = .data[[consequence_col]]
  )

  if (!is.null(class_col)) {
    mapping$shape <- ggplot2::aes(shape = .data[[class_col]])$shape
  }

  if (!is.null(fdr_col)) {
    point_mapping <- ggplot2::aes(
      alpha = .data[[fdr_col]] <= fdr_threshold
    )
  } else {
    point_mapping <- NULL
  }

  p <- ggplot2::ggplot(data, mapping) +
    ggplot2::geom_point(point_mapping, size = point_size) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
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

  if (!is.null(class_col)) {
    p <- p + ggplot2::scale_shape_manual(values = shapes)
  }

  if (!is.null(fdr_col)) {
    p <- p + ggplot2::scale_alpha_discrete(
      labels = c(
        paste0("FDR > ", fdr_threshold),
        paste0("FDR <= ", fdr_threshold)
      )
    )
  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(
      facet_col,
      scales = "free_x",
      nrow = facet_nrow
    )
  }

  p
}

# Position Effect Plot Functions -----------------------------------------------

#' Save position effect ratio plots to files
#'
#' Saves multiple position effect ratio plots to PNG files.
#'
#' @param pep_list A named list of data frames from
#'   [calculate_position_effect()].
#' @param output_dir Directory to save plots
#'   (default: "results/plasmid_ratio").
#' @param timepoint Character string specifying the timepoint
#'   (default: "Day4").
#' @param width Plot width in pixels (default: 900).
#' @param height Plot height in pixels (default: 600).
#'
#' @return Invisibly returns NULL. Side effect: writes PNG files.
#'
#' @examples
#' \dontrun{
#' plot_position_effects(pep_list, output_dir = "results/plasmid_ratio")
#' }
#'
#' @seealso [position_effect_plot()], [calculate_position_effect()]
#' @export
plot_position_effects <- function(pep_list,
                                  output_dir = "results/plasmid_ratio",
                                  timepoint = "Day4",
                                  width = 900,
                                  height = 600) {
  .save_plots_to_dir(
    pep_list, output_dir, "_ratio", width, height,
    function(x) print(position_effect_plot(x, timepoint = timepoint))
  )
}
