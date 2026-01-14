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
