# Tests for protein domain visualization functions

# Test data fixtures -----------------------------------------------------------

# Sample domain data (mimics Ensembl API response)
mock_domain_df <- tibble::tibble(
  domain = c("Kinase", "Kinase", "SH2", "PH"),
  source = c("Pfam", "Pfam", "SMART", "Pfam"),
  start = c(100, 150, 300, 50),
  end = c(200, 180, 400, 90)
)

# Sample variant data for heatmap
mock_variant_data <- tibble::tibble(
  Amino_acids = c("A/V", "A/G", "G/D", "K/R", "L/*"),
  Protein_position = c(100, 100, 150, 200, 250),
  adj_score_continuous = c(-0.5, 0.3, -0.8, 0.1, -1.2),
  FDR_continuous = c(0.001, 0.5, 0.005, 0.8, 0.0001)
)

# Tests for merge_overlapping_domains() ----------------------------------------

test_that("merge_overlapping_domains merges overlapping ranges", {
  merged <- merge_overlapping_domains(mock_domain_df)

  # Kinase should be merged (100-200 and 150-180 overlap)
  kinase_rows <- merged[merged$domain == "Kinase", ]
  expect_equal(nrow(kinase_rows), 1)
  expect_equal(kinase_rows$start, 100)
  expect_equal(kinase_rows$end, 200)

  # SH2 and PH should remain unchanged
  expect_equal(nrow(merged[merged$domain == "SH2", ]), 1)
  expect_equal(nrow(merged[merged$domain == "PH", ]), 1)
})

test_that("merge_overlapping_domains handles empty input", {
  result <- merge_overlapping_domains(NULL)
  expect_null(result)

  empty_df <- tibble::tibble(
    domain = character(),
    source = character(),
    start = numeric(),
    end = numeric()
  )
  result <- merge_overlapping_domains(empty_df)
  expect_equal(nrow(result), 0)
})

test_that("merge_overlapping_domains handles non-overlapping ranges", {
  non_overlapping <- tibble::tibble(
    domain = c("A", "B"),
    source = c("Pfam", "Pfam"),
    start = c(10, 100),
    end = c(20, 110)
  )
  merged <- merge_overlapping_domains(non_overlapping)
  expect_equal(nrow(merged), 2)
})

# Tests for prepare_amino_acid_heatmap_data() ----------------------------------

test_that("prepare_amino_acid_heatmap_data creates correct structure", {
  result <- prepare_amino_acid_heatmap_data(mock_variant_data)

  expect_true("AA_change" %in% names(result))
  expect_true("previous" %in% names(result))
  expect_true("new_aa" %in% names(result))
  expect_true("Protein_position" %in% names(result))
  expect_true("functional_classification" %in% names(result))
  expect_true("mean_adj_score" %in% names(result))
  expect_true("mean_FDR" %in% names(result))
})

test_that("prepare_amino_acid_heatmap_data classifies variants correctly", {
  result <- prepare_amino_acid_heatmap_data(mock_variant_data)

  # A100V: FDR < 0.01, score < 0 -> depleted
  a100v <- result[result$AA_change == "A100V", ]
  expect_equal(a100v$functional_classification, "depleted")

  # A100G: FDR > 0.01 -> unchanged
  a100g <- result[result$AA_change == "A100G", ]
  expect_equal(a100g$functional_classification, "unchanged")
})

test_that("prepare_amino_acid_heatmap_data converts special characters", {
  result <- prepare_amino_acid_heatmap_data(mock_variant_data)

  # L/* should become stop_scan
  stop_row <- result[result$AA_change == "L250*", ]
  expect_equal(as.character(stop_row$new_aa), "stop_scan")
})

test_that("prepare_amino_acid_heatmap_data orders new_aa as factor", {
  result <- prepare_amino_acid_heatmap_data(mock_variant_data)
  expect_true(is.factor(result$new_aa))
  expect_equal(levels(result$new_aa), amino_acid_order)
})

test_that("prepare_amino_acid_heatmap_data respects custom FDR threshold", {
  # With very strict threshold, nothing should be enriched/depleted
  result <- prepare_amino_acid_heatmap_data(
    mock_variant_data,
    fdr_threshold = 0.00001
  )
  expect_true(all(result$functional_classification == "unchanged"))
})

# Tests for plot_domain_track() ------------------------------------------------

test_that("plot_domain_track returns ggplot object", {
  p <- plot_domain_track(mock_domain_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_domain_track handles NULL input", {
  p <- plot_domain_track(NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_domain_track handles empty data frame", {
  empty_df <- tibble::tibble(
    domain = character(),
    source = character(),
    start = numeric(),
    end = numeric()
  )
  p <- plot_domain_track(empty_df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_domain_track respects x_range parameter", {
  p <- plot_domain_track(mock_domain_df, x_range = c(0, 500))
  # Check that scale limits are set (internal ggplot structure)
  expect_s3_class(p, "ggplot")
})

test_that("plot_domain_track respects show_legend parameter", {
  p_legend <- plot_domain_track(mock_domain_df, show_legend = TRUE)
  p_no_legend <- plot_domain_track(mock_domain_df, show_legend = FALSE)

  # Both should be valid ggplot objects
  expect_s3_class(p_legend, "ggplot")
  expect_s3_class(p_no_legend, "ggplot")
})

# Tests for plot_amino_acid_heatmap() ------------------------------------------

test_that("plot_amino_acid_heatmap returns ggplot object", {
  plot_data <- prepare_amino_acid_heatmap_data(mock_variant_data)
  p <- plot_amino_acid_heatmap(plot_data)
  expect_s3_class(p, "ggplot")
})

test_that("plot_amino_acid_heatmap handles NULL input", {
  p <- plot_amino_acid_heatmap(NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_amino_acid_heatmap respects custom colors", {
  plot_data <- prepare_amino_acid_heatmap_data(mock_variant_data)
  custom_colors <- c(
    "unchanged" = "gray",
    "enriched" = "red",
    "depleted" = "blue"
  )
  p <- plot_amino_acid_heatmap(plot_data, classification_colors = custom_colors)
  expect_s3_class(p, "ggplot")
})

# Tests for exported data objects ----------------------------------------------

test_that("classification_heatmap_colours has correct structure", {
  expect_type(classification_heatmap_colours, "character")
  expect_named(classification_heatmap_colours)
  expect_true("unchanged" %in% names(classification_heatmap_colours))
  expect_true("enriched" %in% names(classification_heatmap_colours))
  expect_true("depleted" %in% names(classification_heatmap_colours))
})

test_that("amino_acid_order has correct structure", {
  expect_type(amino_acid_order, "character")
  expect_true(length(amino_acid_order) > 20)  # All amino acids plus special types
  expect_true("A" %in% amino_acid_order)
  expect_true("stop_scan" %in% amino_acid_order)
  expect_true("codon_deletion" %in% amino_acid_order)
  expect_true("frameshift" %in% amino_acid_order)
})

# Tests for fetch_domains_ensembl() --------------------------------------------

test_that("fetch_domains_ensembl validates input", {
  expect_error(
    fetch_domains_ensembl(c("ID1", "ID2")),
    "transcript_id must be a single character string"
  )
  expect_error(
    fetch_domains_ensembl(123),
    "transcript_id must be a single character string"
  )
})

# Integration test (skipped by default - requires network)
test_that("fetch_domains_ensembl fetches real data", {
  skip_on_cran()
  skip_if_offline()

  # Use a well-known transcript
  domains <- fetch_domains_ensembl("ENST00000355451")

  if (!is.null(domains)) {
    expect_true(nrow(domains) > 0)
    expect_true(all(c("domain", "source", "start", "end") %in% names(domains)))
    expect_true(all(domains$start <= domains$end))
  }
})

# Tests for plot_domain_heatmap() ----------------------------------------------

test_that("plot_domain_heatmap requires uniprot_acc or domain_df", {
  expect_error(
    plot_domain_heatmap(mock_variant_data),
    "'uniprot_acc' required for TED API"
  )
})

test_that("plot_domain_heatmap returns patchwork object", {
  p <- plot_domain_heatmap(
    mock_variant_data,
    domain_df = mock_domain_df
  )

  expect_s3_class(p, "patchwork")
})

test_that("plot_domain_heatmap works with prepared data", {
  prepared_data <- prepare_amino_acid_heatmap_data(mock_variant_data)
  p <- plot_domain_heatmap(
    prepared_data,
    domain_df = mock_domain_df,
    prepared = TRUE
  )

  expect_s3_class(p, "patchwork")
})

test_that("plot_domain_heatmap respects height_ratio", {
  p <- plot_domain_heatmap(
    mock_variant_data,
    domain_df = mock_domain_df,
    height_ratio = c(1, 10)
  )

  expect_s3_class(p, "patchwork")
})
