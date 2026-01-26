# Test SGEResults S7 class

# Helper function to create a mock rlog object for testing
# DESeq2's rlog() needs realistic data, so we create minimal valid objects
create_mock_rlog <- function(n_rows = 10, n_samples = 4) {
  # Create count data with enough variation for DESeq2
set.seed(42)
  count_data <- matrix(
    rpois(n_rows * n_samples, lambda = 100),
    nrow = n_rows,
    dimnames = list(
      paste0("SEQ", seq_len(n_rows)),
      paste0("sample", seq_len(n_samples))
    )
  )

  col_data <- data.frame(
    condition = factor(rep(c("Day4", "Day7"), each = n_samples / 2)),
    row.names = paste0("sample", seq_len(n_samples))
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)

  # Use blind=FALSE and fitType="mean" for small datasets
  DESeq2::rlog(dds, blind = FALSE, fitType = "mean")
}

test_that("SGEResults can be created with valid inputs", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  # Create mock data
  results_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    baseMean = runif(10, 50, 200),
    log2FoldChange = rnorm(10)
  )

  contrast_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    log2FoldChange_Day7_vs_Day4 = rnorm(10),
    pvalue_Day7_vs_Day4 = runif(10, 0, 0.1)
  )

  rld <- create_mock_rlog()

  # Create SGEResults object
  sge_obj <- SGEResults(
    results = results_df,
    rlog = rld,
    contrast_summary = contrast_df,
    metadata = list(alpha = 0.05)
  )

  expect_true(inherits(sge_obj, "sgeasy::SGEResults"))
  expect_equal(nrow(sge_obj@results), 10)
  expect_equal(nrow(sge_obj@contrast_summary), 10)
  expect_equal(sge_obj@metadata$alpha, 0.05)
})

test_that("SGEResults validator rejects missing SEQUENCE in results", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  # Results without SEQUENCE column
  results_df <- data.frame(
    baseMean = runif(10, 50, 200),
    log2FoldChange = rnorm(10)
  )

  contrast_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    log2FoldChange_Day7_vs_Day4 = rnorm(10)
  )

  rld <- create_mock_rlog()

  expect_error(
    SGEResults(
      results = results_df,
      rlog = rld,
      contrast_summary = contrast_df
    ),
    "results must contain a SEQUENCE column"
  )
})

test_that("SGEResults validator rejects missing SEQUENCE in contrast_summary", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  results_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    baseMean = runif(10, 50, 200)
  )

  # contrast_summary without SEQUENCE
  contrast_df <- data.frame(
    log2FoldChange_Day7_vs_Day4 = rnorm(10)
  )

  rld <- create_mock_rlog()

  expect_error(
    SGEResults(
      results = results_df,
      rlog = rld,
      contrast_summary = contrast_df
    ),
    "contrast_summary must contain a SEQUENCE column"
  )
})

test_that("SGEResults validator rejects mismatched row counts", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  results_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    baseMean = runif(10, 50, 200)
  )

  # contrast_summary with different number of rows
  contrast_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:5),
    log2FoldChange_Day7_vs_Day4 = rnorm(5)
  )

  rld <- create_mock_rlog()

  expect_error(
    SGEResults(
      results = results_df,
      rlog = rld,
      contrast_summary = contrast_df
    ),
    "results and contrast_summary must have same number of rows"
  )
})

test_that("SGEResults validator rejects invalid rlog type", {
  results_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    baseMean = runif(10, 50, 200)
  )

  contrast_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    log2FoldChange_Day7_vs_Day4 = rnorm(10)
  )

  # Pass a plain matrix instead of DESeqTransform
  fake_rlog <- matrix(1:40, nrow = 10)

  expect_error(
    SGEResults(
      results = results_df,
      rlog = fake_rlog,
      contrast_summary = contrast_df
    ),
    "rlog must be a DESeqTransform object"
  )
})

test_that("SGEResults print method runs without error", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("cli")

  results_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    baseMean = runif(10, 50, 200),
    log2FoldChange = rnorm(10)
  )

  contrast_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:10),
    log2FoldChange_Day7_vs_Day4 = rnorm(10),
    pvalue_Day7_vs_Day4 = runif(10, 0, 0.1)
  )

  rld <- create_mock_rlog()

  sge_obj <- SGEResults(
    results = results_df,
    rlog = rld,
    contrast_summary = contrast_df,
    metadata = list(alpha = 0.05)
  )

  # Print should not error and should return invisibly
  # cli output goes to stderr/message, so we capture with cli's test helper
  expect_no_error(print(sge_obj))
  expect_invisible(print(sge_obj))
})
