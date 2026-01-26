# Test SGEResults S7 class

test_that("SGEResults can be created with valid inputs", {
  # Create mock data
  results_df <- data.frame(
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    baseMean = c(100, 200, 300),
    log2FoldChange = c(-1, 0, 1)
  )

  contrast_df <- data.frame(
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    log2FoldChange_Day7_vs_Day4 = c(-1, 0, 1),
    pvalue_Day7_vs_Day4 = c(0.01, 0.5, 0.02)
  )

  # Create a mock DESeqTransform-like object
  # For testing, we need to skip this test if DESeq2 isn't available
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  # Create minimal DESeq2 objects for testing
  count_data <- matrix(
    c(100, 200, 300, 110, 190, 310),
    nrow = 3,
    dimnames = list(c("ATCG", "GCTA", "TTAA"), c("sample1", "sample2"))
  )
  col_data <- data.frame(
    condition = factor(c("Day4", "Day7")),
    row.names = c("sample1", "sample2")
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  rld <- DESeq2::rlog(dds)

  # Create SGEResults object
  sge_obj <- SGEResults(
    results = results_df,
    rlog = rld,
    contrast_summary = contrast_df,
    metadata = list(alpha = 0.05)
  )

  expect_s3_class(sge_obj, "SGEResults")
  expect_equal(nrow(sge_obj@results), 3)
  expect_equal(nrow(sge_obj@contrast_summary), 3)
  expect_equal(sge_obj@metadata$alpha, 0.05)
})

test_that("SGEResults validator rejects missing SEQUENCE in results", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  # Results without SEQUENCE column
  results_df <- data.frame(
    baseMean = c(100, 200, 300),
    log2FoldChange = c(-1, 0, 1)
  )

  contrast_df <- data.frame(
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    log2FoldChange_Day7_vs_Day4 = c(-1, 0, 1)
  )

  # Create mock rlog
  count_data <- matrix(
    c(100, 200, 300, 110, 190, 310),
    nrow = 3,
    dimnames = list(c("ATCG", "GCTA", "TTAA"), c("sample1", "sample2"))
  )
  col_data <- data.frame(
    condition = factor(c("Day4", "Day7")),
    row.names = c("sample1", "sample2")
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  rld <- DESeq2::rlog(dds)

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
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    baseMean = c(100, 200, 300)
  )

  # contrast_summary without SEQUENCE
  contrast_df <- data.frame(
    log2FoldChange_Day7_vs_Day4 = c(-1, 0, 1)
  )

  count_data <- matrix(
    c(100, 200, 300, 110, 190, 310),
    nrow = 3,
    dimnames = list(c("ATCG", "GCTA", "TTAA"), c("sample1", "sample2"))
  )
  col_data <- data.frame(
    condition = factor(c("Day4", "Day7")),
    row.names = c("sample1", "sample2")
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  rld <- DESeq2::rlog(dds)

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
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    baseMean = c(100, 200, 300)
  )

  # contrast_summary with different number of rows
  contrast_df <- data.frame(
    SEQUENCE = c("ATCG", "GCTA"),
    log2FoldChange_Day7_vs_Day4 = c(-1, 0)
  )

  count_data <- matrix(
    c(100, 200, 300, 110, 190, 310),
    nrow = 3,
    dimnames = list(c("ATCG", "GCTA", "TTAA"), c("sample1", "sample2"))
  )
  col_data <- data.frame(
    condition = factor(c("Day4", "Day7")),
    row.names = c("sample1", "sample2")
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  rld <- DESeq2::rlog(dds)

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
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    baseMean = c(100, 200, 300)
  )

  contrast_df <- data.frame(
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    log2FoldChange_Day7_vs_Day4 = c(-1, 0, 1)
  )

  # Pass a plain matrix instead of DESeqTransform
  fake_rlog <- matrix(1:6, nrow = 3)

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
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    baseMean = c(100, 200, 300),
    log2FoldChange = c(-1, 0, 1)
  )

  contrast_df <- data.frame(
    SEQUENCE = c("ATCG", "GCTA", "TTAA"),
    log2FoldChange_Day7_vs_Day4 = c(-1, 0, 1),
    pvalue_Day7_vs_Day4 = c(0.01, 0.5, 0.02)
  )

  count_data <- matrix(
    c(100, 200, 300, 110, 190, 310),
    nrow = 3,
    dimnames = list(c("ATCG", "GCTA", "TTAA"), c("sample1", "sample2"))
  )
  col_data <- data.frame(
    condition = factor(c("Day4", "Day7")),
    row.names = c("sample1", "sample2")
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ condition
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  rld <- DESeq2::rlog(dds)

  sge_obj <- SGEResults(
    results = results_df,
    rlog = rld,
    contrast_summary = contrast_df,
    metadata = list(alpha = 0.05)
  )

  # Print should not error and should return invisibly
  expect_output(print(sge_obj), "SGE Analysis Results")
  expect_invisible(print(sge_obj))
})
