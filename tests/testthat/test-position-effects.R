# Test per-timepoint position effect correction

# =============================================================================
# calculate_all_position_effects tests
# =============================================================================

test_that("calculate_all_position_effects returns named list with expected timepoints", {
  # Minimal fake data matching the structure expected by calculate_position_effect
  timepoints <- c("Day4", "Day7")

  dataframe <- data.frame(
    SEQUENCE = rep(paste0("SEQ", 1:3), each = 3),
    NAME = rep(paste0("NAME", 1:3), each = 3),
    COUNT = c(100, 200, 150, 120, 180, 160, 90, 210, 170),
    condition = rep(c("Day0", "Day4", "Day7"), 3),
    targeton_id = "TGT1"
  )

  annotation <- data.frame(
    Seq = paste0("SEQ", 1:3),
    Targeton_ID = "TGT1",
    vcf_pos = c(100, 200, 300)
  )

  result <- calculate_all_position_effects(dataframe, annotation, timepoints)

  expect_type(result, "list")
  expect_named(result, timepoints)
  expect_length(result, 2)

  # Each element should be a data frame with ratio column
  expect_s3_class(result[["Day4"]], "data.frame")
  expect_s3_class(result[["Day7"]], "data.frame")
  expect_true("ratio" %in% names(result[["Day4"]]))
  expect_true("ratio" %in% names(result[["Day7"]]))
})


# =============================================================================
# .build_positional_norm_matrix tests
# =============================================================================

test_that("per-timepoint normalization matrix has different pos_effects per condition", {
  # Create a named list of pos_effect data frames with different values
  pos_effect_list <- list(
    Day4 = data.frame(
      SEQUENCE = paste0("SEQ", 1:3),
      pos_effect = c(0.5, 1.0, -0.3)
    ),
    Day7 = data.frame(
      SEQUENCE = paste0("SEQ", 1:3),
      pos_effect = c(0.2, -0.5, 0.8)
    )
  )

  count_matrix <- matrix(
    c(100, 200, 150, 120, 180, 160, 90, 210, 170, 130, 190, 140),
    nrow = 3, ncol = 4,
    dimnames = list(
      paste0("SEQ", 1:3),
      c("S1", "S2", "S3", "S4")
    )
  )

  size_factors <- c(S1 = 1.0, S2 = 1.1, S3 = 0.9, S4 = 1.2)

  condition_data <- data.frame(
    condition = factor(c("Day4", "Day4", "Day7", "Day7")),
    row.names = c("S1", "S2", "S3", "S4")
  )

  norm_matrix <- sgeasy:::.build_positional_norm_matrix(
    pos_effect_data = pos_effect_list,
    count_matrix = count_matrix,
    size_factors = size_factors,
    condition_data = condition_data,
    reference_level = "Day0"
  )

  expect_equal(dim(norm_matrix), dim(count_matrix))

  # Day4 pos_effects centered: c(0.5, 1.0, -0.3) - median(0.5) = c(0, 0.5, -0.8)
  day4_centered <- c(0.5, 1.0, -0.3) - median(c(0.5, 1.0, -0.3))
  # Day7 pos_effects centered: c(0.2, -0.5, 0.8) - median(0.2) = c(0, -0.7, 0.6)
  day7_centered <- c(0.2, -0.5, 0.8) - median(c(0.2, -0.5, 0.8))

  # Day4 samples (S1, S2) should use centered Day4 pos_effects
  expect_equal(norm_matrix[1, "S1"], size_factors["S1"] * 2^day4_centered[1], ignore_attr = TRUE)
  expect_equal(norm_matrix[2, "S2"], size_factors["S2"] * 2^day4_centered[2], ignore_attr = TRUE)

  # Day7 samples (S3, S4) should use centered Day7 pos_effects
  expect_equal(norm_matrix[1, "S3"], size_factors["S3"] * 2^day7_centered[1], ignore_attr = TRUE)
  expect_equal(norm_matrix[3, "S4"], size_factors["S4"] * 2^day7_centered[3], ignore_attr = TRUE)

  # Day4 and Day7 columns should have DIFFERENT pos_effects for same gene
  expect_false(
    identical(norm_matrix[, "S1"] / size_factors["S1"],
              norm_matrix[, "S3"] / size_factors["S3"])
  )
})


test_that("reference columns have pos_effect = 0 in normalization matrix", {
  # Single data frame input
  pos_effect_df <- data.frame(
    SEQUENCE = paste0("SEQ", 1:3),
    pos_effect = c(0.5, 1.0, -0.3)
  )

  count_matrix <- matrix(
    c(100, 200, 150, 120, 180, 160, 90, 210, 170),
    nrow = 3, ncol = 3,
    dimnames = list(
      paste0("SEQ", 1:3),
      c("S1", "S2", "S3")
    )
  )

  size_factors <- c(S1 = 1.0, S2 = 1.1, S3 = 0.9)

  condition_data <- data.frame(
    condition = factor(c("Day0", "Day4", "Day7")),
    row.names = c("S1", "S2", "S3")
  )

  norm_matrix <- sgeasy:::.build_positional_norm_matrix(
    pos_effect_data = pos_effect_df,
    count_matrix = count_matrix,
    size_factors = size_factors,
    condition_data = condition_data,
    reference_level = "Day0"
  )

  # Reference column (S1, Day0) should have norm = size_factor only
  expect_equal(norm_matrix[, "S1"], rep(size_factors["S1"], 3), ignore_attr = TRUE)

  # Non-reference columns should have centered positional offsets applied
  # c(0.5, 1.0, -0.3) - median(0.5) = c(0, 0.5, -0.8)
  centered <- c(0.5, 1.0, -0.3) - median(c(0.5, 1.0, -0.3))
  expect_equal(norm_matrix[1, "S2"], size_factors["S2"] * 2^centered[1], ignore_attr = TRUE)
  expect_equal(norm_matrix[2, "S3"], size_factors["S3"] * 2^centered[2], ignore_attr = TRUE)
})


test_that("reference columns have pos_effect = 0 with list input", {
  pos_effect_list <- list(
    Day4 = data.frame(
      SEQUENCE = paste0("SEQ", 1:3),
      pos_effect = c(0.5, 1.0, -0.3)
    )
  )

  count_matrix <- matrix(
    c(100, 200, 150, 120, 180, 160),
    nrow = 3, ncol = 2,
    dimnames = list(
      paste0("SEQ", 1:3),
      c("S1", "S2")
    )
  )

  size_factors <- c(S1 = 1.0, S2 = 1.1)

  condition_data <- data.frame(
    condition = factor(c("Day0", "Day4")),
    row.names = c("S1", "S2")
  )

  norm_matrix <- sgeasy:::.build_positional_norm_matrix(
    pos_effect_data = pos_effect_list,
    count_matrix = count_matrix,
    size_factors = size_factors,
    condition_data = condition_data,
    reference_level = "Day0"
  )

  # Reference column should be pure size factors
  expect_equal(norm_matrix[, "S1"], rep(size_factors["S1"], 3), ignore_attr = TRUE)

  # Non-reference should have centered pos_effect applied
  # c(0.5, 1.0, -0.3) - median(0.5) = c(0, 0.5, -0.8)
  centered <- c(0.5, 1.0, -0.3) - median(c(0.5, 1.0, -0.3))
  expect_equal(norm_matrix[1, "S2"], size_factors["S2"] * 2^centered[1], ignore_attr = TRUE)
})


test_that("variants not in pos_effect_data get offset of 0", {
  pos_effect_df <- data.frame(
    SEQUENCE = c("SEQ1"),
    pos_effect = c(0.5)
  )

  count_matrix <- matrix(
    c(100, 200, 120, 180),
    nrow = 2, ncol = 2,
    dimnames = list(
      c("SEQ1", "SEQ2"),  # SEQ2 not in pos_effect_data
      c("S1", "S2")
    )
  )

  size_factors <- c(S1 = 1.0, S2 = 1.1)

  condition_data <- data.frame(
    condition = factor(c("Day0", "Day4")),
    row.names = c("S1", "S2")
  )

  norm_matrix <- sgeasy:::.build_positional_norm_matrix(
    pos_effect_data = pos_effect_df,
    count_matrix = count_matrix,
    size_factors = size_factors,
    condition_data = condition_data,
    reference_level = "Day0"
  )

  # After centering: c(0.5, 0) - median(c(0.5, 0)) = c(0.25, -0.25)
  centered <- c(0.5, 0) - median(c(0.5, 0))

  # SEQ2 (missing from pos_effect_data) gets centered offset
  expect_equal(norm_matrix["SEQ2", "S2"], size_factors["S2"] * 2^centered[2], ignore_attr = TRUE)

  # SEQ1 should have its centered offset
  expect_equal(norm_matrix["SEQ1", "S2"], size_factors["S2"] * 2^centered[1], ignore_attr = TRUE)
})


# =============================================================================
# Per-targeton empirical null tests
# =============================================================================

test_that("per-targeton medians differ between targetons", {
  data <- data.frame(
    Targeton_ID = rep(c("TGT1", "TGT2"), each = 10),
    Consequence = rep(c("synonymous_variant", "missense_variant"), 10),
    log2FoldChange_condition_Day7_vs_Day4 = c(
      # TGT1: neutral variants have higher LFCs
      rnorm(5, mean = 2.0, sd = 0.1), rnorm(5, mean = 3.0, sd = 0.5),
      # TGT2: neutral variants have lower LFCs
      rnorm(5, mean = -1.0, sd = 0.1), rnorm(5, mean = 0.5, sd = 0.5)
    ),
    lfcSE_condition_Day7_vs_Day4 = rep(0.3, 20)
  )

  result <- sgeasy::calculate_median_scores(data, per_targeton = TRUE)

  median_col <- "median_log2FoldChange_condition_Day7_vs_Day4"
  expect_true(median_col %in% names(result))

  tgt1_median <- unique(result[[median_col]][result$Targeton_ID == "TGT1"])
  tgt2_median <- unique(result[[median_col]][result$Targeton_ID == "TGT2"])

  # Medians should be different between targetons
  expect_length(tgt1_median, 1)
  expect_length(tgt2_median, 1)
  expect_true(abs(tgt1_median - tgt2_median) > 1)
})


test_that("empirical null SD is used when per_targeton = TRUE", {
  set.seed(42)
  n_per_tgt <- 50

  data <- data.frame(
    Targeton_ID = rep(c("TGT1", "TGT2"), each = n_per_tgt),
    Consequence = rep(c(rep("synonymous_variant", 30),
                        rep("missense_variant", 20)), 2),
    SEQUENCE = paste0("SEQ", 1:(2 * n_per_tgt)),
    log2FoldChange_condition_Day7_vs_Day4 = c(
      rnorm(n_per_tgt, mean = 0.5, sd = 0.8),
      rnorm(n_per_tgt, mean = -0.3, sd = 0.4)
    ),
    lfcSE_condition_Day7_vs_Day4 = rep(0.1, 2 * n_per_tgt)
  )

  # With per_targeton = TRUE: should use empirical SD
  result_pt <- sgeasy::recalculate_screen_statistics(
    data, per_targeton = TRUE
  )

  # With per_targeton = FALSE: should use lfcSE
  result_global <- sgeasy::recalculate_screen_statistics(
    data, per_targeton = FALSE
  )

  # The adj_scores should be different
  expect_false(
    identical(
      result_pt$adj_score_condition_Day7_vs_Day4,
      result_global$adj_score_condition_Day7_vs_Day4
    )
  )

  # With per_targeton, synonymous adj_score SD should be closer to 1
  syn_scores_pt <- result_pt$adj_score_condition_Day7_vs_Day4[
    result_pt$Consequence == "synonymous_variant"
  ]
  syn_scores_global <- result_global$adj_score_condition_Day7_vs_Day4[
    result_global$Consequence == "synonymous_variant"
  ]

  # Empirical null should produce scores with SD much closer to 1
  # lfcSE=0.1 is deliberately too small, so global scores are inflated
  expect_true(abs(sd(syn_scores_pt, na.rm = TRUE) - 1) <
                abs(sd(syn_scores_global, na.rm = TRUE) - 1))
})


test_that(".compute_null_sd returns per-targeton SD", {
  set.seed(123)
  data <- data.frame(
    Targeton_ID = rep(c("TGT1", "TGT2"), each = 20),
    Consequence = rep(c(rep("synonymous_variant", 15),
                        rep("missense_variant", 5)), 2),
    adj_lfc_condition_Day7_vs_Day4 = c(
      rnorm(20, sd = 0.5),
      rnorm(20, sd = 1.5)
    )
  )

  result <- sgeasy:::.compute_null_sd(
    data, "condition_Day7_vs_Day4",
    neutral_variants = c("synonymous_variant")
  )

  expect_true("Targeton_ID" %in% names(result))
  expect_true("null_sd_condition_Day7_vs_Day4" %in% names(result))
  expect_equal(nrow(result), 2)

  # TGT2 should have larger null SD than TGT1 (sd=1.5 vs sd=0.5)
  tgt1_sd <- result$null_sd_condition_Day7_vs_Day4[result$Targeton_ID == "TGT1"]
  tgt2_sd <- result$null_sd_condition_Day7_vs_Day4[result$Targeton_ID == "TGT2"]
  expect_true(tgt2_sd > tgt1_sd)
})
