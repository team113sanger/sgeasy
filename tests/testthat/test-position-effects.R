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

  # Day4 samples (S1, S2) should use Day4 pos_effects
  expect_equal(norm_matrix[1, "S1"], size_factors["S1"] * 2^0.5, ignore_attr = TRUE)
  expect_equal(norm_matrix[2, "S2"], size_factors["S2"] * 2^1.0, ignore_attr = TRUE)

  # Day7 samples (S3, S4) should use Day7 pos_effects
  expect_equal(norm_matrix[1, "S3"], size_factors["S3"] * 2^0.2, ignore_attr = TRUE)
  expect_equal(norm_matrix[3, "S4"], size_factors["S4"] * 2^0.8, ignore_attr = TRUE)

  # Day4 and Day7 columns should have DIFFERENT pos_effects for same gene
  # (unless pos_effects happen to be equal)
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

  # Non-reference columns should have positional offsets applied
  expect_equal(norm_matrix[1, "S2"], size_factors["S2"] * 2^0.5, ignore_attr = TRUE)
  expect_equal(norm_matrix[2, "S3"], size_factors["S3"] * 2^1.0, ignore_attr = TRUE)
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

  # Non-reference should have pos_effect applied
  expect_equal(norm_matrix[1, "S2"], size_factors["S2"] * 2^0.5, ignore_attr = TRUE)
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

  # SEQ2 should have no positional offset (2^0 = 1)
  expect_equal(norm_matrix["SEQ2", "S2"], size_factors["S2"], ignore_attr = TRUE)

  # SEQ1 should have the offset

  expect_equal(norm_matrix["SEQ1", "S2"], size_factors["S2"] * 2^0.5, ignore_attr = TRUE)
})
