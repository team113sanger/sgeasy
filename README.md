# sgeasy

Downstream analysis tools for Saturation Genome Editing (SGE) screens.

## Installation

```r
# Install from GitHub (when available)
# devtools::install_github("your-org/sgeasy")

# Or install locally
devtools::install()
```
## Overview

sgedown provides a complete workflow for analyzing SGE screen data:

- **Normalization**: Size factor estimation using control oligos (synonymous/intronic variants)
- **Differential Analysis**: DESeq2-based differential abundance across time points
- **Visualization**: Quality control plots and sample comparisons

## Quick Start

```r
library(sgeasy)

# Create normalization matrix from neutral variants
norm_matrix <- create_normalization_matrix(annotated_data)

# Create count matrix for all variants
count_matrix <- create_count_matrix(annotated_data)

# Run differential abundance analysis
results <- run_differential_analysis(
  count_matrix = count_matrix,
  normalization_matrix = norm_matrix,
  sample_metadata = metadata
)

# Access results
deseq_results <- results$results
rlog_data <- results$rlog
summary_table <- results$contrast_summary
```

## Documentation

See the package vignette for a complete workflow tutorial:

```r
vignette("sge-analysis-workflow", package = "sgedown")
```
