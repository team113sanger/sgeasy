# Function to estimate size factors (typically run for control oligos only)
estimate_control_size_factors <- function (countData = countData, colData = colData, design = design, minRowSum = 10, ref = NULL ) {
  logger::log_info( "Estimating control size factors...")
  dds <- DESeqDataSetFromMatrix( countData = countData, 
                                 colData = colData, 
                                 design = as.formula( design ) )
  
  dds <- dds[ rowSums( counts( dds ) ) > minRowSum, ]
  
  if ( ! is.null( ref ) ) {
    dds$condition <- relevel( dds$condition, ref = ref )
  }
  
  control_size_factors <- sizeFactors( estimateSizeFactors( dds ) )
  
  return( control_size_factors )
}


create_normalisation_matrix <- function(df, 
                                        neutral_variants = c("intron_variant", "synonymous_variant")){
dplyr::filter(df, Consequence %in% neutral_variants) |> 
dplyr::select(SEQUENCE, starts_with("count")) |> 
tibble::column_to_rownames("SEQUENCE") |>
as.matrix()
}


create_input_matrix <- function(df){
    df |>
    dplyr::select(SEQUENCE, starts_with("count")) |> 
    tibble::column_to_rownames("SEQUENCE") |>
    as.matrix()
}


filter_matrices <- function(mat, min_counts = 10, max_counts = Inf) {
    mat |>
    rowwise() |>
    mutate(total_counts = sum(c_across(contains("pel")), na.rm = TRUE)) |> 
    dplyr::filter(between(total_counts, min_counts, max_counts)) |>
    dplyr::select(-contains("_0_"))
}

ggpairs_ext <- function(data, mapping, pts=list(), smt=list(), ...){
  ggplot(data = data, mapping = mapping, ...) +
    do.call(geom_point, pts) +
    do.call(geom_smooth, smt)
}

## Unique DEG table creation

create_unique_deg_tab <- function(table, name){
    deg(table, "raw") |> 
    as_tibble(rownames = "SEQUENCE") |> 
    rename_with(~paste0(.x, "_", name), .cols = -SEQUENCE)
}



get_control_size_factors <- function(matrix, condition_data, 
                                    ref = "Day4", 
                                    column_name = "condition"){
    logger::log_info( "...")
    control_size_factors <- estimate_control_size_factors(
                    countData = matrix,
                    colData = condition_data,
                    design = glue("~ {column_name}"),
                    minRowSum = 10,
                    ref = ref)
    return(control_size_factors)
}



diff_abundance_analysis <- function(matrix, norm, screen_metadata){
  logger::log_info("Running DESeq2 differential abundance analysis on timepoint contrasts")

  condition_data <- data.frame(condition = 
  screen_metadata$condition[match(colnames(matrix), 
                                  paste0("count_", screen_metadata$supplier_name))],
                                  duration = screen_metadata$duration[match(colnames(matrix), 
                                  paste0("count_", screen_metadata$supplier_name))]
) |> mutate(condition = factor(condition, levels = c("Day4", "Day7", "Day10", "Day15"))) 

    size_factors <- get_control_size_factors(norm, condition_data)
    dds <- DESeqDataSetFromMatrix(countData = matrix, 
            colData = condition_data, 
            design = ~condition)
    dds <- estimateSizeFactors(dds)


    logger::log_info("Running DESeq2 differential abundance analysis on continous data for Rate estimation")
    ## Update dds size factors with control size factors
    sizeFactors(dds) <- size_factors[[1]]
    dds <- DESeq(dds)
    rld <- rlog(dds)
    res <- results(dds) |> 
    as_tibble(rownames = "SEQUENCE")



    z_score <- assay(rld) |> 
    as.matrix() |> 
    t() |> 
    scale() |> 
    t() |>
    as_tibble(rownames = "SEQUENCE")
    colnames(z_score) <- paste0(colnames(z_score), "_z_score")
    colnames(z_score)[1] <- "SEQUENCE"

    available_contrasts <- resultsNames(dds)[2:length(resultsNames(dds))]
    # Shirnkage type set to normal, apeglm would be better for RNASeq - but have selected no shrinkage is table
    table_wald <- degComps(dds, combs = "condition", contrast = available_contrasts, alpha = 0.05, skip = FALSE, type = "normal", pairs = FALSE, fdr = "default")
  
    first_contrast <- names(table_wald)[1]

    all_contrast_summary <- imap(table_wald,create_unique_deg_tab) |>
    purrr::reduce(left_join, by = "SEQUENCE") |> 
    dplyr::rename(baseMean = paste0("baseMean_", first_contrast)) |>
    dplyr::select(-starts_with("baseMean_"))

    # all_contrast_summary |> write_tsv("DEG_summary.tsv")

    all_contrast_summary <- all_contrast_summary |> 
    left_join(z_score, by = "SEQUENCE")


    # Continuous analysis 

    continous_size_factors <- get_control_size_factors(norm, 
                                condition_data, 
                                column_name = "duration", 
                                ref = NULL)
    continous_dds <- DESeqDataSetFromMatrix(countData = matrix, 
            colData = condition_data, 
            design = ~duration)
    continous_dds <- estimateSizeFactors(continous_dds)
    sizeFactors(continous_dds) <- continous_size_factors[[1]]
    continous_dds <- DESeq(continous_dds)
    wald_continuous <- degComps(continous_dds, 
                                combs = "duration", alpha = 0.05, 
                                skip = FALSE, type = "normal", 
                                pairs = FALSE, fdr = "default")
    
    rate <- wald_continuous[[1]] |> 
    as_tibble(rownames = "SEQUENCE") |> 
    rename_with(~paste0(.x, "_", "continous"), .cols = -SEQUENCE) |>
    dplyr::select(-starts_with("baseMean_"))

    all_contrast_summary <- left_join(all_contrast_summary, rate, by = "SEQUENCE")



return(list(res, rld, all_contrast_summary))
}


scatter_plots <- function(dataset, sample_list){
    subset <- dataset |> 
    as_tibble(rownames = "SEQUENCE") |>
    dplyr::select(all_of(sample_list)) 
    print(ggpairs(subset, 
    lower = list(continuous = wrap(ggpairs_ext, 
                   pts=list(size=0.4,colour="black"), 
                   smt=list(method="lm", se=F, size=0.2, 
                   colour="blue"))), 
                   title = "Regularized-log Transformed Read Count") + 
                   theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))

}

#function to add text to column names
appendDataFrameColumns<-function(df, prefix='', suffix='', sep='')
{
  colnames(df) <- paste(prefix, colnames(df), suffix, sep=sep)
  
  return(df)
}



calculate_median_scores <- function(df, 
                                    neutral_variants = c("intron_variant", "synonymous_variant")){
  
  scores <- df |> 
  dplyr::filter(Consequence %in% neutral_variants) |>
  mutate(
    across(
      starts_with("log2FoldChange"),
      ~ median(.x, na.rm = TRUE),
      .names = "median_{.col}"
    )) |> 
    select(contains("median")) |> 
    slice_head(n=1)

    df <- bind_cols(df, scores)
    return(df)
}


add_adj_columns <- function(data, suffix) {
    lfc_col <- paste0("log2FoldChange_", suffix)
    se_col <- paste0("lfcSE_", suffix)
    median_col <- paste0("median_",lfc_col)
    data <- data |>
      mutate(
        "adj_lfc_{suffix}" := .data[[lfc_col]] - .data[[median_col]],
        "adj_score_{suffix}" := .data[[paste0("adj_lfc_", suffix)]] / .data[[se_col]],
        "pval_{suffix}" := pnorm(abs(.data[[paste0("adj_score_", suffix)]]), lower.tail = FALSE) * 2,
        "FDR_{suffix}" := p.adjust(.data[[paste0("pval_", suffix)]], method = "BH")
      )
      return(data)
  }

adjust_all_contrasts <- function(data) {                                                                
  suffixes <- names(data) |>                                                                                   
  str_subset("^log2FoldChange_") |>                                                                            
  str_remove("^log2FoldChange_")    
                                                                                                               
  output <- purrr::reduce(suffixes, \(d, s) add_adj_columns(d, s), .init = data)                                           
    return(output)
  }           

classify_variants <- function(df){
  df |> 
  mutate(variant_class = case_when(
    Consequence %in% "synonymous_variant" ~ "synonymous",
    Consequence %in% "missense_variant" ~ "missense",
    Consequence %in% "stop_gained" ~ "nonsense",
    Consequence %in% "intron_variant" ~ "intronic",
    TRUE ~ "other"
  ))
}

add_classification <- function(data, suffix, fdr_threshold = 0.01) {
    fdr_col <- paste0("FDR_", suffix)
    lfc_col <- paste0("adj_lfc_", suffix)

    data |>
      mutate(
        "functional_classification_{suffix}" := case_when(
          .data[[fdr_col]] < fdr_threshold & .data[[lfc_col]] < 0 ~ "depleted",
          .data[[fdr_col]] < fdr_threshold & .data[[lfc_col]] > 0 ~ "enriched",
          .data[[fdr_col]] >= fdr_threshold ~ "unchanged",
          TRUE ~ NA_character_
        )
      )
  }

classify_all_contrasts <- function(data, fdr_threshold = 0.01) {
suffixes <- names(data) |>
    str_subset("^FDR_") |>
    str_remove("^FDR_")

purrr::reduce(suffixes, \(d, s) add_classification(d, s, fdr_threshold), .init = data)
}



slim_consequence_fun <- function(df) {
    df <- df |>
      mutate(
        slim_consequence = case_when(
          # More specific conditions first (with mutator requirements)
          str_detect(Consequence, "inframe_deletion") & mutator == "custom" ~ "clinical_inframe_deletion",
          str_detect(Consequence, "inframe_insertion") & mutator == "custom" ~ "clinical_inframe_insertion",
          str_detect(Consequence, "stop_retained_variant") & mutator %in% c("snvre", "snv", "custom") ~ "synonymous",
          str_detect(Consequence, "start_lost") & mutator %in% c("snvre", "snv", "custom", "ala") ~ "start_lost",
          str_detect(Consequence, "stop_lost") & mutator %in% c("snvre", "snv", "custom", "ala") ~ "stop_lost",
          str_detect(mutator, "inframe") ~ "codon_deletion",

          # Simple pattern matches
          str_detect(Consequence, "stop_gained") ~ "stop_gained",
          str_detect(Consequence, "frameshift_variant") ~ "frameshift",
          str_detect(Consequence, "splice_acceptor_variant") ~ "splice_acceptor",
          str_detect(Consequence, "splice_donor_variant") ~ "splice_donor",
          str_detect(Consequence, "missense_variant") ~ "missense",
          str_detect(Consequence, "synonymous_variant") ~ "synonymous",
          str_detect(Consequence, "intron_variant") ~ "intron",
          str_detect(Consequence, "UTR_variant") ~ "UTR",

          TRUE ~ NA_character_
        )
      )
      return(df)
  }

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

  classification_shapes <- c("unchanged" = 0, "depleted" = 6, "enriched" = 3)


# Function to create one plot
  plot_condition <- function(data, suffix, id) {
    score_col <- paste0("adj_score_", suffix)
    class_col <- paste0("functional_classification_", suffix)
    fdr_col <- paste0("FDR_", suffix)

    ggplot(data, aes(
      x = vcf_pos,
      y = .data[[score_col]],
      shape = .data[[class_col]],
      colour = slim_consequence,
      fill = slim_consequence
    )) +
      geom_point(aes(alpha = .data[[fdr_col]] <= 0.01), size = 1.5) +
      scale_colour_manual(values = consequence_colours) +
      scale_fill_manual(values = consequence_colours) +
      scale_shape_manual(values = classification_shapes) +
      scale_alpha_discrete(labels = c("FDR > 0.01", "FDR <= 0.01")) +
      xlab("GRCh38 coordinate") +
      ylab(paste0("z-score (", suffix, ")")) +
      ggtitle(paste("z-score Plot:", id, "-", suffix)) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = "black", size = 10),
        legend.title = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        legend.margin = margin()
      )
  }

  # Generate all plots
  plot_all_conditions <- function(data, id) {
    suffixes <- names(data) |>
      str_subset("^adj_score_") |>
      str_remove("^adj_score_")

    plots <- map(suffixes, ~ plot_condition(data, .x, id))
    names(plots) <- paste0(suffixes)
    plots
  }


recalculate_screen_statistics <- function(screen){
    screen <- screen |>
    calculate_median_scores() |> 
    adjust_all_contrasts() |> 
    classify_all_contrasts()
    
    return(screen)
}


write_waterfall_plots <- function(plot_set, dest = "results"){
    dir.create(dest)
    plot_set |>
    list_flatten() |>
    iwalk(~ {
        png(glue::glue("{dest}/{.y}.png"), width = 900, height = 600)
        print(.x)
        dev.off()})
}

post_process <- function(dataframe, targeton_id, annotation_file){
    dataframe <- dataframe |>
                mutate(Targeton_ID = targeton_id) |>
                left_join(annotation_file, by = c("SEQUENCE" = "Seq", "Targeton_ID" = "Targeton_ID")) |>
                recalculate_screen_statistics() |>
                slim_consequence_fun()
}