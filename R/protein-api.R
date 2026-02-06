# Protein Domain API Functions
# ============================
# Functions for fetching protein domains and pLDDT scores from
# Ensembl, InterPro, TED, and AlphaFold REST APIs.

#' Fetch protein domains from Ensembl REST API
#'
#' Retrieves protein domain annotations for a given transcript from the
#' Ensembl REST API.
#'
#' @param transcript_id Character. Ensembl transcript ID (e.g., "ENST00000355451").
#' @param species Character. Species name (default: "human"). Currently unused
#'   but reserved for future multi-species support.
#' @param domain_sources Character vector. Filter to specific domain sources
#'   (e.g., c("Pfam", "SMART")). If NULL (default), returns all sources.
#' @param timeout Numeric. Request timeout in seconds (default: 30).
#'
#' @return A tibble with columns: domain, source, start, end.
#'   Returns NULL if no domains are found.
#'
#' @details
#' The function first looks up the protein ID associated with the transcript,
#' then fetches all protein features annotated as domains.
#'
#' @examples
#' \dontrun{
#' domains <- fetch_domains_ensembl("ENST00000355451")
#' domains <- fetch_domains_ensembl("ENST00000355451",
#'                                   domain_sources = c("Pfam", "SMART"))
#' }
#'
#' @export
fetch_domains_ensembl <- function(transcript_id,
                                  species = "human",
                                  domain_sources = NULL,
                                  timeout = 30) {
  # Validate input
  if (!is.character(transcript_id) || length(transcript_id) != 1) {
    stop("transcript_id must be a single character string")
  }

  # Get protein ID from transcript
  url <- paste0(
    "https://rest.ensembl.org/lookup/id/", transcript_id,
    "?content-type=application/json;expand=1"
  )

  response <- httr::GET(url, httr::timeout(timeout))

  if (httr::status_code(response) != 200) {
    stop("Failed to fetch transcript info. Check transcript ID.")
  }

  transcript_info <- jsonlite::fromJSON(
    httr::content(response, "text", encoding = "UTF-8")
  )
  protein_id <- transcript_info$Translation$id

  if (is.null(protein_id)) {
    stop("No protein associated with this transcript.")
  }

  # Fetch protein features (domains)
  url_features <- paste0(
    "https://rest.ensembl.org/overlap/translation/",
    protein_id,
    "?content-type=application/json;feature=protein_feature"
  )

  response_features <- httr::GET(url_features, httr::timeout(timeout))

  if (httr::status_code(response_features) != 200) {
    stop("Failed to fetch protein features.")
  }

  features <- jsonlite::fromJSON(
    httr::content(response_features, "text", encoding = "UTF-8")
  )

  if (length(features) == 0) {
    logger::log_info("No domains found for this protein.")
    return(NULL)
  }

  # Parse domain information
  domain_df <- features |>
    dplyr::select(
      domain = "description",
      source = "type",
      "start",
      "end"
    ) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$start)

  # Filter by domain sources if specified
  if (!is.null(domain_sources)) {
    domain_df <- domain_df |>
      dplyr::filter(.data$source %in% domain_sources)
  }


  tibble::as_tibble(domain_df)
}

#' Fetch protein domains from InterPro REST API
#'
#' Retrieves protein domain annotations for a given UniProt accession from the
#' InterPro REST API. InterPro aggregates domains from multiple sources
#' including Pfam, SMART, CDD, PROSITE, and others.
#'
#' @param uniprot_acc Character. UniProt accession (e.g., "P04637" for p53).
#' @param domain_sources Character vector. Filter to specific domain sources
#'   (e.g., c("pfam", "smart")). If NULL (default), returns all sources.
#'   Available sources include: pfam, smart, cdd, prosite, panther, etc.
#' @param entry_types Character vector. Filter to specific InterPro entry types
#'   (e.g., c("domain", "family")). If NULL (default), returns only "domain".
#' @param timeout Numeric. Request timeout in seconds (default: 30).
#'
#' @return A tibble with columns: domain, source, start, end, interpro_id.
#'   Returns NULL if no domains are found.
#'
#' @details
#' The InterPro API provides comprehensive domain annotations aggregated from
#' multiple member databases. Each domain entry includes the InterPro accession
#' and mappings to the underlying member database signatures.
#'
#' @examples
#' \dontrun{
#' # Fetch all domains for p53
#' domains <- fetch_domains_interpro("P04637")
#'
#' # Fetch only Pfam domains
#' domains <- fetch_domains_interpro("P04637", domain_sources = "pfam")
#' }
#'
#' @export
fetch_domains_interpro <- function(uniprot_acc,
                                   domain_sources = NULL,
                                   entry_types = "domain",
                                   timeout = 30) {
  # Validate input
  if (!is.character(uniprot_acc) || length(uniprot_acc) != 1) {
    stop("uniprot_acc must be a single character string")
  }

  # Build API URL
  url <- paste0(
    "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/",
    uniprot_acc,
    "?page_size=200"
  )

  response <- httr::GET(url, httr::timeout(timeout))

  if (httr::status_code(response) == 404) {
    logger::log_info("UniProt accession not found in InterPro.")
    return(NULL)
  }

  if (httr::status_code(response) != 200) {
    stop("Failed to fetch InterPro data. Status: ", httr::status_code(response))
  }

  data <- jsonlite::fromJSON(
    httr::content(response, "text", encoding = "UTF-8"),
    flatten = FALSE
  )

  if (is.null(data$results) || length(data$results) == 0) {
    logger::log_info("No domains found for this protein.")
    return(NULL)
  }

  # Parse results into domain data frame
  results <- lapply(seq_len(nrow(data$results)), function(i) {
    entry <- data$results[i, ]
    meta <- entry$metadata

    # Filter by entry type
    if (!is.null(entry_types) && !meta$type %in% entry_types) {
      return(NULL)
    }

    # Get protein locations
    proteins <- entry$proteins[[1]]
    if (is.null(proteins) || nrow(proteins) == 0) {
      return(NULL)
    }

    locations <- proteins$entry_protein_locations[[1]]
    if (is.null(locations) || length(locations) == 0) {
      return(NULL)
    }

    # Extract fragments from all locations
    fragments_list <- lapply(seq_len(nrow(locations)), function(j) {
      frags <- locations$fragments[[j]]
      if (is.null(frags) || nrow(frags) == 0) {
        return(NULL)
      }
      frags
    })

    fragments <- do.call(rbind, fragments_list)
    if (is.null(fragments) || nrow(fragments) == 0) {
      return(NULL)
    }

    # Get member database sources
    member_dbs <- names(meta$member_databases)
    source_str <- if (length(member_dbs) > 0) {
      paste(member_dbs, collapse = ",")
    } else {
      "interpro"
    }

    tibble::tibble(
      domain = meta$name,
      source = source_str,
      start = fragments$start,
      end = fragments$end,
      interpro_id = meta$accession
    )
  })

  domain_df <- do.call(rbind, results)

  if (is.null(domain_df) || nrow(domain_df) == 0) {
    logger::log_info("No domains found matching criteria.")
    return(NULL)
  }

  # Filter by domain sources if specified
  if (!is.null(domain_sources)) {
    domain_sources <- tolower(domain_sources)
    domain_df <- domain_df |>
      dplyr::filter(
        vapply(
          strsplit(.data$source, ","),
          function(x) any(tolower(x) %in% domain_sources),
          logical(1)
        )
      )

    if (nrow(domain_df) == 0) {
      logger::log_info("No domains found from specified sources.")
      return(NULL)
    }
  }

  domain_df <- domain_df |>
    dplyr::distinct() |>
    dplyr::arrange(.data$start)

  tibble::as_tibble(domain_df)
}

#' Fetch protein domains from TED (The Encyclopedia of Domains)
#'
#' Retrieves consensus structural domain annotations for a given UniProt
#' accession from the TED API. TED provides domain boundaries derived from
#' AlphaFold structures using consensus of Chainsaw, Merizo, and UniDoc.
#'
#' @param uniprot_acc Character. UniProt accession (e.g., "Q8IVD9").
#' @param consensus_levels Character vector. Filter to specific consensus
#'   levels: "high", "medium", or "low". If NULL (default), returns all.
#' @param timeout Numeric. Request timeout in seconds (default: 30).
#'
#' @return A tibble with columns: domain, source, start, end, ted_id,
#'   consensus_level, cath_label, plddt.
#'   Returns NULL if no domains are found.
#'
#' @details
#' TED (The Encyclopedia of Domains) provides structural domain annotations
#' derived from AlphaFold models. Domains are identified using a consensus
#' of three segmentation methods (Chainsaw, Merizo, UniDoc). High consensus
#' means all three methods agree; medium means two agree.
#'
#' @examples
#' \dontrun{
#' # Fetch all TED domains
#' domains <- fetch_domains_ted("Q8IVD9")
#'
#' # Fetch only high-confidence domains
#' domains <- fetch_domains_ted("Q8IVD9", consensus_levels = "high")
#' }
#'
#' @export
fetch_domains_ted <- function(uniprot_acc,
                              consensus_levels = NULL,
                              timeout = 30) {

  # Validate input
  if (!is.character(uniprot_acc) || length(uniprot_acc) != 1) {
    stop("uniprot_acc must be a single character string")
  }

  # Build API URL
  url <- paste0(
    "https://ted.cathdb.info/api/v1/uniprot/summary/",
    uniprot_acc,
    "?skip=0&limit=200"
  )

  response <- httr::GET(url, httr::timeout(timeout))

  if (httr::status_code(response) == 404) {
    logger::log_info("UniProt accession not found in TED.")
    return(NULL)
  }

  if (httr::status_code(response) != 200) {
    stop("Failed to fetch TED data. Status: ", httr::status_code(response))
  }

  data <- jsonlite::fromJSON(
    httr::content(response, "text", encoding = "UTF-8")
  )

  if (is.null(data$data) || length(data$data) == 0) {
    logger::log_info("No domains found for this protein in TED.")
    return(NULL)
  }

  domains <- data$data

  # Parse chopping field to get start/end positions
  chopping_parsed <- strsplit(domains$chopping, "-")
  starts <- vapply(chopping_parsed, function(x) as.integer(x[1]), integer(1))
  ends <- vapply(chopping_parsed, function(x) as.integer(x[2]), integer(1))

  # Build domain name from CATH label or TED ID
  domain_names <- ifelse(
    domains$cath_label != "-",
    paste0("CATH:", domains$cath_label),
    paste0("TED:", sub(".*_(TED\\d+)$", "\\1", domains$ted_id))
  )

  domain_df <- tibble::tibble(
    domain = domain_names,
    source = "TED",
    start = starts,
    end = ends,
    ted_id = domains$ted_id,
    consensus_level = domains$consensus_level,
    cath_label = domains$cath_label,
    plddt = domains$plddt
  )

  # Filter by consensus level if specified
  if (!is.null(consensus_levels)) {
    domain_df <- domain_df |>
      dplyr::filter(.data$consensus_level %in% consensus_levels)

    if (nrow(domain_df) == 0) {
      logger::log_info("No domains found with specified consensus levels.")
      return(NULL)
    }
  }

  domain_df <- domain_df |>
    dplyr::arrange(.data$start)

  tibble::as_tibble(domain_df)
}

#' Fetch per-residue pLDDT scores from AlphaFold Database
#'
#' Retrieves per-residue pLDDT (predicted Local Distance Difference Test)
#' confidence scores from the AlphaFold Database for a given UniProt accession.
#'
#' @param uniprot_acc Character. UniProt accession (e.g., "Q8IVD9").
#' @param timeout Numeric. Request timeout in seconds (default: 30).
#'
#' @return A tibble with columns: position, plddt, confidence_category.
#'   Returns NULL if no AlphaFold model is available.
#'
#' @details
#' pLDDT scores range from 0 to 100:
#' \itemize{
#'   \item Very high (>90): High accuracy for backbone and sidechains
#'   \item Confident (70-90): Correct backbone, some sidechain uncertainty
#'   \item Low (50-70): Low confidence
#'   \item Very low (<50): Often disordered or flexible regions
#' }
#'
#' Confidence categories: "H" (high), "M" (medium), "L" (low), "D" (disordered).
#'
#' @examples
#' \dontrun{
#' plddt <- fetch_plddt_alphafold("Q8IVD9")
#'
#' # Plot pLDDT along protein
#' library(ggplot2)
#' ggplot(plddt, aes(x = position, y = plddt)) +
#'   geom_line() +
#'   geom_hline(yintercept = c(50, 70, 90), linetype = "dashed", alpha = 0.5)
#' }
#'
#' @export
fetch_plddt_alphafold <- function(uniprot_acc, timeout = 30) {
  # Validate input
  if (!is.character(uniprot_acc) || length(uniprot_acc) != 1) {
    stop("uniprot_acc must be a single character string")
  }

  # First get the model metadata to find the confidence file URL
  meta_url <- paste0(
    "https://alphafold.ebi.ac.uk/api/prediction/",
    uniprot_acc
  )

  response <- httr::GET(meta_url, httr::timeout(timeout))

  if (httr::status_code(response) == 404) {
    logger::log_info("UniProt accession not found in AlphaFold Database.")
    return(NULL)
  }

  if (httr::status_code(response) != 200) {
    stop(
      "Failed to fetch AlphaFold metadata. Status: ",
      httr::status_code(response)
    )
  }

  meta <- jsonlite::fromJSON(
    httr::content(response, "text", encoding = "UTF-8")
  )

  # Handle case where API returns a list (take first entry)
  if (is.data.frame(meta)) {
    confidence_url <- meta$plddtDocUrl[1]
  } else if (is.list(meta)) {
    confidence_url <- meta[[1]]$plddtDocUrl
  } else {
    stop("Unexpected API response format")
  }

  if (is.null(confidence_url) || is.na(confidence_url)) {
    logger::log_info("No confidence data available for this protein.")
    return(NULL)
  }

  # Fetch the confidence JSON
  conf_response <- httr::GET(confidence_url, httr::timeout(timeout))

  if (httr::status_code(conf_response) != 200) {
    stop(
      "Failed to fetch confidence data. Status: ",
      httr::status_code(conf_response)
    )
  }

  conf_data <- jsonlite::fromJSON(
    httr::content(conf_response, "text", encoding = "UTF-8")
  )

  # Build tibble from the parallel arrays
  tibble::tibble(
    position = conf_data$residueNumber,
    plddt = conf_data$confidenceScore,
    confidence_category = conf_data$confidenceCategory
  )
}

#' Merge overlapping domain ranges
#'
#' Merges overlapping ranges for domains with the same name and source
#' using IRanges::reduce() for efficient interval operations.
#'
#' @param domain_df A data frame with columns: domain, source, start, end.
#'
#' @return A tibble with merged ranges (same columns as input).
#'
#' @examples
#' \dontrun{
#' domains <- fetch_domains_ensembl("ENST00000355451")
#' merged <- merge_overlapping_domains(domains)
#' }
#'
#' @export
merge_overlapping_domains <- function(domain_df) {
  if (is.null(domain_df) || nrow(domain_df) == 0) {
    return(domain_df)
  }

  merged <- domain_df |>
    dplyr::group_by(.data$domain, .data$source) |>
    dplyr::summarise(
      ranges = list(IRanges::reduce(IRanges::IRanges(.data$start, .data$end))),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      start = list(IRanges::start(.data$ranges)),
      end = list(IRanges::end(.data$ranges))
    ) |>
    dplyr::select(-"ranges") |>
    tidyr::unnest(cols = c("start", "end")) |>
    dplyr::ungroup()

  tibble::as_tibble(merged)
}
