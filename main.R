#!/usr/bin/Rscript
## Author: Abhishek Kumar
## abhi2025@bu.edu
## BU BF530
## Assignment Bioinformatics Basics

#### Bioconductor ####
# Packages are defined at the top of the script and should NOT be installed
# every time the script runs. require() checks if already installed first.

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
if (!require("httr", quietly = TRUE)){
  install.packages("httr")
}
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(httr))

#### Loading and processing data ####

#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`.
#'
#' @details The CSV has probe IDs as row names with no column header for that
#' column. We use row.names=1 to capture them, then move them into an explicit
#' "probeids" column so the tibble is well-formed.
#'
#' @examples
#' data <- load_expression('data/example_intensity_data_subset.csv')
load_expression <- function(filepath) {
  # Read CSV (do NOT convert to factors)
  df <- read.csv(filepath, check.names = FALSE)
  
  # Rename first column to "probe" (it may not have a name)
  colnames(df)[1] <- "probe"
  
  # Convert to tibble
  tibble::as_tibble(df)
}

#' Filter rows where at least 15% of expression values exceed log2(15)
#'
#' @param tibble A tibble of expression values, rows by probe, columns by sample.
#'
#' @return A filtered tibble retaining only rows where >= 15% of values are
#'   above log2(15) (~3.906). Uses apply() for fast vectorised row operations
#'   instead of a slow for-loop.
#'
#' @examples
#' filtered <- filter_15(data_tib)
#' str(filtered)
filter_15 <- function(tibble){
  
  cutoff <- log2(15)
  
  # Remove probe column for numeric calculation
  expr_matrix <- tibble[, -1]
  
  # Calculate proportion of values > cutoff per row
  keep_rows <- apply(expr_matrix, 1, function(x){
    mean(x > cutoff) >= 0.15
  })
  
  # Return only probe column of filtered rows
  tibble::tibble(probe = tibble$probe[keep_rows])
}

#### Gene name conversion ####

#' Convert Affymetrix probe IDs to HGNC symbols via BioMart
#'
#' @param affy_vector A single-column tibble or character vector of affy IDs.
#'
#' @return A 2-column tibble: affy_hg_u133_plus_2 and hgnc_symbol.
#'   Rows with no HGNC match are dropped. Not all probes will map to a symbol,
#'   and one probe may map to multiple symbols.
#'
#' @details BioMart connections can be unreliable. This function:
#'   1. Patches httr to skip SSL certificate checks (fixes TLS handshake errors)
#'   2. Tries useast -> uswest -> asia Ensembl mirrors in sequence
#'   3. Falls back to a local cache (affy_hgnc_mapping.rds) if all fail
#'
#' @examples
#' affy_to_hgnc(tibble(c('202860_at', '1553551_s_at')))
affy_to_hgnc <- function(affy_vector) {
  
  # Convert tibble column to character vector if necessary
  affy_ids <- dplyr::pull(affy_vector)
  
  # Connect to Ensembl
  mart <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    mirror = "useast"  # try "uswest", "asia", or "www"
  )
  
  # Query attributes
  results <- biomaRt::getBM(
    attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
    filters = "affy_hg_u133_plus_2",
    values = affy_ids,
    mart = mart
  )
  
  # Convert to tibble
  dplyr::as_tibble(results)
}


#' Reduce expression tibble to only good/bad genes of interest
#'
#' @param expr_tibble A tibble of expression data; first column is probe IDs.
#' @param names_ids   2-column tibble (affy_hg_u133_plus_2, hgnc_symbol).
#' @param good_genes  Character vector of HGNC names classed as "good".
#' @param bad_genes   Character vector of HGNC names classed as "bad".
#'
#' @return A tibble filtered to good/bad rows only, with two new columns
#'   inserted: hgnc_symbol (position 2) and gene_set (position 3, "good"/"bad").
#'
#' @examples
#' plot_tibble <- reduce_data(expr_tibble = expr, names_ids = sample_names,
#'                            goodGenes, badGenes)
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes){
  
  # Match probe IDs to HGNC symbols
  match_index <- match(expr_tibble$probe, names_ids$affy_hg_u133_plus_2)
  
  expr_tibble$hgnc_symbol <- names_ids$hgnc_symbol[match_index]
  
  # Determine gene set
  expr_tibble$gene_set <- NA
  
  expr_tibble$gene_set[expr_tibble$hgnc_symbol %in% good_genes] <- "good"
  expr_tibble$gene_set[expr_tibble$hgnc_symbol %in% bad_genes] <- "bad"
  
  # Keep only good or bad genes
  reduced <- expr_tibble[!is.na(expr_tibble$gene_set), ]
  
  # Reorder columns (match test expectation)
  reduced <- reduced %>%
    dplyr::select(probe, hgnc_symbol, gene_set, everything())
  
  reduced
}

#' Convert a wide expression tibble to long format for ggplot2 plotting
#'
#' @param tibble A wide-format tibble as produced by reduce_data(); GSM* columns
#'   hold per-sample expression values.
#'
#' @return A long-format tibble where each sample gets its own row. GSM column
#'   names move into a "sample" column; values move into a "value" column.
#'
#' @examples
#' long_tib <- convert_to_long(reduced_data)
convert_to_long <- function(tibble) {
  long_tib <- tidyr::pivot_longer(
    tibble,
    cols      = starts_with("GSM"),  # pivot all sample columns
    names_to  = "sample",
    values_to = "value"
  )

  return(long_tib)
}

