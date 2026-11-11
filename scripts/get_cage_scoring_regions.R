#! /usr/bin/env Rscript

# This script creates regions for cage scoring - which covers the entire range of EPD and PLS, then adds a buffer 200bp to either ends

library(universalmotif)
library(GenomicRanges)
library(dplyr)

epd_prom_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/epd_promoters_grangeslist.rds")
prox_prom_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/proximal_promoters_grangeslist.rds")

# Make sure the lists are matching
all_identical <- function(...) {
  x <- list(...)
  res <- c()
  for (i in seq.int(from = 2, to = length(x))) {
    res[i - 1] <- all(identical(x[1], x[i]))
  }
  return(all(res))
}
res <- all_identical(names(prox_prom_l), names(epd_prom_l))
if (res) {
  message("Loading of all GRangesLists successfull")
} else {
  stop("The gene names across the lists are not identical")
}

# Create a range for scoring the CAGE signal promoter from all the EPD (with 60 and 10 bp expansion) and PLS (ENCODE) elements.
cage_scoring_region_no_strand_l <- mendoapply(range, epd_prom_l, prox_prom_l, MoreArgs = list(ignore.strand = TRUE))
# We need the strand information! and also need to expand the region by 200 bp at the start and the end.
cage_scoring_region_l <- mendoapply(
  FUN = function(x, y) {
    strand(x) <- unique(strand(y))
    start(x) <- start(x) - 200
    end(x) <- end(x) + 200
    return(x)
  },
  cage_scoring_region_no_strand_l,
  epd_prom_l
)

# Save it
saveRDS(cage_scoring_region_l, file = "./full_prom_to_test/cage_scoring_regions_grangerslist.rds")
