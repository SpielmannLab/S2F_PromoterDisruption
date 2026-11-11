#! /usr/bin/env Rscript

# This script creates a full promoter per gene - which covers the entire range of EPD, pELS and PLS, and CAGE peaks.
# Also creates the different regions which needs to be disrupted and not.

library(universalmotif)
library(GenomicRanges)
library(dplyr)

epd_prom_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/epd_promoters_grangeslist.rds")
cage_prom_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/cage_promoters_grangeslist.rds")
prox_prom_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/proximal_promoters_grangeslist.rds")
prox_enh_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/proximal_enhancers_grangeslist.rds")

# Make sure the lists are matching
all_identical <- function(...) {
  x <- list(...)
  res <- c()
  for (i in seq.int(from = 2, to = length(x))) {
    res[i - 1] <- all(identical(x[1], x[i]))
  }
  return(all(res))
}
res <- all_identical(names(prox_prom_l), names(prox_enh_l), names(cage_prom_l), names(epd_prom_l))
if (res) {
  message("Loading of all GRangesLists successfull")
} else {
  stop("The gene names across the lists are not identical")
}

# Create full promoter from all the promoter and enhancer elements. Ignore strand at this point.
full_prom_l <-  mendoapply(range, cage_prom_l, epd_prom_l, prox_prom_l, prox_enh_l, MoreArgs = list(ignore.strand = TRUE))

### NOW, for each gene, define regions to vary and not to vary ###

# Regions of promoters to be varied
regions_to_vary <- mendoapply(c, cage_prom_l, epd_prom_l, prox_prom_l, prox_enh_l) |>
  endoapply(reduce, ignore.strand = TRUE) |>
  unlist()
mcols(regions_to_vary) <- data.frame(shuffle = TRUE)
regions_to_vary_l <- regions_to_vary |>
  split(f = names(regions_to_vary))

### Regions of promoters not to be varied
regions_not_to_vary <- mendoapply(setdiff, full_prom_l, regions_to_vary_l) |>
  unlist()
mcols(regions_not_to_vary) <- data.frame(shuffle = FALSE)
regions_not_to_vary_l <- regions_not_to_vary |>
  split(f = names(regions_not_to_vary))

# Save them all
dir.create("full_prom_to_test")
saveRDS(full_prom_l, file = "./full_prom_to_test/full_promoters_to_test_grangerslist.rds")
saveRDS(regions_to_vary_l, file = "./full_prom_to_test/regions_to_vary_to_test_grangerslist.rds")
saveRDS(regions_not_to_vary_l, file = "./full_prom_to_test/regions_not_to_vary_to_test_grangerslist.rds")
