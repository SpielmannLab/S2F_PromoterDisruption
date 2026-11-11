#!/usr/bin/env Rscript
library(universalmotif)
library(GenomicRanges)
library(dplyr)

n_random <- 4 # Not yet implemented! This will let create multiple ALT allele for multiple trials.
# set.seed(121) # Optional

k_let <- 2 # Do you want to maintain the frequency of dinucleotide, trinucleotide, or k-let nucleotide?
# e.g., k=2 will maintain the dinucleotide frequency

regions_to_vary_l <- readRDS(file = "./full_prom_with_seqs/regions_to_vary_to_test_grangerslist.rds")
regions_not_to_vary_l <- readRDS(file = "./full_prom_with_seqs/regions_not_to_vary_to_test_grangerslist.rds")
cage_scoring_region_l <- readRDS(file = "./full_prom_to_test/cage_scoring_regions_grangerslist.rds")
full_prom_l <- readRDS("./full_prom_to_test/full_promoters_to_test_grangerslist.rds")
rep_trans_l <- readRDS(file = "./gene_seggregated_enhs_proms_trans/representative_transcripts_grangeslist.rds")

# Make sure the lists are matching
all_identical <- function(...) {
  x <- list(...)
  res <- c()
  for (i in seq.int(from = 2, to = length(x))) {
    res[i - 1] <- all(identical(x[1], x[i]))
  }
  return(all(res))
}
res <- all_identical(names(regions_to_vary_l), names(regions_not_to_vary_l), names(rep_trans_l), names(full_prom_l), names(cage_scoring_region_l))
if (res) {
  message("Loading of all GRangesLists successfull")
} else {
  stop("The gene names across the lists are not identical")
}

regions_not_to_vary <- unlist(regions_not_to_vary_l, use.names = FALSE)
regions_to_vary <- unlist(regions_to_vary_l, use.names = FALSE)

# Now create the variants
regions_to_vary$alt_seq <- shuffle_sequences(regions_to_vary$ref_seq, k = k_let)
regions_not_to_vary$alt_seq <- regions_not_to_vary$ref_seq

all_regions <- c(regions_to_vary, regions_not_to_vary) |> sort()
all_regions_l <- split(all_regions, f = names(all_regions))

# This is a reduce function, but carries over the sequence in the metadata column
combine_regions <- function(gr) {
  sorted_gr <- sort(gr)
  reduced_gr <- reduce(sorted_gr)
  mcols(reduced_gr) <- data.frame(ref_seq = as.character(unlist(sorted_gr$ref_seq)),
    alt_seq = as.character(unlist(sorted_gr$alt_seq)))
  return(reduced_gr)
}

prom_variants <- endoapply(FUN = combine_regions, X = all_regions_l) |>
  unlist()
prom_variants$variant_id <- paste0(
  seqnames(prom_variants),
  ":",
  start(prom_variants),
  "-",
  end(prom_variants),
  "(",
  names(prom_variants),
  ")"
)

# Now add the cage scoring regions
cage_scoring_region <- unlist(cage_scoring_region_l)
rep_trans <- unlist(rep_trans_l)

# Make sure the orders are the same
res <- all_identical(names(prom_variants), names(cage_scoring_region), names(rep_trans))
if (res) {
  message("Ready to write to file")
} else {
  stop("The gene names across the GRanges objects are not identical")
}
prom_variants$cage_score_start <- start(cage_scoring_region)
prom_variants$cage_score_end <- end(cage_scoring_region)
prom_variants$cage_score_strand <- strand(cage_scoring_region)
prom_variants$rna_score_start <- start(rep_trans)
prom_variants$rna_score_end <- end(rep_trans)
prom_variants$rna_score_strand <- strand(rep_trans)
prom_variants$ENSEMBL <- names(cage_scoring_region)

write.table(prom_variants, file = "./prom_variants_to_test.tsv", sep = "\t", quote = FALSE)
