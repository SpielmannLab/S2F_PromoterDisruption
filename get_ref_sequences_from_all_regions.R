library(universalmotif)
library(GenomicRanges)
library(dplyr)
library(BSgenome)

# Get the reference sequences for all the regions (to disrupt and not to disrupt).

# Load the regions to disrupt and not to disrupt
regions_to_vary_l <- readRDS("./full_prom_to_test/regions_to_vary_to_test_grangerslist.rds")
regions_not_to_vary_l <- readRDS("./full_prom_to_test/regions_not_to_vary_to_test_grangerslist.rds")

# Define the genome - note this has to be separately installed
hg38 <- getBSgenome("hg38")

# Extract sequences of the regions from the genome.
regions_not_to_vary <- unlist(regions_not_to_vary_l, use.names = FALSE)
mcols(regions_not_to_vary) <- data.frame(ENSEMBL = names(regions_not_to_vary))
sequences <- Biostrings::getSeq(hg38, regions_not_to_vary)
regions_not_to_vary$ref_seq <- sequences

regions_to_vary <- unlist(regions_to_vary_l, use.names = FALSE)
mcols(regions_to_vary) <- data.frame(ENSEMBL = names(regions_to_vary))
sequences <- Biostrings::getSeq(hg38, regions_to_vary)
regions_to_vary$ref_seq <- sequences

regions_not_to_vary_l <- regions_not_to_vary |>
  split(f = factor(regions_not_to_vary$ENSEMBL))
regions_to_vary_l <- regions_to_vary |>
  split(f = factor(regions_to_vary$ENSEMBL))

# Save them all
saveRDS(regions_to_vary_l, file = "./full_prom_with_seqs/regions_to_vary_to_test_grangerslist.rds")
saveRDS(regions_not_to_vary_l, file = "./full_prom_with_seqs/regions_not_to_vary_to_test_grangerslist.rds")
