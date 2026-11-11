#! /usr/bin/env Rscript

# Exports a GRanges object of the following format:
#     seqnames    ranges strand |      SYMBOL         ENSEMBL  refseq_mrna
#        <Rle> <IRanges>  <Rle> | <character>     <character>  <character>
# [1]     chr1    959256      - |       NOC2L ENSG00000188976    NM_015658

library(GenomicRanges)
library(dplyr)

cage_peaks <- read.table("./input/hg38_fantom5_cage_peaks.bed", sep = "\t", skip = 1) |>
  setNames(c("seq", "start", "end", "CAGE_ID", "score", "strand", "thick_start", "thick_end")) |>
  dplyr::filter(grepl(CAGE_ID, pattern = "@")) |>
  dplyr::mutate(SYMBOL = gsub(CAGE_ID, pattern = ".*@", replacement = "")) |>
  select(-thick_start, -thick_end, -score)

epd_cross_references <- read.table("./input/epd_cross_references.txt", sep = "\t", quote = "") |>
  setNames(c("ENSEMBL", "SYMBOL", "TXNAME", "description")) |>
  dplyr::select(-description) |>
  # NM_ in the NCBI system is only coding trascripts. NR_ is for non-coding.
  dplyr::filter(grepl(ENSEMBL, pattern = "ENSG.*") & grepl(TXNAME, pattern = "NM_.*")) |> # Remove where identifiers do not match norms
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# Combine the different EPD data and convert to a GRanges object
cage <- inner_join(cage_peaks, epd_cross_references, by = join_by(SYMBOL), unmatched = "drop")
cage_granges <- GRanges(seqnames = cage$seq, ranges = IRanges(start = cage$start, end = cage$end), strand = cage$strand)
mcols(cage_granges, use.names = TRUE) <- dplyr::select(cage, SYMBOL, ENSEMBL, TXNAME, CAGE_ID)

# Save as RDS object
saveRDS(cage_granges, "./input/cage_peaks_with_gene_annotation.rds")

