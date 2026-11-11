#! /usr/bin/env Rscript

# Exports a GRanges object of the following format:
#     seqnames    ranges strand |      SYMBOL         ENSEMBL  refseq_mrna
#        <Rle> <IRanges>  <Rle> | <character>     <character>  <character>
# [1]     chr1    959256      - |       NOC2L ENSG00000188976    NM_015658

library(GenomicRanges)
library(dplyr)

# The EPD data can be downloaded from: https://epd.expasy.org/ftp/epdnew/H_sapiens/006/
epd_promoters <- read.table("./input/epd_promoter_ucsc.txt", sep = "\t") |>
  setNames(c("EPD_ID", "genome_build", "seq", "strand", "start")) |>
  dplyr::filter(genome_build == "hg38") |>
  dplyr::mutate(SYMBOL = gsub(EPD_ID, pattern = "_.*", replacement = "")) |>
  dplyr::select(-genome_build)
epd_cross_references <- read.table("./input/epd_cross_references.txt", sep = "\t", quote = "") |>
  setNames(c("ENSEMBL", "SYMBOL", "TXNAME", "description")) |>
  dplyr::select(-description) |>
  # NM_ in the NCBI system is only coding trascripts. NR_ is for non-coding.
  dplyr::filter(grepl(ENSEMBL, pattern = "ENSG.*") & grepl(TXNAME, pattern = "NM_.*")) |> # Remove where identifiers do not match norms
  dplyr::distinct(SYMBOL, .keep_all = TRUE)
epd_promoter_type <- read.table("./input/epd_promoter_coordinate.txt", sep = "\t") |>
  setNames(c("EPD_ID", "seq_ncbi", "start", "strand", "species", "EPD_TYPE")) |>
  dplyr::select(-c(seq_ncbi, start, strand, species)) |>
  filter(EPD_TYPE != "unclassified")

# Add the type of promoter to the list of promoters
epd_promoters <- inner_join(epd_promoters, epd_promoter_type, by = join_by(EPD_ID))

# Combine the different EPD data and convert to a GRanges object
epd <- inner_join(epd_promoters, epd_cross_references, by = join_by(SYMBOL))
message(paste0("Of ", dim(epd_promoters)[1], " EPD promoters, ", dim(epd)[1], " have proper transcript and gene annotations"))
epd_granges <- GRanges(seqnames = epd$seq, ranges = IRanges(start = epd$start, end = epd$start), strand = epd$strand)
mcols(epd_granges, use.names = TRUE) <- dplyr::select(epd, SYMBOL, ENSEMBL, TXNAME, EPD_ID, EPD_TYPE)

# Save as RDS object
saveRDS(epd_granges, "./input/epd_promoters_with_gene_annotation.rds")

