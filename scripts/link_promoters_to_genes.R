#! /usr/bin/env Rscript

library(txdbmaker) # Create TxDb from UCSC gene tracks
library(AnnotationHub)
library(rtracklayer) # Enables importing bedfiles etc
library(GenomicRanges) # For interesction of regions
library(BSgenome)
library(ggplot2)

# This script links the EPD TSS/promoters to pELS and PLS (from ENCODE).
# It creates GRangesList object for all the elements grouped by ENSEMBL IDs.

# Allocate core, promoters, TSS, proximal promoters, and proximal enhancers to refseq select transcripts/promoters #####
prox_prom_enh <- readRDS("../input/encode_promoters.rds")
prox_prom <- prox_prom_enh[prox_prom_enh$cCRE_class == "Promoter"]
prox_enh <- prox_prom_enh[prox_prom_enh$cCRE_class == "Proximal enhancer"]
rep_trans <- readRDS("../input/representative_transcripts.rds")
epd_tss <- readRDS("../input/epd_promoters_with_gene_annotation.rds")
cage_tss <- readRDS("../input/cage_peaks_with_gene_annotation.rds")
genes_to_use <- "custom_gene_list" # options are "random" or "custom_genes_list"
custom_gene_list <- c("ENSG00000100336", "ENSG00000076344", "ENSG00000204920") # APOL4 (Fig4), RGS11 (Supp Fig 8), ZNF155 (Supp Fig 8)

# Step 1. Allocate PLS to nearest gene(s) based on all TSS in EPD database.
dists <- distanceToNearest(prox_prom, epd_tss, ignore.strand = TRUE, select = "all")
# lets plot a histogram of distances between prox_prom and nearest genes
p1 <- ggplot(as.data.frame(mcols(dists)), aes(x = distance + 1)) +
  geom_density() +
  scale_x_log10(labels = function(x) {
    format(x, scientific = FALSE)
  }) +
  ggtitle("Before filtering")
# remove all PLS with distance more than 1000 bp from the nearerst EPD-TSS. Based on the histogram
dists <- dists[mcols(dists)$distance < 1000]
p2 <- ggplot(as.data.frame(mcols(dists)), aes(x = distance + 1)) +
  geom_density() +
  scale_x_log10(labels = function(x) {
    format(x, scientific = FALSE)
  }) +
  ggtitle("After filtering")
ggsave(plot = cowplot::plot_grid(p1, p2, nrow = 2), file = "./plots/distance_distribtion_between_PLS_and_TSS.pdf")
prox_prom <- prox_prom[from(dists)]
new_metadata_columns <- mcols(epd_tss[to(dists)])
mcols(prox_prom) <- new_metadata_columns

# Step 2. Allocate pELS to gene within 2kb based on all TSS in EPD databases
# SCREEN says PLS is within ±2 kilobases of an experimental or predicted TSS. But how they do distance is not documented.

dists <- distanceToNearest(prox_enh, epd_tss, ignore.strand = TRUE, select = "all")
# lets plot a histogram of distances between prox_enh and nearest genes
p1 <- ggplot(as.data.frame(mcols(dists)), aes(x = distance + 1)) +
  geom_density() +
  scale_x_log10(labels = function(x) {
    format(x, scientific = FALSE)
  }, breaks = c(1000, 2000, 3000, 10000, 20000, 30000)) +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  ggtitle("Before filtering")
# remove all PLS with distance more than 2kb from the nearerst EPD-TSS. Based on the histogram
dists <- dists[mcols(dists)$distance < 2000]
p2 <- ggplot(as.data.frame(mcols(dists)), aes(x = distance + 1)) +
  geom_density() +
  scale_x_log10(labels = function(x) {
    format(x, scientific = FALSE)
  }, breaks = c(1000, 2000, 3000, 10000, 20000, 30000)) +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  ggtitle("After filtering")
ggsave(plot = cowplot::plot_grid(p1, p2, nrow = 2), file = "./plots/distance_distribtion_between_pELS_and_TSS.pdf")
prox_enh <- prox_enh[from(dists)]
new_metadata_columns <- mcols(epd_tss[to(dists)])
mcols(prox_enh) <- new_metadata_columns

# Step 3. Expand EPD promoters 50bp upstream and 10bp downstream.
epd_prom <- promoters(epd_tss, upstream = 50, downstream = 10)

# Step 4. Expand CAGE Peak to a minimum of 100 bases
cage_prom <- resize(cage_tss, width = pmax(100, width(cage_tss)), fix = "center")

# Keep features related to genes common across all datasets
common_genes <- purrr::reduce(.f = intersect, .x = list(rep_trans$ENSEMBL, epd_prom$ENSEMBL, prox_enh$ENSEMBL, prox_prom$ENSEMBL, cage_prom$ENSEMBL))

# Select a set of genes to use at random and (optionally) add the list of custom genes to it (if they are in common_genes)
# Also subset to smaller set of genes
set.seed(11235)
common_genes_to_test <- sample(common_genes, size = 2000L, replace = FALSE)
available_in_requested_genes <- intersect(custom_gene_list, common_genes)
common_genes_to_test <- union(common_genes_to_test, available_in_requested_genes)

rep_trans <- rep_trans[rep_trans$ENSEMBL %in% common_genes_to_test]
rep_trans_l <- rep_trans |>
  split(f = factor(rep_trans$ENSEMBL))
cage_prom <- cage_prom[cage_prom$ENSEMBL %in% common_genes_to_test]
cage_prom_l <- cage_prom |>
  split(f = cage_prom$ENSEMBL)
epd_prom <- epd_prom[epd_prom$ENSEMBL %in% common_genes_to_test]
epd_prom_l <- epd_prom |>
  split(f = epd_prom$ENSEMBL)
prox_enh <- prox_enh[prox_enh$ENSEMBL %in% common_genes_to_test]
prox_enh_l <- prox_enh |>
  split(f = prox_enh$ENSEMBL)
prox_prom <- prox_prom[prox_prom$ENSEMBL %in% common_genes_to_test]
prox_prom_l <- prox_prom |>
  split(f = prox_prom$ENSEMBL)

# Make sure the list names are identical
all_identical <- function(...) {
  x <- list(...)
  res <- c()
  for (i in seq.int(from = 2, to = length(x))) {
    res[i - 1] <- all(identical(x[1], x[i]))
  }
  return(all(res))
}
res <- all_identical(names(prox_prom_l), names(prox_enh_l), names(epd_prom_l), names(rep_trans_l))
if (res) {
  message("Conversion of all GRanges to GRangesList successfull")
} else {
  stop("The gene names across the lists are not identical")
}

dir.create("./gene_seggregated_enhs_proms_trans")
saveRDS(rep_trans_l, file = "./gene_seggregated_enhs_proms_trans/representative_transcripts_grangeslist.rds")
saveRDS(epd_prom_l, file = "./gene_seggregated_enhs_proms_trans/epd_promoters_grangeslist.rds")
saveRDS(cage_prom_l, file = "./gene_seggregated_enhs_proms_trans/cage_promoters_grangeslist.rds")
saveRDS(prox_prom_l, file = "./gene_seggregated_enhs_proms_trans/proximal_promoters_grangeslist.rds")
saveRDS(prox_enh_l, file = "./gene_seggregated_enhs_proms_trans/proximal_enhancers_grangeslist.rds")
