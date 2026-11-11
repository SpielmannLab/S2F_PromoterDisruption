#!/usr/bin/env Rscript

# Output is a R GenomicRanges object containing representaive transcripts (MANE, Refseq SELECT) with the following format saved as rds file
#                     seqnames   ranges        strand |   TXID        TXNAME           ENSEMBL
#                     <Rle>      <IRanges>     <Rle>  |   <integer>   <character>      <character>
#  NM_001005484.2     chr1       65019-65518      +   |   1           NM_001005484.2   ENSG00000186092

library(txdbmaker) # Create TxDb from UCSC gene tracks
# library(AnnotationHub) #
library(rtracklayer) # Enables importing bedfiles etc
library(GenomicRanges) # For interesction of regions
library(BSgenome)
library(biomaRt)

# Collect refseqselect/mane transcripts from UCSC #####
refseq_select <- makeTxDbFromUCSC(genome = "hg38", tablename = "ncbiRefSeqSelect") # This is RefSeq Select and MANE

# Remove non chromosomal scaffolds -> this works. (https://bioc.r-universe.dev/GenomicFeatures/doc/GenomicFeatures.html)
seqlevels(refseq_select) <- grep(seqlevels(refseq_select), invert = TRUE, pattern = "_", value = TRUE)

# Get the transcripts per gene
rep_transcripts <- transcripts(refseq_select, columns = c("TXNAME", "TXID"))

# Note this actually contains duplicates, because of genes in the ?pseudoautosomal region between X and Y chromosomes
# Optionally one can remove both copies here:
remove_duplicated_genes_in_XnY <- TRUE
if (remove_duplicated_genes_in_XnY) {
  dup_ids <- duplicated(rep_transcripts$TXNAME) & (as.vector(seqnames(rep_transcripts)) %in% c("chrY", "chrX"))
  dup_transcripts <- rep_transcripts$TXNAME[dup_ids]
  rep_transcripts <- rep_transcripts[!(rep_transcripts$TXNAME %in% dup_transcripts)]
}

# Use biomart to convert NCBI refseq transcript ids to ensembl genes
# Connect to Ensembl database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
tx_names <- gsub(rep_transcripts$TXNAME, pattern = "\\..*", replacement = "")

# Query the conversion table
biomart_results <- getBM(
  attributes = c("refseq_mrna", "ensembl_transcript_id", "external_gene_name", "chromosome_name", "ensembl_gene_id"),
  filters = "refseq_mrna",
  values = tx_names,
  mart = ensembl
)
biomart_results <- biomart_results %>% filter(!grepl(chromosome_name, pattern = "_"))
if (any(duplicated(biomart_results$refseq_mrna))) {
  warning("There are duplicated NCBI refseq transcript ids")
}
gene_ids <- biomart_results$ensembl_gene_id[match(tx_names, biomart_results$refseq_mrna)]
rep_transcripts$ENSEMBL <- gene_ids

# Remove NA
rep_transcripts <- rep_transcripts[!is.na(rep_transcripts$ENSEMBL)]
rep_transcripts$TXNAME <- gsub(rep_transcripts$TXNAME, pattern = "\\..*", replacement = "")

# This now contains promoter regions of Reseq select transcripts with ensmbl gene_ids - no duplicates
saveRDS(rep_transcripts, "./input/representative_transcripts.rds")
