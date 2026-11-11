#!/usr/bin/env Rscript
library(rtracklayer) # Enables importing bedfiles etc

# Exports a GRanges object of the following format:
# GRanges object with 296996 ranges and 1 metadata column:
#            seqnames            ranges strand |        cCRE_class
#               <Rle>         <IRanges>  <Rle> |       <character>
#        [1]     chr1       10034-10250      * | Proximal enhancer

bb_file <- BigBedFile("../input/encodeCcreRegistry.bb") # Dowloaded from https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode4/ccre/encodeCcreRegistry.bb")
ccres <- import.bb(bb_file, selection = BigBedSelection(bb_file, colnames = "cCRE_class"))
encode_promoters <- ccres[ccres$cCRE_class %in% c("Promoter", "Proximal enhancer")]
saveRDS(encode_promoters, "./input/encode_promoters.rds")

