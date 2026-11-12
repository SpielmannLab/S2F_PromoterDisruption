# Readme

This repository contains the following files:

- The scripts used to generate the data and figures presented in our [manuscript](https://doi.org/10.64898/2026.11.06.764851)
- Lists of promoters/genes that were analysed in the manuscript and the predicted scores for the three promoter disruptions.

This is the general procedure:

1. First get representative transcripts, ENCODE promoters, and EPD core promoters using the scripts in utils/*
2. To link promoters to genes: `scripts/link_promoters_to_genes.R`
3. Define the promoter elements and region: `scripts/get_full_promoters_with_regions.R`
4. Define the region for CAGE scoring: `scripts/get_cage_scoring_regions.R`
5. Get the reference sequences for all the promoter regions: `scripts/get_ref_sequences_from_all_regions.R`
6. Create the ALT sequence and a table of promoters with variants to test: `scripts/create_full_promoter_variants_dinuc_shuffle.R` `scripts/create_full_promoter_variants_GCadj_random.R` `scripts/create_full_promoter_variants_N_repl.R`
7. Score the effect of the ALT by submission as variants to alpha genome: `scripts/score_cage_n_rna_seq_alphagenome.py`
8. Some plotting scripts
