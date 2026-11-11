# Readme

This repository contains the following files:

- The scripts used to generate the data and figures presented in our [manuscript](https://doi.org/10.64898/2026.11.06.764851)
- Lists of promoters/genes that were analysed in the manuscript.

This is the general procedure:

1. First get representative transcripts, ENCODE promoters, and EPD core promoters using the scripts in utils/*.
2. To link promoters to genes: `link_promoters_to_genes.R`
3. Define the promoter elements and region: `get_full_promoters_with_regions.R`
4. Define the region for CAGE scoring: `get_cage_scoring_regions.R`
5. Get the reference sequences for all the promoter regions: `get_ref_sequences_from_all_regions.R`
6. Create the ALT sequence: `create_full_promoter_variants.R`
7. Score the effect of the ALT by submission as variants to alpha genome: `score_....py`. In this script, there are two version - one where the scoring region is shifted by 200000bp so, the promoter actually is quite far to the edge in the 500kbp scoring window. The other one is where the scoring region is shifted by 50kbp , so the promoter is closer to the center of the 500kbp scoring window. There is actually a third one, where the scoring region is shifted by 200kbp and then the scores are not calculated per exon, but for the entire gene body. Use the more centered verions.
8. Some plotting scripts
