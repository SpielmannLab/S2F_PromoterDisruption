#!/usr/bin/env python
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from alphagenome.data import genome, track_data, transcript
from alphagenome.models import dna_client, variant_scorers
from alphagenome_research.model import dna_model
from dotenv import load_dotenv
from huggingface_hub import login
from tqdm import tqdm

# Add custom functions for scoring
sys.path.append("./scripts_n_settings/")
import my_helper_functions
print("Now loaded all the libraries")

# Log in to huggingface (this works, because we set the MY_TOKEN_PATH env variable separately in the shell)
load_dotenv(dotenv_path=Path(os.environ["MY_TOKEN_PATH"]))
login()
print("Now logged in to HuggingFace")

# Set up a dna_model (which is a server)
model = dna_model.create_from_huggingface("all_folds")
print("Started the dna_model")

# Choose tissues to score
available_ontologies_shared = pd.read_csv(
    "../input/alphagenome_available_ontologies_shared.tsv", sep="\t"
)
selected_ontologies = (
    available_ontologies_shared[
        (available_ontologies_shared.biosample_type == "tissue")
        & (
            ~available_ontologies_shared["Assay title_rnaseq"].str.contains(
                "total RNA-seq"
            )
        )
    ]
    .copy()
    .reset_index(drop=True)
)

# Gather arguments
ontology_required = set(selected_ontologies["ontology_curie"])
# ontology_required = "UBERON:0000955" # For testing just the brain

# Load the variants
vcf = pd.read_csv("prom_variants_to_test.tsv", sep="\t")
exons_df = pd.read_csv("../input/representative_exons.tsv", sep="\t")

# Parse organism specification.
organism = "human"
organism_map = {
    "human": dna_model.Organism.HOMO_SAPIENS,
    "mouse": dna_model.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]

# Choose built in variant scorers
wanted_scorers = [
    "RNA_SEQ"
]
# Do not use the CAGE scorer - not great! We will use custom scorer
# RNA_SEQ scorer is also not great, but lets keep it
all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [all_scorers[key] for key in wanted_scorers]

# Score variants in the VCF file. For visualisation, we will use a separate script.
scores_cage = []
scores_rna_seq = []
variants = []
cage_scoring_intervals = []
rna_seq_scoring_intervals = []

for i, vcf_row in tqdm(vcf.iterrows()):

    # Remove all where the promter region is too large
    if len(vcf_row.ref_seq) > 10000:
        continue

    variant = genome.Variant(
        chromosome=str(vcf_row.seqnames),
        position=int(vcf_row.start),
        reference_bases=vcf_row.ref_seq,
        alternate_bases=vcf_row.alt_seq,
        name=vcf_row.variant_id,
    )

    cage_interval = genome.Interval(
        chromosome=str(vcf_row.seqnames),
        strand=vcf_row.cage_score_strand,
        start=vcf_row.cage_score_start,
        end=vcf_row.cage_score_end,
    )

    rna_seq_interval = genome.Interval(
        chromosome=str(vcf_row.seqnames),
        strand=vcf_row.rna_score_strand,
        start=vcf_row.rna_score_start,
        end=vcf_row.rna_score_end,
    )

    # Predict the track - for scoring ref and alt alleles in CAGE and RNA.
    # Supported are 16KB, 100KB, 500KB and 1MB - Lets keep 500kb.
    interval_to_predict = cage_interval.resize(dna_client.SEQUENCE_LENGTH_500KB).shift(
        50000
    )
    cage_scoring_interval = cage_interval.intersect(interval_to_predict)
    rna_seq_scoring_interval = rna_seq_interval.intersect(interval_to_predict)

    track = model.predict_variant(
        interval=interval_to_predict,
        variant=variant,
        ontology_terms=ontology_required,
        requested_outputs=[dna_model.OutputType.CAGE,
                           dna_model.OutputType.RNA_SEQ],
    )

    # ---- Add my scoring diff-log-sum/mean around at the cage/rna scoring-region ------ #
    # extract the relevant individual tracks
    if cage_interval.strand == "+":
        ref_rna_seq = (
            track.reference.rna_seq.filter_to_positive_strand().slice_by_interval(
                interval=rna_seq_scoring_interval
            )
        )
        alt_rna_seq = (
            track.alternate.rna_seq.filter_to_positive_strand().slice_by_interval(
                interval=rna_seq_scoring_interval
            )
        )
        ref_cage = track.reference.cage.filter_to_positive_strand().slice_by_interval(
            interval=cage_scoring_interval
        )
        alt_cage = track.alternate.cage.filter_to_positive_strand().slice_by_interval(
            interval=cage_scoring_interval
        )
    else:
        ref_rna_seq = (
            track.reference.rna_seq.filter_to_negative_strand().slice_by_interval(
                interval=rna_seq_scoring_interval
            )
        )
        alt_rna_seq = (
            track.alternate.rna_seq.filter_to_negative_strand().slice_by_interval(
                interval=rna_seq_scoring_interval
            )
        )
        ref_cage = track.reference.cage.filter_to_negative_strand().slice_by_interval(
            interval=cage_scoring_interval
        )
        alt_cage = track.alternate.cage.filter_to_negative_strand().slice_by_interval(
            interval=cage_scoring_interval
        )

    # Calculate the mean across the entire transcript (at scored region)
    score_rna_seq_df = ref_rna_seq.metadata.copy().reset_index()
    score_rna_seq_df["variant_id"] = vcf_row.variant_id
    score_rna_seq_df["ENSEMBL"] = vcf_row.ENSEMBL
    # mean of rna-track across the entire scored region. Not ideal
    score_rna_seq_df["ref_mean_rna_seq"] = ref_rna_seq.values.mean(axis=0)
    score_rna_seq_df["alt_mean_rna_seq"] = alt_rna_seq.values.mean(axis=0)
    score_rna_seq_df["diff_log_mean_rna_seq"] = np.log2(
        (score_rna_seq_df["alt_mean_rna_seq"] + 1)
        / (score_rna_seq_df["ref_mean_rna_seq"] + 1)
    )

    # Calculate the sum centered around tss from positive or negative strands
    score_cage_df = ref_cage.metadata.copy().reset_index()
    score_cage_df["variant_id"] = vcf_row.variant_id
    score_cage_df["ENSEMBL"] = vcf_row.ENSEMBL
    score_cage_df["ref_sum_cage"] = ref_cage.values.sum(axis=0)
    score_cage_df["alt_sum_cage"] = alt_cage.values.sum(axis=0)
    score_cage_df["diff_log_sum_cage"] = np.log2(
        (score_cage_df["alt_sum_cage"] + 1) /
        (score_cage_df["ref_sum_cage"] + 1)
    )

    ### Now calculate rna_seq scores only at exonic regios. Note: this is not per exon, but mean across all exons
    exons_of_just_this_gene = (
        my_helper_functions.extract_exon_intervals_list_per_gene_from_exons_df(
            exons_df, vcf_row.ENSEMBL
        )
    )
    # Sort the exons, based on strand
    exons_of_just_this_gene.sort(reverse = (set([i.strand for i in exons_of_just_this_gene]) == set("-")))

    # Get mean across all exons
    score_rna_seq_df["ref_exon_mean_rna_seq"] = my_helper_functions.score_mean_rna_seq_across_intervals(
        ref_rna_seq, exons_of_just_this_gene
    )
    score_rna_seq_df["alt_exon_mean_rna_seq"] = my_helper_functions.score_mean_rna_seq_across_intervals(
        alt_rna_seq, exons_of_just_this_gene
    )
    score_rna_seq_df["diff_log_exon_mean_rna_seq"] = np.log2(
        (score_rna_seq_df["alt_exon_mean_rna_seq"] + 1)
        / (score_rna_seq_df["ref_exon_mean_rna_seq"] + 1)
    )

    # Get mean_per_exon. Also save exonic interval information
    score_rna_seq_df["ref_mean_per_exon_rna_seq"] = my_helper_functions.score_mean_rna_seq_per_interval(
        ref_rna_seq, exons_of_just_this_gene
    )
    score_rna_seq_df["alt_mean_per_exon_rna_seq"] = my_helper_functions.score_mean_rna_seq_per_interval(
        alt_rna_seq, exons_of_just_this_gene
    )
    score_rna_seq_df["exon_starts"] = ",".join([str(i.start) for i in exons_of_just_this_gene])
    score_rna_seq_df["exon_ends"] = ",".join([str(i.end) for i in exons_of_just_this_gene])

    scores_rna_seq.append(score_rna_seq_df)
    scores_cage.append(score_cage_df)

    # also save the variants
    variants.append(variant)
    cage_scoring_intervals.append(cage_scoring_interval)
    rna_seq_scoring_intervals.append(rna_seq_scoring_interval)

scored_intervals_rna_seq_df = pd.DataFrame({
    "chr": [interval.chromosome for interval in rna_seq_scoring_intervals],
    "start": [interval.start for interval in rna_seq_scoring_intervals],
    "end": [interval.end for interval in rna_seq_scoring_intervals],
    "strand": [interval.strand for interval in rna_seq_scoring_intervals]
})
scored_intervals_cage_df = pd.DataFrame({
    "chr": [interval.chromosome for interval in cage_scoring_intervals],
    "start": [interval.start for interval in cage_scoring_intervals],
    "end": [interval.end for interval in cage_scoring_intervals],
})

# Now save the results either as tsv
pd.concat(scores_cage).to_csv("scores_cage.tsv", sep="\t")
pd.concat(scores_rna_seq).to_csv("scores_rna_seq.tsv", sep="\t")
scored_intervals_cage_df.to_csv("scored_intervals_cage.tsv", sep="\t")
scored_intervals_rna_seq_df.to_csv("scored_intervals_rna_seq.tsv", sep="\t")
