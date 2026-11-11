import pandas as pd
from alphagenome.data import genome
import numpy as np

def score_mean_rna_seq_per_interval(rna_track_data, intervals_list: pd.DataFrame):
    """Given track data predicted by alphagenome and a list containing alphagenome intervals, calculate the mean scores at each interval. Only (parts of) intervals that has been scored will be used. Returns a comma separated string for each interval in the same order. Empty string for intervals not scored.
    """
    # First sort the intervals taking into accoung the strand direction
    scored_interval = rna_track_data.interval
    scored_intervals = [
        scored_interval.intersect(interval) for interval in intervals_list
    ]
    strand_of_interest = set([interval.strand for interval in intervals_list])

    if len(strand_of_interest) > 1:
        raise ValueError("Strands in the intervals_list are not all the same")

    if strand_of_interest.pop() == "+":
        rna_track_strand_data = rna_track_data.filter_to_positive_strand()
    else:
        rna_track_strand_data = rna_track_data.filter_to_negative_strand()

    mean_interval_scores = [
        # (rna_track_strand_data.slice_by_interval(scored_interval).values).sum(axis = 0)
        np.concatenate(rna_track_strand_data.slice_by_interval(scored_interval).values).mean(axis = 0)
        if scored_interval is not None
        else np.nan
        for scored_interval in scored_intervals
    ]
    mean_score_per_interval = ",".join([np.array2string(s) if not np.isnan(s) else '' for s in mean_interval_scores])

    return mean_score_per_interval

def score_mean_rna_seq_across_intervals(rna_track_data, intervals_list: pd.DataFrame):
    """Given track data predicted by alphagenome and a list containing alphagenome intervals, provides the mean scores across all intervals. Only (parts of) intervals that were scored will be used.
    """
    scored_interval = rna_track_data.interval
    scored_intervals = [
        scored_interval.intersect(interval) for interval in intervals_list
    ]
    strand_of_interest = set([interval.strand for interval in intervals_list])

    if len(strand_of_interest) > 1:
        raise ValueError("Strands in the intervals_list are not all the same")

    if strand_of_interest.pop() == "+":
        rna_track_strand_data = rna_track_data.filter_to_positive_strand()
    else:
        rna_track_strand_data = rna_track_data.filter_to_negative_strand()

    interval_scores = [
        # (rna_track_strand_data.slice_by_interval(scored_interval).values).sum(axis = 0)
        (rna_track_strand_data.slice_by_interval(scored_interval).values)
        for scored_interval in scored_intervals
        if scored_interval is not None
    ]
    interval_scores_mean = np.concatenate(interval_scores, axis = 0).mean(axis = 0)

    return interval_scores_mean

def extract_exon_intervals_list_per_gene_from_exons_df(
    exons_df: pd.DataFrame, ensembl: str
):
    """Provided a pandas dataframe containing an exon per entry, with the following columns: `seqnames`, `start`, `end`, `strand`, `ENSEMBL`
    extract only those for a given `ensembl` id. Returns a list containing alphagenome intervals.
    """
    exons_df_subset = exons_df[exons_df["ENSEMBL"] == ensembl]
    exon_intervals = [
        genome.Interval(
            chromosome=exon["seqnames"],
            start=exon["start"],
            end=exon["end"],
            strand=exon["strand"],
        )
        for idx, exon in exons_df_subset.iterrows()
    ]
    return exon_intervals
