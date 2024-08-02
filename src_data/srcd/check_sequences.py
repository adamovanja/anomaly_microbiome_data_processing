"""
Module to get summary statistics of sequences of each study
after each processing step
"""

import argparse
import os
import re
import warnings
from statistics import mean

import numpy as np
import pandas as pd
import qiime2 as q2
from qiime2.plugins import demux
from skbio import DNA
from srcd.utils import extract_file

# ignoring FutureWarning for demux summarize action
warnings.simplefilter(action="ignore", category=FutureWarning)


def parse_arguments():
    """Parse command line arguments provided"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--path2md", type=str, required=True)
    return parser.parse_args()


def _get_nt_float(str):
    """Extract digit before 'nts' string as float"""
    str_nb = re.search(r"(\d*) nts", str).group(1)
    return float(str_nb)


def _extract_counts_sample(viz, read_type):
    """Extract counts from demux summarize visualisation artifact"""
    # get table with number of sequences
    nb_seqs = extract_file(viz, "per-sample-fastq-counts.tsv")
    summary_seqs = nb_seqs.describe()
    # extract
    # ... sample count & seq nb mean: averaged between reverse and forward reads
    sample_count = summary_seqs.loc[["count"], :].mean(axis=1).values[0]
    mean_nb_seq_p_sample = round(
        summary_seqs.loc[["mean"], :].mean(axis=1).values[0], 6
    )

    # # get total sequences & assert above calculated mean
    # # ! these are not unique sequences - hence not comparable to counts in
    # # ! denoise and cluster steps
    # overview_tab = extract_file(viz, "overview.html")
    # total_nb_sequences = overview_tab[0].iloc[4, 1:].values.mean()
    # assert mean_nb_seq_p_sample == overview_tab[0].iloc[2, 1:].values.mean()

    # get tables with length of reads
    if read_type == "paired":
        fwd_len, rev_len = extract_file(viz, "quality-plot.html")
        # ... length of reads: averaged between reverse and forward reads
        fwd_nts = _get_nt_float(fwd_len.iloc[4, 1])
        rev_nts = _get_nt_float(rev_len.iloc[4, 1])
        median_len_sequences = mean([fwd_nts, rev_nts])
    else:
        fwd_len = extract_file(viz, "quality-plot.html")[0]
        # ... length of reads: averaged between reverse and forward reads
        median_len_sequences = _get_nt_float(fwd_len.iloc[4, 1])
    return (
        sample_count,
        np.NaN,
        mean_nb_seq_p_sample,
        np.NaN,
        median_len_sequences,
    )


def _extract_counts_feature(repseq, ft):
    """
    Extract counts of a FeatureData[Sequence] and a FeatureTable[Frequency]
    artifact. Counts to be retrieved include: sample count, total number of sequences,
    mean number of sequences per sample, mean number of unique sequences per sample and
    median length of sequences.
    """

    # process feature table:
    ft_df = ft.view(pd.DataFrame)

    sample_count = float(ft_df.shape[0])
    total_nb_sequences = float(ft_df.shape[1])

    # this count does NOT count unique sequences
    mean_nb_seq_p_sample = ft_df.sum(axis=1).mean()

    # count only unique sequences
    mean_nb_unique_seq_p_sample = ft_df.apply(
        lambda row: (row != 0.0).sum(), axis=1
    ).mean()

    # process representative sequences
    seq_ser = repseq.view(pd.Series)
    seq_lengths = seq_ser.apply(lambda seq: len(DNA(seq)))
    median_len_sequences = seq_lengths.median()

    return (
        sample_count,
        total_nb_sequences,
        mean_nb_seq_p_sample,
        mean_nb_unique_seq_p_sample,
        median_len_sequences,
    )


def check_sequences_per_study(study_name, study_id, read_type, md_df, path2data):
    """
    Get summary statistics of sequences of study after each processing step
    """
    print(f"Summarising {study_name} ...")
    # initialise summary dataframe
    ls_steps = ["raw", "trim", "denoise", "cluster"]
    ls_cols = [
        "sample_count",
        "total_nb_sequences",
        "mean_nb_seq_p_sample",
        "mean_nb_unique_seq_p_sample",
        "median_len_sequences",
    ]
    df_all = pd.DataFrame(index=ls_steps, columns=ls_cols)

    # ! Step1) raw reads
    path2raw_seq = os.path.join(path2data, study_id, f"{read_type.lower()}_reads.qza")
    raw_reads = q2.Artifact.load(path2raw_seq)
    # filter by study_cohort_name
    (filt_raw_reads,) = demux.actions.filter_samples(
        raw_reads, q2.Metadata(md_df), where=f"study_cohort_name='{study_name}'"
    )
    (sum_raw,) = demux.actions.summarize(data=filt_raw_reads)
    df_all.loc["raw", :] = _extract_counts_sample(sum_raw, read_type)

    # ! Step2) trimmed reads
    path2trim_seq = os.path.join(path2data, f"trimmed_{study_name}_summary.qzv")
    sum_trim = q2.Visualization.load(path2trim_seq)
    df_all.loc["trim", :] = _extract_counts_sample(sum_trim, read_type)

    # ! Step3) denoised reads
    # # ... get FeatureData[Sequences]
    asv_rep_seqs = q2.Artifact.load(os.path.join(path2data, f"repseq_{study_name}.qza"))
    # ... get FeatureTable[Frequency]
    asv = q2.Artifact.load(os.path.join(path2data, f"asv_{study_name}.qza"))
    # ... get counts:
    df_all.loc["denoise", :] = _extract_counts_feature(asv_rep_seqs, asv)

    # ! Step4) clustered reads
    # ... get FeatureData[Sequences]
    otu_rep_seqs = q2.Artifact.load(
        os.path.join(path2data, f"otu_seq_{study_name}.qza")
    )
    # ... get FeatureTable[Frequency]
    otu = q2.Artifact.load(os.path.join(path2data, f"otu_table_{study_name}.qza"))
    # ... get counts:
    df_all.loc["cluster", :] = _extract_counts_feature(otu_rep_seqs, otu)

    return df_all


def check_sequences(path2md):
    """
    Run check_sequences_per_study for each study present in path2md file
    """
    # read metadata to get study_name, study_id and library values
    md_df = pd.read_csv(path2md, sep="\t", index_col=0, dtype="str")
    ls_study_pid = (
        md_df[["study_cohort_name", "bioproject_id", "library_layout"]]
        .drop_duplicates()
        .values.tolist()
    )
    path2data = os.path.dirname(path2md)
    path2save = os.path.join(path2data, "check_seqs")
    if not os.path.exists(path2save):
        os.makedirs(path2save)

    # check sequences
    for study_name, study_id, read_type in ls_study_pid:
        path2out = os.path.join(path2save, f"stats_{study_name}.csv")

        if not os.path.isfile(path2out):
            df_stats = check_sequences_per_study(
                study_name, study_id, read_type, md_df, path2data
            )
            df_stats.to_csv(path2out)


if __name__ == "__main__":
    args = parse_arguments()
    print(f"Checking sequences of studies saved in {args.path2md}")
    check_sequences(path2md=args.path2md)
