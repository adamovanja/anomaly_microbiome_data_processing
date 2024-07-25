"""
Module to cluster sequences of each study with closed-reference
clustering
"""
import argparse
import os
import subprocess

import pandas as pd
import qiime2 as q2
from qiime2.plugins import vsearch
from srcd.denoise_sequences import save_artifact


def parse_arguments():
    """parse CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--path2md", type=str, required=True)
    parser.add_argument("--path2seq", type=str, required=True)
    parser.add_argument("--threads", default=6, type=int)
    return parser.parse_args()


def load_reference_seqs(path2ref):
    """Loads or downloads reference sequences"""
    if not os.path.isfile(path2ref):
        url = "https://data.qiime2.org/2024.5/common/silva-138-99-seqs-515-806.qza"
        cmd_wget = ["wget", "-nv", "-O", path2ref, url]
        subprocess.run(cmd_wget)

    return q2.Artifact.load(path2ref)


def cluster_sequences_one_study(study_cohort, path2seq, threads):
    """Cluster sequences of one study"""
    print(f"Clustering: {study_cohort}...")
    path_otu_table = os.path.join(path2seq, f"otu_table_{study_cohort}")

    if os.path.isfile(path_otu_table):
        print(f"Clustered otu table already exists in: {path_otu_table}.")
    else:
        # read ASV sequence table
        asv = q2.Artifact.load(os.path.join(path2seq, f"asv_{study_cohort}.qza"))
        # read representative sequences table
        repseq = q2.Artifact.load(os.path.join(path2seq, f"repseq_{study_cohort}.qza"))
        # read reference sequences
        ref = load_reference_seqs(
            os.path.join(path2seq, "silva_138_99_seqs_515_806.qza")
        )

        # closed-reference clustering
        (
            otu_table,
            otu_seq,
            unmatched_seq,
        ) = vsearch.actions.cluster_features_closed_reference(
            sequences=repseq,
            table=asv,
            reference_sequences=ref,
            perc_identity=0.97,
            threads=threads,
        )
        # save outputs
        save_artifact(otu_table, path_otu_table, "otu table")

        save_artifact(
            otu_seq,
            os.path.join(path2seq, f"otu_seq_{study_cohort}.qza"),
            "otu sequences",
        )
        save_artifact(
            unmatched_seq,
            os.path.join(path2seq, f"otu_umn_{study_cohort}.qza"),
            "unmatched sequences",
        )


def cluster_sequences(path2md, path2seq, threads):
    """Cluster denoised sequences saved in path2seq of all study
    subcohorts included in path2md"""
    # read metadata of all studies
    df_md = pd.read_csv(path2md, sep="\t", index_col=0, dtype="str")

    for study_cohort in df_md["study_cohort_name"].unique():
        cluster_sequences_one_study(study_cohort, path2seq, threads)


if __name__ == "__main__":
    args = parse_arguments()
    print(f"Clustering sequences of studies present in metadata {args.path2md}")

    cluster_sequences(
        path2md=args.path2md,
        path2seq=args.path2seq,
        threads=args.threads,
    )
