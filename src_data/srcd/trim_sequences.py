"""
Module to find trim sequences of each study & display demux summary
to choose truncation length
"""

import argparse
import glob
import os
import re

import pandas as pd
import qiime2 as q2
from qiime2.plugins import cutadapt, demux
from srcd import meta_processor as mproc


def parse_arguments():
    """Parse command line arguments provided"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--path2md", type=str, required=True)
    parser.add_argument("--path2seq", type=str, required=True)
    parser.add_argument("--threads", default=6, type=int)
    return parser.parse_args()


def get_study_reads(study, df_md, path2seq):
    """Get study sequence reads, library layout, primers of this study
    & study metadata"""
    # select metadata of this one study
    df_study = df_md[df_md["study_name"] == study].copy()
    bioproject_id = df_study["bioproject_id"].unique().tolist()[0]

    # define library layout to select correct sequences
    type_of_reads = df_study["library_layout"].unique().tolist()
    assert len(type_of_reads) == 1
    type_of_reads = type_of_reads[0]

    # get only sequences that belong to this study
    if type_of_reads == "SINGLE":
        study_reads = q2.Artifact.load(
            os.path.join(path2seq, bioproject_id, "single_reads.qza")
        )
    else:
        study_reads = q2.Artifact.load(
            os.path.join(path2seq, bioproject_id, "paired_reads.qza")
        )

    # identify & trim primers
    primers = df_study["exp_primer"].unique()
    assert len(primers) == 1
    primers = primers[0]

    return study_reads, type_of_reads, primers, df_study


def summarize_reads(study_reads, path2save):
    """Summarize reads in study_reads and save in file path2save"""
    (sum_reads,) = demux.actions.summarize(data=study_reads)
    sum_reads.save(path2save)
    print(f"Saved demux summary in: {path2save}")


def trim_sequences(path2md, path2seq, threads):
    """Trims sequences in path2seq of all studies included in path2md"""
    # read metadata of all studies
    df_md = pd.read_csv(path2md, sep="\t", index_col=0, dtype="str")

    for study in df_md["study_name"].unique():
        print(f"Trimming: {study}...")
        path2save_seq = os.path.join(path2seq, f"trimmed_{study}*.qza")

        if glob.glob(path2save_seq):
            print(f"Trimmed sequences already exist: {path2save_seq}.")
        else:
            study_reads, type_of_reads, primers, df_study = get_study_reads(
                study, df_md, path2seq
            )

            if type_of_reads == "PAIRED":
                fwd_primer, rev_primer = tuple(re.findall(r"\[(.*?)\]", primers))
                # great Q2forum post on parameters and explanation:
                # https://forum.qiime2.org/t/cutadapt-adapter-vs-front/15450/2?u=adamova
                # (5’) 515primer->fwd-sequence->index (3’) …
                # (3’) index<-rev-sequence<-806primer (5’)
                # front_f = 515, front_r = 806
                (study_trim_reads,) = cutadapt.actions.trim_paired(
                    demultiplexed_sequences=study_reads,
                    cores=threads,
                    front_f=[fwd_primer],  # 515f: search in fwd read
                    front_r=[rev_primer],  # 806r: search in rev read
                    match_read_wildcards=True,
                    match_adapter_wildcards=True,
                )
            elif type_of_reads == "SINGLE":
                primer = tuple(re.findall(r"\[(.*?)\]", primers))[0]
                (study_trim_reads,) = cutadapt.actions.trim_single(
                    demultiplexed_sequences=study_reads,
                    cores=threads,
                    front=[primer],
                    match_read_wildcards=True,
                    match_adapter_wildcards=True,
                )
            # Save trimmed reads & demux summary for manual extraction of
            # trunc length for dada 2. Split by study subcohort if it
            # exists - since trunc_len are
            # not always identical across subcohorts.
            unique_subcohorts = df_study["study_subcohort"].unique()
            # only 1 subcohort:
            if len(unique_subcohorts) == 1:
                # save reads
                path2save_seq = path2save_seq.replace("*", "")
                study_trim_reads.save(path2save_seq)
                print(f"Saved trimmed sequences: {path2save_seq}")
                # create and save summary
                path4summary = path2save_seq.replace(".qza", "_summary.qzv")
                summarize_reads(study_trim_reads, path4summary)
            else:
                # multiple subcohorts available - separately summarize to
                # find accurate trim_length:
                for cohort in unique_subcohorts:
                    # get reads of study subcohort only
                    md_cohort = df_study[df_study["study_subcohort"] == cohort].copy()
                    (seq_cohort,) = demux.actions.filter_samples(
                        demux=study_trim_reads,
                        metadata=mproc.transform_to_q2metadata(md_cohort),
                    )
                    # save and summarize
                    path4cohort = path2save_seq.replace("*.qza", f"_{cohort}.qza")
                    seq_cohort.save(path4cohort)
                    print(
                        f"Saved trimmed sequences for subcohort {cohort}: {path4cohort}"
                    )
                    path4summary = path4cohort.replace(".qza", "_summary.qzv")
                    summarize_reads(seq_cohort, path4summary)


if __name__ == "__main__":
    args = parse_arguments()
    print(
        f"Trimming sequences saved in {args.path2seq} with metadata in {args.path2md}"
    )
    trim_sequences(path2md=args.path2md, path2seq=args.path2seq, threads=args.threads)
