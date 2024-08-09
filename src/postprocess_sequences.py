# Functions to process OTU sequences
import os
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import qiime2 as q2
from qiime2.plugins import diversity, feature_table

plt.rcParams.update({"font.family": "DejaVu Sans"})
plt.style.use("tableau-colorblind10")


def get_merged_ft_frequency(paths: list) -> pd.DataFrame:
    """
    This function loads FeatureTable[Frequency] from paths, merges them into a
    table, and converts the result to a pandas DataFrame.

    :param paths: List of paths to FeatureTable files
    :return: A pandas DataFrame of the merged sequence table
    """
    try:
        frequency_sequences = [q2.Artifact.load(path) for path in paths]
        (merged_frequency_table,) = feature_table.actions.merge(frequency_sequences)
        merged_frequency_df = merged_frequency_table.view(pd.DataFrame)

    except Exception as e:
        raise (f"Error in loading or merging frequency tables: {e}")

    return merged_frequency_df


def get_merged_sequences(paths: list) -> q2.Artifact:
    """
    This function loads FeatureData[Sequence] from paths and merges them into one
    table.

    :param paths: List of paths to sequence files
    :return: A Q2 artifact of Q2 type FeatureData[Sequence]
    """
    try:
        sequences = [q2.Artifact.load(path) for path in paths]
        (merged_table,) = feature_table.actions.merge_seqs(sequences)

    except Exception as e:
        raise (f"Error in loading or merging sequences: {e}")

    return merged_table


def remove_shallow_samples(freq_df, min_freq=1000):
    # Removes samples that have very shallow samples (samples with low seq per
    # sample) Q2 forum advise from JWD: >= 1000 seqs per sample as good
    # threshold to use
    # https://forum.qiime2.org/t/high-sample-size-heterogeneity-and-coda-analysis/24036/2

    check_msg = "Shape before filtering by shallow samples:"
    # transform to Q2 artifact
    freq_table = q2.Artifact.import_data("FeatureTable[Frequency]", freq_df)
    print(check_msg + f" {freq_table.view(pd.DataFrame).shape}")
    # funfact: 1 unique host (256 -> 255) is lost due to this filtering
    (freq_f_table,) = feature_table.actions.filter_samples(
        freq_table, min_frequency=min_freq
    )

    freq_f_df = freq_f_table.view(pd.DataFrame)
    print(check_msg.replace("before", "after") + f" {freq_f_df.shape}")

    # count no columns with all zeros in freq_df
    assert (freq_f_df == 0).all(axis=0).sum() == 0
    # check no more rows with all zeros
    assert (freq_f_df.sum(axis=1) == 0).sum() == 0

    return freq_f_df


def get_relative_abundance(
    ft: pd.DataFrame, feature_prefix: str, no_features: list = []
) -> pd.DataFrame:
    """
    Transform feature table from absolute to relative abundance. Only columns in
    feature_prefix are transformed. If feature_prefix is not set, then all
    features except no_features are transformed.
    """
    if feature_prefix:
        ft_cols = [x for x in ft.columns if x.startswith(feature_prefix)]
    elif len(no_features) > 0:
        ft_cols = [x for x in ft.columns if x not in no_features]
    else:
        ft_cols = ft.columns.tolist()
    ft_rel = ft.copy()
    ft_rel[ft_cols] = ft_rel[ft_cols].div(ft_rel[ft_cols].sum(axis=1), axis=0)

    # round needed as certain 1.0 are represented in different digits 2e-16
    assert ft_rel[ft_cols].sum(axis=1).round(5).eq(1.0).all()

    return ft_rel


def _load_silva_phylo_tree(path_to_get, n_threads):
    path2tree = os.path.join(
        path_to_get, "silva-138-1-ssu-nr99-seqs-rooted-tree-proc-ben.qza"
    )
    if not os.path.isfile(path2tree):
        command = f"../src/build_silva_full_tree.sh {path_to_get} {n_threads}"
        subprocess.run(command, shell=True)

    return q2.Artifact.load(path2tree)


def bootstrapped_alpha_div(
    div_metric, freq_df, rarefaction_depth, path_to_phylo, n_threads, n=100
) -> pd.Series:
    df_alpha = pd.DataFrame(index=freq_df.index)
    freq_art = q2.Artifact.import_data("FeatureTable[Frequency]", freq_df)

    for i in range(0, n):
        # todo: add seed once it's available
        (rarefied_table,) = feature_table.actions.rarefy(
            table=freq_art, sampling_depth=rarefaction_depth
        )

        # calculating alpha diversity
        if div_metric == "faith_pd":
            tree_art = _load_silva_phylo_tree(path_to_phylo, n_threads)
            (i_alpha_tab,) = diversity.actions.alpha_phylogenetic(
                table=rarefied_table, metric=div_metric, phylogeny=tree_art
            )
        else:
            (i_alpha_tab,) = diversity.actions.alpha(
                table=rarefied_table, metric=div_metric
            )
        i_alpha_ser = i_alpha_tab.view(pd.Series)
        i_alpha_ser.name = f"{div_metric}_{i}"
        df_alpha = df_alpha.merge(
            i_alpha_ser.to_frame(), left_index=True, right_index=True
        )

    # calculate mean over all iterations
    alpha_all = df_alpha.mean(axis=1)
    alpha_all.name = f"div_alpha_{div_metric}"
    return alpha_all
