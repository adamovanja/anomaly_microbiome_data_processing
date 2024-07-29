import argparse
import os
import warnings

import pandas as pd
from srcd.postprocess_sequences import (
    bootstrapped_alpha_div,
    get_merged_ft_frequency,
    get_merged_sequences,
    get_relative_abundance,
    remove_shallow_samples,
)
from srcd.utils import filter_md_by_ft
from srcd.utils_tax import perform_taxonomic_classification

warnings.filterwarnings("ignore", category=FutureWarning)


def filter_metadata(path_md):
    # read metadata
    md_df = pd.read_csv(path_md, sep="\t", index_col=0)
    md_df = md_df.loc[md_df["study_subcohort"].isin(["abx", "karelia", "t1d"])]

    # filter out hosts from t1d that had abx event (no detailed info available
    # these can only be used for noabx learning: 7 hosts = 140 samples)
    excl_t1d_samples = md_df.loc[
        ((md_df["study_subcohort"] == "t1d") & (md_df["abx_ever"] == True))
    ].index.tolist()
    md_df = md_df.loc[~md_df.index.isin(excl_t1d_samples)]

    md_df.sort_values(
        [
            "host_id",
            "age_days",
        ],
        ascending=[True, True],
        inplace=True,
    )
    # impute missing columns of t1d subcohort for consistency
    md_df.loc[md_df["study_subcohort"] == "t1d", "abx_any_cumcount"] = 0.0
    md_df.loc[md_df["study_subcohort"] == "t1d", "abx_max_count_ever"] = 0.0

    # select only md columns of interest to this usecase
    cols_to_keep = [
        "age_days",
        "host_id",
        # static covariate
        "delivery_mode",
        # optional static covariates
        "sex",
        "geo_location_name",
        # dynamic covariate
        "diet_milk",
        "diet_weaning",
        # abx abnormality indicators
        "age_months_rounded05",  # reference for all abx_ columns
        "age_months_rounded1",  # for plotting
        # abx any
        "abx_any_last_t_dmonths",
        "abx_any_last_dur_days",
        "abx_any_cumcount",
        "abx_any_cumdur_days",
        # todo: remove broad + narrow grouping here
        # abx broad
        "abx_broad_last_t_dmonths",
        "abx_broad_last_dur_days",
        "abx_broad_cumcount",
        "abx_broad_cumdur_days",
        # abx narrow
        "abx_narrow_last_t_dmonths",
        "abx_narrow_last_dur_days",
        "abx_narrow_cumcount",
        "abx_narrow_cumdur_days",
        # abx total count overall (used for plotting)
        "abx_max_count_ever",
    ]
    md_df = md_df[cols_to_keep]
    return md_df, excl_t1d_samples


def filter_sequences(path_to_data, excl_t1d_samples):
    path_otu = [
        os.path.join(path_to_data, f"otu_table_vatanen19_{subcohort}.qza")
        for subcohort in ["abx", "karelia", "t1d"]
    ]
    freq_df = get_merged_ft_frequency(path_otu)

    # remove abx samples from t1d subcohort
    freq_df = freq_df.loc[~freq_df.index.isin(excl_t1d_samples)]

    # remove shallow samples
    freq_df = remove_shallow_samples(freq_df)

    # identify and remove mitochondrial reads
    seq_paths = [x.replace("table", "seq") for x in path_otu]
    seq_table = get_merged_sequences(seq_paths)
    file_tax_classifier = "silva-138.1-ssu-nr99-515f-806r-classifier.qza"
    df_taxonomy = perform_taxonomic_classification(
        path_to_data, file_tax_classifier, seq_table
    )
    mitochondria_reads = df_taxonomy[
        df_taxonomy["Taxon"].str.contains("Mitochondria")
    ].index.tolist()

    print(f"Shape before mitochondria removal: {freq_df.shape}")
    freq_df = freq_df.drop(columns=mitochondria_reads)
    print(f"Shape after: {freq_df.shape}")

    print("If keeping all OTUs:\n")
    print(f"nb samples: {freq_df.shape[0]}")
    print(f"nb features: {freq_df.shape[1]}")
    return freq_df, df_taxonomy


def find_enteropathogens(df_taxonomy):
    # find enteropathogenic family reads
    f_enterobacteriaceae = df_taxonomy[
        df_taxonomy["Taxon"].str.contains("f__Enterobacteriaceae")
    ].index.tolist()

    # find genus level enterophathogenic reads
    g_groups_enteropathogenic = [
        "g__Clostridium",
        "g__Salmonella",
        "g__Clostridioides",
        "g__Escherichia",
        "g__Shigella",
        "g__Streptococcus",
    ]

    g_enteropathogens = []
    for genus in g_groups_enteropathogenic:
        g_enteropathogens += df_taxonomy[
            df_taxonomy["Taxon"].str.contains(genus)
        ].index.tolist()
    g_enteropathogens = list(set(g_enteropathogens))

    return {"family": f_enterobacteriaceae, "genus": g_enteropathogens}


def calculate_bootstrapped_alpha_metrics(md_df, freq_df, path_to_data, n=500):
    # rarefaction depth was previously evaluated with this command:
    # diversity.actions.alpha_rarefaction(freq_art, max_depth=10000, steps=20,
    # metadata=transform_to_q2metadata(md_df))
    rarefaction_depth = 1000
    alpha_all = pd.DataFrame(index=freq_df.index)
    ls_div_metrics = ["shannon", "observed_features", "faith_pd"]
    for div_metric in ls_div_metrics:
        alpha_mean_div = bootstrapped_alpha_div(
            div_metric, freq_df, rarefaction_depth, path_to_data, n=n
        )
        alpha_all = alpha_all.merge(
            alpha_mean_div.to_frame(), left_index=True, right_index=True
        )
    # add to metadata
    md_df = md_df.merge(alpha_all, left_index=True, right_index=True)
    return md_df


def get_relative_abund_n_sampling_depth(md_df, freq_df):
    # get sampling depth into md_df
    freq_df["sampling_depth"] = freq_df.sum(axis=1)
    md_df = pd.merge(
        md_df, freq_df[["sampling_depth"]], left_index=True, right_index=True
    )
    freq_df.drop(columns=["sampling_depth"], inplace=True)

    # get relative abundances
    freq_df = get_relative_abundance(freq_df, None)

    return md_df, freq_df


def create_feature_table(tag, path_to_data="../data/raw_data"):
    # filename of resulting datasets of this notebook
    path_to_ft = "../data/original_data"
    if not os.path.exists(path_to_ft):
        os.makedirs(path_to_ft)

    # get and process metadata
    md_df, excl_t1d_samples = filter_metadata(
        path_md=os.path.join(path_to_data, f"metadata_proc_v{tag}.tsv")
    )

    # get and process sequences
    freq_df, df_taxonomy = filter_sequences(path_to_data, excl_t1d_samples)

    # filter samples from freq_df in md_df
    print(f"Shape before filtering: {md_df.shape}")
    md_df = filter_md_by_ft(md_df, freq_df)
    print(f"Shape after filtering: {md_df.shape}")
    print(f"Unique host count: {md_df.host_id.nunique()}")

    # find entropathogenic reads
    enteropathogens = find_enteropathogens(df_taxonomy)

    # calculate and add bootstrapped diversity metrics
    md_df = calculate_bootstrapped_alpha_metrics(md_df, freq_df, path_to_data)

    # get sampling depth
    md_df, freq_df = get_relative_abund_n_sampling_depth(md_df, freq_df)

    # merge and save final datasets
    for name, ls_entero in enteropathogens.items():
        new_col_name = f"rel_abd_enteropathogens_{name}"
        ft_new = pd.DataFrame(index=freq_df.index)
        ft_new[new_col_name] = freq_df[ls_entero].sum(axis=1)

        # merge with md_df
        print(ft_new.shape)
        print(md_df.shape)
        ft_merged = pd.merge(ft_new, md_df, left_index=True, right_index=True)
        print(ft_merged.shape)

        # resort columns
        cols = ft_merged.columns.tolist()
        cols.remove(new_col_name)
        ft_merged = ft_merged[cols + [new_col_name]].copy()

        # save
        output_filename = os.path.join(
            path_to_ft, f"ft_vat19_anomaly_v{tag}_entero_{name}.tsv"
        )
        ft_merged.to_csv(output_filename, sep="\t")
        print(f"Saved {name} data to {output_filename}")


if __name__ == "__main__":
    # get user inputs
    parser = argparse.ArgumentParser(description="Create feature table.")
    parser.add_argument(
        "--tag",
        type=str,
        required=True,
        help="Tag of metadata and sequences to be used.",
    )
    args = parser.parse_args()

    create_feature_table(tag=args.tag)
