"""Module to process metadata"""

import os
from typing import Tuple

import numpy as np
import pandas as pd
import qiime2 as q2
from qiime2.plugins import fondue

# 30.437 is avg. number of days per month
DAYS_PER_MONTH = 30.437


def fetch_metadata(
    ids, email, n_jobs, path2md, path2failed
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch metadata of `ids` and store the artifacts in `path2` location.

    Args:
        ids (q2.sdk.result.Artifact): Ids to fetch metadata for.
        email (str): Email.
        n_jobs (int): Number of jobs for q2fondue get-metadata.
        path2md (str): Path to store metadata in.
        path2failed (str): Path to store failed ids in.

    Returns:
        pd.DataFrame: Dataframe of fetched metadata
        pd.DataFrame: Dataframe of failed run IDs.
    """
    if not os.path.isfile(path2md):
        # fetch metadata
        (
            meta_md,
            failed_runs,
        ) = fondue.actions.get_metadata(
            accession_ids=ids, email=email, n_jobs=n_jobs, log_level="INFO"
        )
        # save for future reuse
        meta_md.save(path2md)
        failed_runs.save(path2failed)
    else:
        meta_md = q2.Artifact.load(path2md)
        failed_runs = q2.Artifact.load(path2failed)
        print(f'Metadata was read from file "{path2md}"')

    return meta_md.view(pd.DataFrame), failed_runs.view(pd.DataFrame)


def select_study_data(study_selected, df_meta, dropna_cols=True):
    """Function that returns data with ID provided in
    `study_selected` from `df`.

    Args:
        study_selected (list): List of bioproject IDs to select.
        df_meta (pd.DataFrame): Dataframe with metadata of all studies.
        dropna_cols (bool, optional): Flag to indicate if columns
        with all NaN values should be removed. Defaults to True.

    Returns:
        pd.DataFrame: Dataframe with only metadata from study selected.
    """
    """Function that only returns data of study_selected"""
    print(f"BioProject ID: {study_selected}")

    # select only this study
    df_study = df_meta[df_meta["Bioproject ID"].isin(study_selected)].copy(deep=True)
    # drop columns with all NaNs
    if dropna_cols:
        df_study.dropna(axis=1, how="all", inplace=True)
    print(f"Shape of study dataset: {df_study.shape}")

    return df_study


def correct_weaning(df, host_id_col):
    """Once diet_weaning is true it should be always stay true per host"""
    cols_sort = ["study_name", host_id_col, "age_days"]
    df = df.sort_values(cols_sort).copy()

    # check values
    values_weaning = df["diet_weaning"].unique()
    invalid_values = [
        x for x in values_weaning if x not in [True, False] and not pd.isna(x)
    ]
    if len(invalid_values) > 0:
        raise Warning("diet_weaning column has incorrect values")
    # set weaning to True once it is True
    # helper column needed since all np.NaN values are counted as False when boolean
    df["diet_weaning_h"] = df["diet_weaning"].astype("boolean")
    df["diet_weaning_h"] = df["diet_weaning_h"].fillna(False)
    df["cumsum_weaning"] = (
        df[[host_id_col, "diet_weaning_h"]].groupby([host_id_col]).cumsum()
    )
    # replace weaning in these columns:
    df.loc[df["cumsum_weaning"] > 0, "diet_weaning"] = True
    # remove helper column
    df.drop(columns=["cumsum_weaning", "diet_weaning_h"], inplace=True)

    # df["diet_weaning"].value_counts(dropna=False)
    return df


def _replace_no_milk(series):
    """
    Set "no milk" entries to NaN in case milk diet is given in following months
    per host
    """
    replace = False
    for i in reversed(range(len(series))):
        if series.iat[i] in ["bd", "fd", "mixed"]:
            replace = True
        if series.iat[i] == "no milk" and replace:
            series.iat[i] = np.NaN
    return series


def correct_nomilk_diet(df, host_id_col):
    """
    Remove wrong "no milk" entries
    """
    df["diet_milk"] = df.groupby(host_id_col)["diet_milk"].transform(_replace_no_milk)
    return df


def _replace_after_finished(group, column, const_val):
    """
    once `column` value is equal to `const_val` it remains this value
    """
    if const_val in group[column].values:
        first_finished = group[column].eq(const_val).idxmax()
        group.loc[first_finished:, column] = const_val
    return group


def refine_weaning(df):
    """
    Transforms True/False values of diet_weaning and defines "finished" weaning state
    """
    # ensure correct sorting
    cols = ["study_name", "host_id", "age_days"]
    df = df.sort_values(cols)

    # correct weaning: once True always True
    df = correct_weaning(df, "host_id")

    # True/False to "yes"/"no"
    df = df.replace({"diet_weaning": {True: "yes", False: "no"}})

    # define "finished" weaning state
    cond_finished = np.logical_and(
        df["diet_milk"] == "no milk", df["diet_weaning"] == "yes"
    )
    df.loc[cond_finished, "diet_weaning"] = "finished"

    # once weaning was finished it remains finished
    df = df.groupby("host_id", group_keys=False).apply(
        lambda x: _replace_after_finished(x, "diet_weaning", "finished")
    )

    return df


def _strict_bd_distinction(group):
    """
    Once previous value is "mixed" or "fd", all following entries of "bd" are no
    longer allowed
    """
    prev_val = None
    replace_bd = False
    for idx, row in group.iterrows():
        if row["diet_milk"] in ["mixed", "fd"]:
            replace_bd = True
            prev_val = row["diet_milk"]
        if replace_bd and row["diet_milk"] == "bd":
            group.loc[idx, "diet_milk"] = prev_val
        elif not replace_bd:
            prev_val = row["diet_milk"]
    return group


def _fill_missing_values(df, colname):
    """Fills missing values in `colname` in case there are known previous and
    later values that match"""
    # Forward fill
    df_ffill = df.groupby("host_id")[colname].ffill()
    # Backward fill
    df_bfill = df.groupby("host_id")[colname].bfill()
    # Fill NaN values only if both forward and backward values are the same
    df[colname] = df[colname].fillna(df_ffill.where(df_ffill == df_bfill))
    return df


def refine_milk(df):
    """
    Ensures that once diet_milk is "no milk" it remains constant
    """
    # ensure correct sorting
    cols = ["study_name", "host_id", "age_days"]
    df = df.sort_values(cols)

    # once diet_milk is "no milk" it remains constant
    df = df.groupby("host_id", group_keys=False).apply(
        lambda x: _replace_after_finished(x, "diet_milk", "no milk")
    )

    # once host has "fd" or "mixed" - "bd" is no longer allowed
    df = df.groupby("host_id", group_keys=False).apply(_strict_bd_distinction)

    # fill missing entries between known values
    df = _fill_missing_values(df, "diet_milk")

    return df


def _check_nan_after_value(series, ls_values):
    """
    in series verify if after values in ls_values NaN occurs
    """
    flag = False
    for value in series:
        if value in ls_values:
            flag = True
        if pd.isna(value) and flag:
            return True
    return False


def _is_in_order(series):
    """Check if the ordered categorical data is sorted, ignoring NaNs"""
    return series.dropna().is_monotonic_increasing


def _verify_correct_order(df_check, colname, order):
    """
    Check that values in `colname` have order as defined in `order`
    """
    df_check[colname] = df_check[colname].astype(str)

    df_check[colname] = pd.Categorical(
        df_check[colname],
        categories=order,
        ordered=True,
    )
    is_ordered = df_check.groupby("host_id")[colname].apply(_is_in_order)
    if sum(is_ordered) != len(is_ordered):
        hosts_w_wrong_order = is_ordered[~is_ordered].index
        raise ValueError(
            "{colname} column is not ordered correctly for these hosts:"
            f" {hosts_w_wrong_order}."
        )


def check_weaning(df):
    """Check diet_weaning column for correct values."""
    df_check = df.copy()
    # no NaNs are allowed after "yes" or "finished"
    has_nan_after_value = df_check.groupby("host_id")["diet_weaning"].apply(
        _check_nan_after_value, ["yes", "finished"]
    )
    if sum(has_nan_after_value) > 0:
        raise ValueError('NaN after "yes" or "finished" in diet_weaning column.')

    # verify order of column
    _verify_correct_order(df_check, "diet_weaning", ["no", "yes", "finished"])


def check_milk(df):
    # "no milk" is only allowed if infant is weaned at the same time
    df_check = df.copy()
    assert (
        df_check.loc[
            np.logical_and(
                df_check["diet_milk"] == "no milk",
                df_check["diet_weaning"] == "no",
            )
        ].shape[0]
        == 0
    )

    # "no milk" is not allowed prior to 3 months
    assert (
        df_check[
            np.logical_and(
                df_check["diet_milk"] == "no milk",
                df_check["age_months_rounded05"] < 3,
            )
        ].shape[0]
        == 0
    )

    # no NaNs are allowed after "no milk"
    has_nan_after_value = df_check.groupby("host_id")["diet_milk"].apply(
        _check_nan_after_value, ["no milk"]
    )
    if sum(has_nan_after_value) > 0:
        raise ValueError(
            "NaN after 'no milk' in diet_milk column for these hosts:"
            f" {has_nan_after_value[has_nan_after_value].index}."
        )

    # verify order of column
    # create new column representing stages of mixed
    df_check["diet_stage"] = df_check["diet_milk"].replace(
        {"bd": "bd", "mixed": "mixed1", "fd": "fd", "no milk": "no milk"}
    )
    df_check.loc[
        (df_check["diet_stage"] == "mixed1")
        & (df_check["diet_stage"].shift(-1) == "fd"),
        "diet_stage",
    ] = "mixed2"
    order = ["bd", "mixed1", "fd", "mixed2", "no milk"]

    _verify_correct_order(df_check, "diet_milk", order)


def post_process_md(df_md):
    # ensure correct sorting
    cols = ["study_name", "host_id", "age_days"]
    df_md = df_md.sort_values(cols)
    # set run ID as index
    df_md.set_index("Run ID", inplace=True)
    df_md.index.name = "id"
    df_md.index = df_md.index.astype(str)

    # replace empty space with _
    df_md.columns = df_md.columns.str.replace(" ", "_")
    # rename
    df_md.rename(
        columns={"Insdc_center_name_[sample]": "insdc_center_name"}, inplace=True
    )
    # make all columns names lowercase
    df_md.columns = df_md.columns.str.lower()
    # transform values to lower case values of specified columns
    ls_lowercase = ["sex", "delivery_mode"]

    for col in ls_lowercase:
        df_md[col] = df_md[col].str.lower()

    # transform strings and digits to boolean values
    df_md = df_md.replace({"diag_diarrhea_7d_prior": {"Yes": True, "No": False}})

    ls_bool = [
        "diag_diarrhea_7d_prior",
        "exp_bead_beating",
        "abx_ever",
        "abx_7d_prior",
        "diag_t1d_at_sampling",
        "diag_seroconv_at_sampling",
    ]

    for col in ls_bool:
        df_md = df_md.replace({col: {1: True, 0: False}})

    # included non-rounded age in months
    df_md["age_months"] = df_md["age_days"] / DAYS_PER_MONTH

    # included 0.5- and fully 1.0 rounded age in months for all studies
    df_md["age_months_rounded05"] = (df_md["age_days"] / DAYS_PER_MONTH * 2).round() / 2
    df_md["age_months_rounded1"] = (df_md["age_days"] / DAYS_PER_MONTH).round()

    # exclude studies that did not sequence v4 region (oyedemi22 excluded)
    print(f"Shape before V4 selection: {df_md.shape}")
    df_md = df_md[df_md["exp_target_subfragment"] == "V4"].copy()
    print(f"Shape after V4 selection: {df_md.shape}")

    df_md["study_name"].value_counts()

    # adjust weaning
    df_md = correct_weaning(df_md, "host_id")
    df_md = correct_nomilk_diet(df_md, "host_id")

    df_md = refine_weaning(df_md)

    # sanity check of diet_weaning column
    check_weaning(df_md)

    # adjust diet_milk
    df_md = refine_milk(df_md)

    # sanity check of diet_milk column
    check_milk(df_md)

    # Create study_cohort columns (easier for separate denoising per subcohort
    # afterwards)
    df_md["study_cohort_name"] = df_md["study_name"] + "_" + df_md["study_subcohort"]

    # subramanian subcohort is only 1 - hence remove
    df_md.loc[
        df_md["study_cohort_name"] == "subramanian14_Healthy Twins & Triplets",
        "study_cohort_name",
    ] = np.nan
    df_md.loc[df_md["study_cohort_name"].isna(), "study_cohort_name"] = df_md[
        "study_name"
    ]

    df_md["study_cohort_name"].value_counts(dropna=False)
    return df_md


def save_file(df_md, path2data, tag):
    path2md = os.path.join(path2data, "metadata.qza")
    path2md_proc = path2md.replace(".qza", f"_proc_v{tag}.tsv")

    print(f"Saved processed metadata to: {path2md_proc}")
    df_md.to_csv(path2md_proc, sep="\t")

    # save unique runIDs of each projectIDs from this metadata file for sequence
    # fetching in d_
    unique_ids = df_md["bioproject_id"].unique()
    path2runids = os.path.join(path2data, "runids")
    if not os.path.exists(path2runids):
        os.makedirs(path2runids)

    for id in unique_ids:
        run_ids = df_md[df_md["bioproject_id"] == id].index.tolist()
        df_run_ids = pd.DataFrame({"id": run_ids}).squeeze()

        path2ids = os.path.join(path2runids, f"{id}.tsv")
        print(f"Saved unique project IDs to: {path2ids}")
        df_run_ids.to_csv(path2ids, sep="\t", index=False)


def transform_to_q2metadata(df):
    """Function that replaces column types not supported in Q2 with strings and
    returns respective Q2 metadata"""
    df = df.replace({True: "True", False: "True"}).copy()

    # convert boolean to string
    for col in df.columns:
        values_col = df[col].unique().tolist()
        if "True" in values_col or "False" in values_col:
            df[col] = df[col].astype("str")

    # convert date to string
    if "collection_date" in df.columns:
        df["collection_date"] = df["collection_date"].astype("str")

    # transform to metadata artifact
    q2_md = q2.Metadata(df)

    return q2_md
