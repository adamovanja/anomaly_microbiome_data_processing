"""Module to process abx metadata of all subcohorts covered in vatanen19"""

import itertools

import numpy as np
import pandas as pd


def get_abx_karelia(path2karelia: str) -> pd.DataFrame:
    """
    Read and process antibiotic (abx) metadata for KARELIA subcohort.

    Parameters:
    path2karelia (str): Path to the KARELIA metadata file.

    Returns:
    pd.DataFrame: A dataframe containing the abx metadata for the KARELIA subcohort.
    """
    karelia = pd.read_excel(path2karelia, sheet_name="Basics and medications")

    # only select participant and AAB columns
    age_cols = [x for x in karelia.columns if "_AAB" in x]
    karelia_s = karelia[["Participant"] + age_cols].copy()

    # Create a list of dataframes using a list comprehension
    cols_abx = [
        "abx_name",
        "abx_reason",
        "abx_start_age_months",
        "abx_duration_days",
    ]
    ls_abx_per_host = [
        # reshape each row to 25x4
        pd.DataFrame(
            karelia_s.iloc[i, 1:].values.reshape(25, 4),
            # assign col names
            columns=cols_abx,
        ).assign(
            # add host id in each row
            host_id=karelia_s.iloc[i, 0]
        )
        for i in range(len(karelia_s))
    ]

    # Concatenate all dataframes in the list into a single dataframe
    df_abx_all = pd.concat(ls_abx_per_host, ignore_index=True)

    # Remove rows where all abx columns are NaN
    df_abx_all = df_abx_all.dropna(subset=cols_abx, how="all")

    # replace "NK" entries with np.NaN
    df_abx_all = df_abx_all.replace("NK", np.NaN).copy()

    # transform types
    for col in ["abx_start_age_months", "abx_duration_days"]:
        df_abx_all[col] = df_abx_all[col].astype(float)

    for col in ["abx_name", "abx_reason"]:
        df_abx_all[col] = df_abx_all[col].astype(str)

    # reorganize columns to new order
    df_abx_all = df_abx_all[
        [
            "host_id",
            "abx_start_age_months",
            "abx_duration_days",
            "abx_name",
            "abx_reason",
        ]
    ].copy()
    return df_abx_all


def get_abx_cohabx(path2abx: str) -> pd.DataFrame:
    abx_subcohort = pd.read_excel(path2abx, sheet_name="Antibiotics")

    # rename columns
    abx_subcohort.rename(
        columns={
            "Subject": "host_id",
            "Age (months)": "abx_start_age_months",
            "Duration (days)": "abx_duration_days",
            "Antibiotic type": "abx_name",
            "Symptoms": "abx_reason",
        },
        inplace=True,
    )
    abx_subcohort.drop(columns=["Antibiotic code"], inplace=True)

    # replace "Not known" entries with np.NaN
    abx_subcohort = abx_subcohort.replace("Not known", np.NaN).copy()
    abx_subcohort = abx_subcohort.replace("Not_known", np.NaN).copy()

    # transform types
    for col in ["abx_start_age_months", "abx_duration_days"]:
        abx_subcohort[col] = abx_subcohort[col].astype(float)

    for col in ["abx_name", "abx_reason"]:
        abx_subcohort[col] = abx_subcohort[col].astype(str)
    return abx_subcohort


def map_abx_spectrum(abx_name: str) -> str:
    """
    Maps an antibiotic name to its corresponding antibiotic spectrum ('broad' or
    'narrow'). Groupings stems from LLM phind.com with the following query:
    "Please group the antibiotics in the below Python list according to their
    spectrum of activity (broad vs. narrow spectrum) using sources such as WHO,
    FDA, CDC and peer-reviewed publications from pubmed."

    Parameters: abx_name (str): The name of the antibiotic to map.

    Returns: str: The corresponding antibiotic spectrum ('broad' or 'narrow')
    for the input antibiotic name.
    """
    # ! grouping of abx could obviously be refined with expert knowledge

    broad_spectrum = [
        "Amoxicillin",
        "Amoxicillin and clavulanic acid",
        "Ampicillin",
        "Azithromycin",
        "Azitromycin",
        "Cefaclor",
        "Cefadroxil",
        "Cefalexin",
        "Cefazolin",
        "Cefexime",
        "Cefotaxim",
        "Cefotaxime",
        "Cefprozil",
        "Ceftriaxone",
        "Cefuroxime",
        "Clarithromycin",
        "Sulfamethoxazole and trimethoprim",
        "Trimethoprim",
        "Trimethoprim and sulfadiazine",
        "Trimetoprim",
        "Trimetoprime and sulfadiazine",
    ]

    narrow_spectrum = [
        "Benzylpenicillin",
        "Furazidin",
        "Gentamicin",
        "Isoniazid",
        "Midecamycin",
        "Netilmicin",
        "Nifuroxazide",
        "Nitrofurantoin",
        "Penicillin G",
        "Phenoxymethylpenicillin",
        "Systemic antibiotic NAS",
    ]
    abx_name = abx_name.strip()
    if abx_name in narrow_spectrum:
        return "narrow"
    elif abx_name in broad_spectrum:
        return "broad"
    else:
        return "unknown"


def add_abx_spectrum_cumcount(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds cumulative count columns for the 'broad' and 'narrow' antibiotic
    spectra to the input df.

    Parameters: df (pd.DataFrame): DataFrame containing antibiotic course data.

    Returns: pd.DataFrame: DataFrame with additional columns for the cumulative
    count of courses for the 'broad' and 'narrow' antibiotic spectra.
    """
    df = df.sort_values(["host_id", "abx_start_age_months"])

    df["spectrum_cumcount"] = df.groupby(["host_id", "abx_spectrum"]).cumcount() + 1
    for spec in ["broad", "narrow"]:
        col_name = f"abx_{spec}_cumcount"

        # create new column for this abx spectrum
        df_spec = df.loc[
            df["abx_spectrum"] == spec,
            ["host_id", "abx_start_age_months", "spectrum_cumcount"],
        ].copy()
        df_spec.rename(columns={"spectrum_cumcount": col_name}, inplace=True)
        df = pd.merge(df, df_spec, how="left", on=["host_id", "abx_start_age_months"])

        # fill NaN values with last available value by host_id group
        df[col_name] = df.groupby("host_id")[col_name].ffill()
        # replace remaining NaN values with 0
        df[col_name].fillna(0, inplace=True)
    df.drop(columns=["spectrum_cumcount"], inplace=True)
    return df


def add_abx_spectrum_cumduration(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds cumulative duration columns for the 'broad' and 'narrow' antibiotic
    spectra to the input DataFrame.

    Parameters: df (pd.DataFrame): DataFrame containing antibiotic course data.

    Returns: pd.DataFrame: DataFrame with additional columns for the cumulative
    duration of courses for the 'broad' and 'narrow' antibiotic spectra.
    """
    df = df.sort_values(["host_id", "abx_start_age_months"])

    for spec in ["broad", "narrow"]:
        col_name = f"abx_{spec}_cumdur_days"
        # apply groupby and cumsum on 'abx_duration_days' for this spec
        df_spec = df.loc[
            df["abx_spectrum"] == spec,
            ["host_id", "abx_start_age_months", "abx_duration_days"],
        ].copy()
        df_spec[col_name] = df_spec["abx_duration_days"].groupby(df["host_id"]).cumsum()
        df_spec.drop(columns=["abx_duration_days"], inplace=True)
        df = pd.merge(
            df, df_spec, how="left", on=["host_id", "abx_start_age_months"]
        ).copy()

        # fill NaN values with last available value by host_id group
        df[col_name] = df.groupby("host_id")[col_name].ffill()
        # replace remaining NaN values with 0
        df[col_name].fillna(0, inplace=True)
    return df


def join_unique(x):
    """Returns a sorted list of unique strings"""
    return sorted(list(set(x)))


def process_abx_metadata(path2karelia, path2abx):
    # get raw data
    karelia = get_abx_karelia(path2karelia)
    abx = get_abx_cohabx(path2abx)

    # merge both subcohorts
    df_abx_all = pd.merge(karelia, abx, how="outer")

    # replace np.NaN duration values with mean of abx_duration_days value over
    # both cohorts
    mean_abx_duration_d = round(df_abx_all["abx_duration_days"].mean(), 0)
    df_abx_all["abx_duration_days"] = df_abx_all["abx_duration_days"].fillna(
        mean_abx_duration_d
    )

    # sum duration if same abx prescribed multiple times in same age month
    # period (A) & join the abx_reasons if multiple are given (B) e.g. for A:
    # E033156 in age 14 and 15 months, for B: host "E024240" in age 19 and host
    # "T021613" in age 13
    df_abx_all = (
        df_abx_all.groupby(
            ["host_id", "abx_start_age_months", "abx_name"], as_index=False
        )
        .agg({"abx_duration_days": "sum", "abx_reason": join_unique})
        .copy()
    )

    # map spectrum of each abx prescribed
    df_abx_all["abx_spectrum"] = df_abx_all["abx_name"].apply(map_abx_spectrum)
    assert (
        "unknown"
        not in df_abx_all["abx_spectrum"].value_counts(dropna=False).index.tolist()
    )

    # drop not needed columns
    df_abx_all.drop(columns=["abx_name", "abx_reason"], inplace=True)

    # group same spectrum abx per host and age into one row (funfact host_id E003188
    # has special case: two times broad spectrum abx in same timeframe)
    df_abx_all = (
        df_abx_all.groupby(
            ["host_id", "abx_start_age_months", "abx_spectrum"], as_index=False
        )
        .agg({"abx_duration_days": "sum"})
        .copy()
    )

    # create cumulative count per group in abx_spectrum per host_id
    df_abx_all = add_abx_spectrum_cumcount(df_abx_all)
    df_abx_all = add_abx_spectrum_cumduration(df_abx_all)

    # add max count of abx exposures ever per host_id (indep. of microbiome
    # sampling)
    df_abx_all["abx_max_count_ever"] = df_abx_all.groupby("host_id")[
        "abx_start_age_months"
    ].transform("count")
    return df_abx_all


def check_increasing_abx_cum(
    df: pd.DataFrame, spec2check: list, cum2check: list, age_col: str
) -> None:
    """
    Check if the cumulative count or duration of antibiotics is continuously
    increasing by host_id and age.

    Args:
        df (pd.DataFrame): The input DataFrame containing the antibiotic
        metadata. spec2check (list): A list of antibiotic spectra to check (e.g.
        ["broad", "narrow"]). cum2check (list): A list of cumulative measures to
        check (e.g. ["count", "duration"]). age_col (str): The name of the
        column containing the age information.

    Raises:
        AssertionError: If any of the cumulative measures are not continuously
        increasing by host_id and age.

    Returns:
        None
    """
    for spec, agg in itertools.product(spec2check, cum2check):
        assert (
            df[["host_id", age_col, f"abx_{spec}_{agg}"]]
            .groupby(["host_id", age_col])
            .apply(lambda x: x[f"abx_{spec}_{agg}"].diff().lt(0).sum())
            .sum()
        ) == 0, f"abx_{spec}_{agg} is not continuously increasing by host_id."


def validate_abx_all(md_abx: pd.DataFrame) -> None:
    """
    This function validates the md_abx DataFrame based on several assertions.

    Args:
        md_abx (pd.DataFrame): The DataFrame to validate.

    Raises:
        AssertionError: If any of the assertions fail.
    """
    # assert that in one abx_start_age_months there shouldn't be >1 same type of
    # abx_spectrum entry for the same host_id
    assert (
        md_abx[
            [
                "host_id",
                "abx_start_age_months",
                "abx_spectrum",
                "abx_duration_days",
            ]
        ]
        .groupby(["host_id", "abx_start_age_months", "abx_spectrum"])
        .count()
        != 1
    ).values.sum() == 0, (
        "There are multiple entries for the same host_id, abx_start_age_months"
        " and abx_spectrum."
    )

    # assert that abx_duration_days isn't <=0
    assert (
        md_abx["abx_duration_days"] <= 0
    ).sum() == 0, "Some abx_duration values are negative or zero."

    # assert that abx_*_cumcount and abx_*_cumdur_days are each continuously
    # increasing by host_id
    check_increasing_abx_cum(
        md_abx,
        ["broad", "narrow"],
        ["cumcount", "cumdur_days"],
        "abx_start_age_months",
    )

    # assert that abx_*_cumdur_days is sum of abx_duration_days per abx_spectrum
    # and host_id
    for spec in ["broad", "narrow"]:
        grouped = (
            md_abx[md_abx["abx_spectrum"] == spec]
            .groupby(["host_id", "abx_spectrum"])
            .agg({"abx_duration_days": "sum", f"abx_{spec}_cumdur_days": "last"})
        )

        assert (
            grouped["abx_duration_days"] == grouped[f"abx_{spec}_cumdur_days"]
        ).all(), (
            f"abx_{spec}_cumdur_days is not sum of abx_duration_days per"
            " abx_spectrum and host_id."
        )


def add_abx_info_t1d(path2t1d, md_both):
    """
    Add abx_ever information to t1d subcohort in md_both
    """
    # get abx_ever information from supp. material
    t1d_abx = pd.read_excel(path2t1d, sheet_name="Data")
    col_abx = [x for x in t1d_abx.columns if "antibiotic drug" in x]
    # rename columns
    t1d_abx.rename(
        columns={"ID, E=Espoo, Finland, T=Tartu, Estonia": "host_id"},
        inplace=True,
    )

    t1d_abx = t1d_abx[["host_id"] + col_abx].set_index("host_id")
    no_abx_ever = t1d_abx.count(axis=1) == 0
    t1d_abx_ever = no_abx_ever.to_frame().reset_index()
    t1d_abx_ever.rename(columns={0: "abx_ever_t1d"}, inplace=True)

    # merge with df_both
    md_both.reset_index(inplace=True)
    md_both = md_both.merge(t1d_abx_ever, how="left", on="host_id")
    md_both.set_index("index", inplace=True)

    bool_t1d_subcohort = md_both["study_subcohort"] == "t1d"
    md_both.loc[bool_t1d_subcohort, "abx_ever"] = md_both.loc[
        bool_t1d_subcohort, "abx_ever_t1d"
    ]
    md_both.drop(columns=["abx_ever_t1d"], inplace=True)

    return md_both
