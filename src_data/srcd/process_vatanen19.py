"""Module to process metadata of study by vatanen19"""

import numpy as np
import pandas as pd
import srcd.meta_processor as mproc
from srcd.process_vatanen19_abx import check_increasing_abx_cum

sra_cols_remove = [
    "Organism",
    "Library Source",
    "Library Selection",
    # # needed for duplicate check
    # "Bases",
    # "Spots",
    # "Bytes",
    "Avg Spot Len",
    "Biosamplemodel 1 [sample]",
    "Biosamplemodel 2 [sample]",
    "Collection Date [run]",
    "Collection Date [sample]",
    "Env Biome [run]",
    "Env Biome [sample]",
    "Env Feature [run]",
    "Env Feature [sample]",
    "Env Material [run]",
    "Env Material [sample]",
    "Geo Loc Name [sample]",
    "Host Age [sample]",
    "Host Age [run]",
    "Host Subject Id [sample]",
    "Host [run]",
    "Host [sample]",
    "Lat Lon [run]",
    "Lat Lon [sample]",
    # 'Library Name',
    "Organism [run]",
    "Tax ID",
]


def process_samples(path2supp):
    """Process tab samples from supp material table"""
    stab_samples = pd.read_excel(path2supp, sheet_name="Samples")

    # indication that samples are 16s not WGS - to be used as filter
    # for SRA metadata
    ls_16s_samples = stab_samples["gid_16s"].unique().tolist()

    # choose only 16s sequenced samples
    stab_samples = stab_samples[~stab_samples["gid_16s"].isna()]

    # select only columns that are not covered by SRA metadata and of relevance
    stab_samples = stab_samples[
        ["subjectID", "sampleID", "cohort", "age_at_collection"]
    ]
    stab_samples.rename(
        columns={"cohort": "study_subcohort", "age_at_collection": "age_days"},
        inplace=True,
    )

    return stab_samples, ls_16s_samples


def process_pregnancy(path2supp):
    """Process tab pregnancy from supp material table"""
    stab_birth = pd.read_excel(path2supp, sheet_name="Pregnancy, birth")

    # select columns that were covered by other studies
    stab_birth = stab_birth[["subjectID", "csection", "gender", "HLA_risk_class"]]
    # note: antibiotics information is given from additional material in scripts
    # proc_vatanen_abx

    stab_birth.rename(
        columns={
            "csection": "delivery_mode",
            "gender": "sex",
            "HLA_risk_class": "diag_HLA_risk_class",
        },
        inplace=True,
    )

    stab_birth["delivery_mode"] = stab_birth["delivery_mode"].astype(bool)
    stab_birth["delivery_mode"].replace(
        {True: "Cesarean", False: "Vaginal"}, inplace=True
    )
    return stab_birth


def add_milk_n_diet(path2supp, stab):
    """Process tab milk and diet from supp material table"""
    # ! Sheet milk
    stab_milk = pd.read_excel(path2supp, sheet_name="Mlik")

    # with the entries in the sheet "Early diet" it was assumed that the "milk"
    # sheet entries do not provide any information on weaning

    # this was retrieved as the distinction between
    # "bf_length" and "bf_length_exclusive":
    # "bf_length_exclusive" - displays start of weaning (zeros don't make sense here)
    # "bf_length" - displays end of milk diet

    # for now interpret as (confirmed by authors via email except for A2.2):
    # (A1)age <= "bf_length_exclusive" as excl. breast-feeding diet (marked as "bd")
    # (A2)"bf_length" > age > "bf_length_exclusive" as mixed diet
    # (A2.2) 12m >= age > "bf_length_exclusive" & no "bf_length" -> formula
    # (A3)if weaning >> bf_length -> "fd" diet
    # merge and create fields diet_milk & diet_weaning according to these assumptions
    stab = stab.merge(stab_milk, on="subjectID", how="left")

    stab["diet_milk"] = np.NaN
    stab["diet_milk"] = stab["diet_milk"].astype(str)
    # fill with zero - to allow for A2 to be filled out properly
    stab.fillna({"bf_length_exclusive": 0}, inplace=True)

    # (A1)age <= "bf_length_exclusive" as excl. breast-feeding diet (marked as "bd")
    stab.loc[stab["age_days"] <= stab["bf_length_exclusive"], "diet_milk"] = "bd"

    # (A2) "bf_length" > age > "bf_length_exclusive" as mixed diet
    stab.loc[
        np.logical_and(
            (stab["bf_length_exclusive"] <= stab["age_days"]),
            (stab["age_days"] <= stab["bf_length"]),
        ),
        "diet_milk",
    ] = "mixed"
    # (A2.2) 12m >= age > "bf_length_exclusive" & no "bf_length" -> formula
    # we assume that formula never stops here since we do not have a stop date
    stab.loc[
        np.logical_and(
            # age larger than excl. breastfeeding
            np.logical_and(
                (stab["bf_length_exclusive"] <= stab["age_days"]),
                # no info on bf_length though
                (stab["bf_length"].isna()),
            ),
            # infant is younger than 12 months
            stab["age_days"] <= 365.0,
        ),
        "diet_milk",
    ] = "fd"

    # ! add Sheet diet
    # weaning information is retrieved here from sheet "Early diet" - not from milk
    # (A3) if weaning >> bf_length -> set diet_milk to "fd"
    stab_diet = pd.read_excel(path2supp, sheet_name="Early diet")
    # remove entries from table that are <1 since they are not rounded to 0.5 as
    # others and they look like false entries (infants of <1 month can't have
    # been fed rye, tomatoes, ice cream etc.)
    stab_diet = stab_diet[stab_diet["start_month"] >= 1.0].copy()

    # take earliest start of diet in infants
    stab_diet = stab_diet[["subjectID", "start_month"]].groupby("subjectID").min()
    stab_diet.rename(columns={"start_month": "start_month_weaning"}, inplace=True)
    stab_diet.reset_index(inplace=True)

    stab = stab.merge(stab_diet, on="subjectID", how="left")

    # estimate age in months for diet estimation
    # (30.437 is avg. number of days per month)
    stab["age_months"] = stab["age_days"] / 30.437

    # start of weaning
    stab["diet_weaning"] = np.NaN

    # round start_month_weaning to 0.5
    stab["start_month_weaning"] = round(stab["start_month_weaning"] * 2 / 2)

    # rounding age also to nearest 0.5 for comparison
    stab.loc[
        ((round(stab["age_months"]) * 2 / 2) < stab["start_month_weaning"]),
        "diet_weaning",
    ] = False

    stab.loc[
        ((round(stab["age_months"]) * 2 / 2) >= stab["start_month_weaning"]),
        "diet_weaning",
    ] = True

    # A3: if weaning False AND breastfeeding diet was recorded but
    # diet_milk at time of sampling is NaN -> set diet_milk = 'fd'
    # infants must have eaten something
    fd_condition = np.logical_and(
        stab["diet_weaning"] == False,
        stab["diet_milk"].isna(),
    )
    stab.loc[fd_condition, "diet_milk"] = "fd"

    # no info on formula, so no "no milk" definition possible

    # remove columns we do not need anymore
    stab.drop(
        columns=[
            "bf_length",
            "bf_length_exclusive",
            "milk_first_three_days",
            "start_month_weaning",
            "age_months",
        ],
        inplace=True,
    )
    return stab


def add_diabetes_health(path2supp, stab):
    """Process tab diabetes from supp material table"""
    stab_diab = pd.read_excel(path2supp, sheet_name="Diabetes")

    stab_diab.rename(
        columns={
            "age_second_pos": "diag_seroconv_age_days",
            "age_at_t1d_onset": "diag_t1d_age_days",
        },
        inplace=True,
    )

    # include individual antibodies if available in other studies
    stab_diab = stab_diab[
        ["subjectID", "diag_seroconv_age_days", "diag_t1d_age_days"]
    ].copy(deep=True)

    # do all patients have an assessment of t1d and seroconversion?
    ls_t1d_tested = stab_diab.subjectID.unique()
    ids_never_tested = [x for x in stab.subjectID.unique() if x not in ls_t1d_tested]
    # ids_t1d = stab_diab.loc[
    #     ~stab_diab["diag_t1d_age_days"].isna(), "subjectID"
    # ].unique()

    # ! Merge
    stab = stab.merge(stab_diab, on="subjectID", how="left")

    # ! create diag_t1d field - at time of sampling
    stab["diag_t1d_at_sampling"] = False

    # this subject was never assessed
    stab.loc[stab["subjectID"].isin(ids_never_tested), "diag_t1d_at_sampling"] = np.NaN
    stab.loc[
        stab["age_days"] >= stab["diag_t1d_age_days"], "diag_t1d_at_sampling"
    ] = True

    t1d_during_sampling = (
        stab[["subjectID", "diag_t1d_at_sampling"]]
        .drop_duplicates()
        .diag_t1d_at_sampling.sum()
    )
    print(f"Only {t1d_during_sampling} infant developed T1D during sampling")

    # ! create diag_seroconv field - at time of sampling
    stab["diag_seroconv_at_sampling"] = False

    # this subject was never assessed
    stab.loc[
        stab["subjectID"].isin(ids_never_tested), "diag_seroconv_at_sampling"
    ] = np.NaN
    stab.loc[
        stab["age_days"] >= stab["diag_seroconv_age_days"],
        "diag_seroconv_at_sampling",
    ] = True

    t1d_during_sampling = (
        stab[["subjectID", "diag_seroconv_at_sampling"]]
        .drop_duplicates()
        .diag_seroconv_at_sampling.sum()
    )
    print(f"{t1d_during_sampling} infant developed pos. seroconv during sampling")

    # ! define healthy
    # all infants carry "human leukocyte antigen (HLA) haplotypes conferring
    # increased risk to autoimmune disorders"
    stab["health_status_at_sampling"] = "genetic risk autoimmune disorders"
    stab.loc[
        ~stab["diag_t1d_age_days"].isna(), "health_status_at_sampling"
    ] = "t1d_diagnosis"  # developed T1D

    return stab


def process_supp_metadata(path2supp):
    """Process all data in supp material table"""
    # ! Sheet "Samples"
    stab_samples, ls_16s_samples = process_samples(path2supp)

    # ! Sheet "pregnancy, birth"
    stab_birth = process_pregnancy(path2supp)

    # ! merge: stab_samples + stab_birth
    stab = stab_samples.merge(stab_birth, how="left", on="subjectID")

    # ! add info in sheets "Milk" and "Diet"
    stab = add_milk_n_diet(path2supp, stab)

    # ! add info in sheet "Diabetes"
    stab = add_diabetes_health(path2supp, stab)

    # ! Read growth information
    # little age overlap during stool sampling
    # maybe include in the future
    # stab_growth = pd.read_excel(path2supp, sheet_name="Growth")
    # stab_growth["agedays_haz_last"].describe()
    # print(f"Only {round((stab['age_days']>1016).sum()/stab.shape[0],2)} "
    # f"of samples have haz, waz info during sampling")

    # ! Rename sample_id
    stab.rename(columns={"subjectID": "host_id", "sampleID": "sample_id"}, inplace=True)
    stab["sample_id"] = stab["sample_id"].astype(str)

    return stab, ls_16s_samples


def process_sra_metadata(study_map, meta, study_abbr, ls_16s_samples):
    """Process sra metadata"""
    proj_id = study_map[study_abbr]
    sra_md = mproc.select_study_data(proj_id, meta)

    # remove columns that we don't need + reformat
    sra_md.drop(
        columns=sra_cols_remove,
        inplace=True,
    )

    # rename and reformat some columns
    sra_md.rename(
        columns={
            "Geo Loc Name [run]": "geo_location_name",
            "Host Subject Id [run]": "host_id",
            "Name": "sample_id",
            "Center Name": "Insdc center name [sample]",
        },
        inplace=True,
    )
    # only select 16s samples not WGS:
    # excluding based on instrument
    # since hiseq was used for WGS and miseq for 16s
    # can be verified with ls_16s_samples
    print(sra_md.shape)
    sra_md = sra_md.loc[sra_md["Instrument"] != "Illumina HiSeq 2500"].copy(deep=True)

    print(sra_md.shape)
    print(sra_md[sra_md["Library Name"].isin(ls_16s_samples)].shape)

    sra_md.drop(columns=["Library Name"], inplace=True)

    return sra_md


def add_info_from_publication(df, study_abbr):
    """add information from publication text"""
    # ! subfragment info from official website: https://diabimmune.broadinstitute.org/
    # sequencing info in Gevers 14 > Caporaso et al., 2012
    # online material contains exact primer sequences:
    # https://static-content.springer.com/esm/art%3A10.1038%2Fismej.2012.8/MediaObjects/41396_2012_BFismej20128_MOESM81_ESM.txt
    df["exp_target_subfragment"] = "V4"
    df["exp_primer"] = "515F [GTGCCAGCMGCCGCGGTAA], 806R [GGACTACHVGGGTWTCTAAT]"
    df["exp_bead_beating"] = np.NaN  # not available
    df["exp_target_gene"] = "16s rRNA"

    # ! add study index and name
    df["study_name"] = study_abbr

    # ! add lat & lon info by capital of the three countries
    # (as no other info given)
    # helsinki
    bool_finland = df.geo_location_name == "Finland"
    df.loc[bool_finland, "geo_latitude"] = 60.1699
    df.loc[bool_finland, "geo_longitude"] = 24.9384

    # Moscow
    bool_finland = df.geo_location_name == "Russia"
    df.loc[bool_finland, "geo_latitude"] = 55.7558
    df.loc[bool_finland, "geo_longitude"] = 37.6173

    # Tallinn
    bool_finland = df.geo_location_name == "Estonia"
    df.loc[bool_finland, "geo_latitude"] = 59.4370
    df.loc[bool_finland, "geo_longitude"] = 24.7536

    return df


def check_n_remove_duplicates(md_both):
    """Identify and remove duplicate samples"""
    # get run IDs of duplicates (2 counts for same sample_id)
    ser_count = md_both["sample_id"].value_counts()
    double_ls = ser_count[ser_count > 1].index

    # remove duplicates with smaller number of bytes
    # (in 336 cases same for two duplicates - only 54 slightly different)
    check_cols = ["sample_id", "Bases", "Spots"]
    dupl_md = md_both.loc[
        md_both.sample_id.isin(double_ls), ["Run ID", "Bytes"] + check_cols
    ]
    # sort by sample and bytes in asc order
    dupl_md.sort_values(["sample_id", "Bytes"], ascending=True, inplace=True)
    # run IDs with smaller bytes to remove from md_both:
    ids_to_remove = dupl_md.loc[
        dupl_md[check_cols].duplicated(keep="last"), "Run ID"
    ].unique()
    md_both = md_both[~md_both["Run ID"].isin(ids_to_remove)]

    # drop columns that were kept only for duplicate checks
    ls_cols_delete = ["Bases", "Spots", "Bytes"]
    md_both.drop(columns=ls_cols_delete, inplace=True)

    return md_both


def merge_abx_to_md(md_supp, df_abx, type: str):
    # join the two dataframes fully
    if type in ["broad", "narrow"]:
        df_abx = df_abx.loc[
            df_abx["abx_spectrum"] == type,
            ["host_id", "abx_start_age_months", "abx_duration_days"],
        ].copy()
        ls_nan_cols = [
            "abx_start_age_months",
            "abx_duration_days",
            "diff_age_since_abx",
        ]
    elif type == "any":
        ls_id = ["host_id"]
        ls_cols = [
            "abx_start_age_months",
            "abx_duration_days",
            "abx_broad_cumcount",
            "abx_broad_cumdur_days",
            "abx_narrow_cumcount",
            "abx_narrow_cumdur_days",
        ]
        df_abx = df_abx[ls_id + ls_cols].copy()
        ls_nan_cols = ls_cols + ["diff_age_since_abx"]

    md_supp = md_supp.sort_values(["host_id", "age_days"]).copy()
    md_merged = pd.merge(md_supp, df_abx, how="left", on=["host_id"])
    # select the last abx entry per sample row
    md_merged["diff_age_since_abx"] = (
        md_merged["age_months_rounded05"] - md_merged["abx_start_age_months"]
    )

    # remove entries of rows that have future abx entries matched
    bool_future_abx = md_merged["diff_age_since_abx"] < 0
    md_merged.loc[bool_future_abx, ls_nan_cols] = np.nan
    md_merged = md_merged.sort_values(
        ["host_id", "age_days", "diff_age_since_abx", "abx_duration_days"],
        ascending=[True, True, True, False],
    )

    # pick first entry per host_id and age_months_rounded
    # if same age_months two abx of two spectra were taken - select the longer one
    # if two samples were taken in same month - we are selecting the one closer to the
    # abx entry
    md_merged = (
        md_merged.groupby(["host_id", "age_days"], as_index=False).first().copy()
    )
    # rename columns
    md_merged.rename(
        columns={"diff_age_since_abx": f"abx_{type}_last_t_dmonths"},
        inplace=True,
    )
    md_merged.rename(
        columns={"abx_duration_days": f"abx_{type}_last_dur_days"}, inplace=True
    )
    md_merged.drop(columns=["abx_start_age_months"], inplace=True)

    return md_merged


def merge_supp_w_abx(md_supp, md_abx):
    # prep existing vatanen md_supp
    md_supp = md_supp.sort_values(["host_id", "age_days"])
    hosts_w_abx_info = md_supp.loc[
        md_supp.study_subcohort.isin(["abx", "karelia"]), "host_id"
    ].unique()

    # round age_days in md_supp to months closest to 0.5 (same as
    # abx_start_age_months in md_abx)
    md_supp["age_months_rounded05"] = md_supp["age_days"] / mproc.DAYS_PER_MONTH
    md_supp["age_months_rounded05"] = (md_supp["age_months_rounded05"] * 2).round() / 2

    # merge any type spectrum information to md file: last_dur_days, last_t_dmonths
    md_merged_any = merge_abx_to_md(md_supp, md_abx, "any")
    assert md_merged_any.shape[0] == md_supp.shape[0]

    # merge broad spectrum information to md file: cumcount, cumdur_days,
    # last_dur_days, last_t_dmonths
    md_merged_any_b = merge_abx_to_md(md_merged_any, md_abx, "broad")
    assert md_merged_any.shape[0] == md_merged_any_b.shape[0]

    # merge narrow spectrum information to md file: cumcount, cumdur_days,
    # last_dur_days, last_t_dmonths
    md_merged_any_bn = merge_abx_to_md(md_merged_any_b, md_abx, "narrow")
    assert md_merged_any_b.shape[0] == md_merged_any_bn.shape[0]

    # adding cumcount, cumdur_days for any type spectrum
    for cum in ["cumcount", "cumdur_days"]:
        md_merged_any_bn[f"abx_any_{cum}"] = (
            md_merged_any_bn[f"abx_narrow_{cum}"] + md_merged_any_bn[f"abx_broad_{cum}"]
        )
    # impute certain np.NaN values with 0 - only for abx cohorts (karelia + abx!)
    cols_impute = [
        x
        for x in md_merged_any_bn.columns
        if x.endswith("_cumcount") or x.endswith("cumdur_days")
    ]
    host_imputable = md_merged_any_bn.host_id.isin(hosts_w_abx_info)
    md_merged_any_bn.loc[host_imputable, cols_impute] = md_merged_any_bn.loc[
        host_imputable, cols_impute
    ].fillna(0)

    # add max count of abx exposures ever per host_id (indep. of microbiome
    # sampling)
    abx_ever_df = md_abx[["host_id", "abx_max_count_ever"]].drop_duplicates()
    md_merged_any_bn = md_merged_any_bn.merge(abx_ever_df, on="host_id", how="left")
    md_merged_any_bn["abx_max_count_ever"] = md_merged_any_bn[
        "abx_max_count_ever"
    ].fillna(0)

    # Hotfix: reset missing diet_weaning values from None to np.NaN
    md_merged_any_bn["diet_weaning"] = md_merged_any_bn["diet_weaning"].fillna(np.NaN)

    return md_merged_any_bn


def validate_merged_abx_entries(md_supp: pd.DataFrame) -> None:
    """
    This function validates the md_supp DataFrame based on several assertions:
    - Host IDs with antibiotic info should only be from 'abx' and 'karelia'
    study subcohorts.
    - Cumulative count and cumulative duration in days are continuously
    increasing per host ID.
    - Values in cumulative duration in days columns change only when cumulative
    count changes as well.

    Args:
        md_supp (pd.DataFrame): The DataFrame to validate.

    Raises:
        AssertionError: If any of the assertions fail.
    """
    # assert host_ids with abx info should only be from abx and karelia
    # study_subcohorts
    agg_df = (
        md_supp[["host_id", "abx_any_last_t_dmonths", "study_subcohort"]]
        .groupby("host_id")
        .agg(
            {
                "abx_any_last_t_dmonths": "max",
                "study_subcohort": lambda x: x.unique()[0],
            }
        )
    )

    filtered_df = agg_df[agg_df["abx_any_last_t_dmonths"].notna()]
    assert filtered_df["study_subcohort"].isin(["abx", "karelia"]).all(), (
        "Some rows with abx information do not belong to the study_subcohort"
        " 'abx' or 'karelia'"
    )

    # assert that cumcount and cumdur_days are still continuously increasing per
    # host id
    check_increasing_abx_cum(
        md_supp,
        ["broad", "narrow", "any"],
        ["cumcount", "cumdur_days"],
        "age_months_rounded05",
    )

    # assert that values in cumdur_days columns change only when cumcount
    # changes as well
    for spec in ["broad", "narrow"]:
        for col in ["cumdur_days", "cumcount"]:
            md_supp[f"abx_{spec}_{col}_shift"] = md_supp.groupby("host_id")[
                f"abx_{spec}_{col}"
            ].shift(1, fill_value=0)
            md_supp[f"abx_{spec}_{col}_change"] = (
                md_supp[f"abx_{spec}_{col}_shift"] != md_supp[f"abx_{spec}_{col}"]
            )

        assert (
            md_supp[f"abx_{spec}_cumdur_days_change"]
            == md_supp[f"abx_{spec}_cumcount_change"]
        ).all(), f"for {spec} abx _cumdur_days changes differently than _cumcount"
        md_supp.drop(
            columns=[
                f"abx_{spec}_cumdur_days_shift",
                f"abx_{spec}_cumcount_shift",
                f"abx_{spec}_cumdur_days_change",
                f"abx_{spec}_cumcount_change",
            ],
            inplace=True,
        )

    # ! since some abx were measured without gut sample after, the summed duration
    # in md_abx does not need to equal the last value of cumdur_days (same holds for
    # count in md_abx and cumcount in df_supp). example host_ids where this does not
    # hold: E002681, E002825
