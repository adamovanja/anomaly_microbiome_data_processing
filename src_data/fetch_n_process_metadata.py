import os
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import qiime2 as q2
from srcd import meta_processor as mproc
from srcd import process_vatanen19 as proc_vatanen
from srcd import process_vatanen19_abx as proc_vatanen_abx
from srcd import utils_fetch as uf

warnings.filterwarnings("ignore", category=FutureWarning)


def fetch_n_process_metadata(
    bioproject_id: str, study_map: dict, email: str = "my@mail.com", n_jobs: int = 6
):
    # setup
    tag = datetime.today().strftime("%Y%m%d")
    path2data = "../data/raw_data/"
    if not os.path.exists(path2data):
        os.makedirs(path2data)

    # import bioproject id as NCBIAccessionIDs for q2-fondue
    sra_ids = pd.Series([bioproject_id], name="ID")
    ids = q2.Artifact.import_data("NCBIAccessionIDs", sra_ids)

    # fetch and store metadata w q2-fondue
    # (takes approx. 20min when needing to fetch for the first time)
    meta = uf._fetch_metadata(path2data, ids, email, n_jobs)

    # fetch supplementary material for this study
    # (due to journal's policies sadly requires manually downloading)
    supp_paths = uf._fetch_all_supp_material(path2data)

    # process all metadata of all subcohorts included
    md_supp, ls_16s = proc_vatanen.process_supp_metadata(supp_paths["all"])
    print(f"md_all shape: {md_supp.shape}")

    # get detailed abx intake information for subcohorts "abx" + "karelia"
    md_abx_all = proc_vatanen_abx.process_abx_metadata(
        supp_paths["karelia"], supp_paths["abx"]
    )
    proc_vatanen_abx.validate_abx_all(md_abx_all)

    # saving since this contains not only sample matched abx exposure events
    # ! since only karelia + abx have detailed abx info - only these are saved in ts
    cols2save = ["host_id", "abx_start_age_months", "abx_spectrum", "abx_duration_days"]
    md_abx_all[cols2save].to_csv(
        os.path.join(path2data, f"ts_vat19_abx_v{tag}.tsv"), sep="\t", index=False
    )

    # merge existing supp. material with additional abx information
    md_supp = proc_vatanen.merge_supp_w_abx(md_supp, md_abx_all)
    proc_vatanen.validate_merged_abx_entries(md_supp)

    # get SRA metadata
    study_abbr = next(iter(study_map))
    md_sra = proc_vatanen.process_sra_metadata(study_map, meta, study_abbr, ls_16s)
    md_sra = proc_vatanen.add_info_from_publication(md_sra, study_abbr)

    # merge supp material (incl. abx) and sra metadata
    print(f"Shape before merge: {md_sra.shape}")
    md_both = md_supp.merge(md_sra, on=["sample_id", "host_id"], how="left")
    print(f"Shape after merge: {md_both.shape}")

    # check and remove duplicate sequences
    print(f"Shape before duplicates removal: {md_both.shape}")
    md_both = proc_vatanen.check_n_remove_duplicates(md_both)
    print(f"Shape after duplicates removal: {md_both.shape}")

    # define abx_ever and abx_7d_prior to align with other studies
    # since some hosts have abx exposure only after last microbial gut sample ->
    # flagging these as False for our purpose

    # must set nan to -1 otherwise they get counted as 0
    md_both["abx_any_cumcount"] = md_both["abx_any_cumcount"].fillna(-1)
    md_both["max_abx_w_microbiome"] = md_both.groupby("host_id")[
        "abx_any_cumcount"
    ].transform("max")

    md_both["abx_ever"] = np.NaN
    md_both.loc[md_both["max_abx_w_microbiome"] == 0.0, "abx_ever"] = False
    md_both.loc[md_both["max_abx_w_microbiome"] > 0.0, "abx_ever"] = True
    # revert
    md_both.drop(columns=["max_abx_w_microbiome"], inplace=True)
    md_both["abx_any_cumcount"] = md_both["abx_any_cumcount"].replace(-1, np.NaN)

    assert md_both.loc[md_both["abx_ever"] == True, "host_id"].nunique() == 142
    assert md_both.loc[md_both["abx_ever"] == False, "host_id"].nunique() == 114

    # define abx_7d_prior -> approximated with start 0.5 months prior to sampling
    # since abx info only given on 0.5 monthly rate
    md_both["abx_7d_prior"] = np.NaN
    md_both.loc[md_both["abx_any_last_t_dmonths"] <= 0.5, "abx_7d_prior"] = True
    md_both.loc[md_both["abx_any_last_t_dmonths"] > 0.5, "abx_7d_prior"] = False

    # t1d samples should not have this entry - no detailed info available  for this
    bool_t1d_subcohort = md_both["study_subcohort"] == "t1d"
    t1d_count_samples = bool_t1d_subcohort.sum()
    assert (
        md_both.loc[bool_t1d_subcohort, "abx_7d_prior"].isna().sum()
        == t1d_count_samples
    )

    # add abx_ever information for infants from t1d subcohort
    md_both = proc_vatanen_abx.add_abx_info_t1d(supp_paths["t1d"], md_both)

    # post-process md_both
    md_both = mproc.post_process_md(md_both)

    # save to file
    mproc.save_file(md_both, path2data, tag)


if __name__ == "__main__":
    bioproject_id_vat19 = "PRJNA497734"
    study_map = {"vatanen19": [bioproject_id_vat19]}
    fetch_n_process_metadata(bioproject_id_vat19, study_map)
