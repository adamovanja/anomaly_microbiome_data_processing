"""Functions to fetch metadata and supplementary material"""

import os

import pandas as pd
import qiime2 as q2
import requests
from src import meta_processor as mproc


def _fetch_metadata(path2data, ids, email, n_jobs):
    path2md = os.path.join(path2data, "metadata.qza")
    path2failed = os.path.join(path2data, "metadata_failed_runs.qza")

    if not os.path.isfile(path2md):
        print("Fetching metadata...")
        meta, failed = mproc.fetch_metadata(ids, email, n_jobs, path2md, path2failed)
        assert failed.shape[0] == 0
        print(f"Metadata was fetched and saved to file {path2md}")
    else:
        meta = q2.Artifact.load(path2md)
        meta = meta.view(pd.DataFrame)
        print(f"Metadata was read from file {path2md}")

    # small postprocess
    meta = meta.reset_index()
    meta.rename(columns={"ID": "Run ID"}, inplace=True)
    col_to_remove = ["Platform", "Public"]
    return meta.drop(columns=col_to_remove)


def _make_dirs(path2dir):
    if not os.path.exists(path2dir):
        os.makedirs(path2dir)


def _fetch_n_store_excel_file(url, filedest):
    response = requests.get(url)
    with open(filedest, "wb") as f:
        f.write(response.content)


def _fetching_programmatically_not_allowed(url, filedest):
    raise ValueError(
        f"The metadata of this subcohort can't be fetched "
        f"programmatically due to the site's robots.txt policy. "
        f"Please visit the following URL to download the metadata "
        f'manually into "{filedest}": {url}'
    )


def _fetch_all_supp_material(path2data):
    paths_dict = {}
    dest_suppmat = os.path.join(path2data, "supp_material")
    _make_dirs(dest_suppmat)

    # get overall metadata
    filedest_all = os.path.join(dest_suppmat, "md.xlsx")
    if not os.path.isfile(filedest_all):
        url_all = "https://static-content.springer.com/esm/"
        "art%3A10.1038%2Fs41564-018-0321-5/MediaObjects/"
        "41564_2018_321_MOESM3_ESM.xlsx"
        _fetch_n_store_excel_file(url_all, filedest_all)
    paths_dict["all"] = filedest_all
    # get subcohort abx: Yassour16 metadata
    path_abx = os.path.join(dest_suppmat, "abx_md")
    _make_dirs(path_abx)
    filedest_abx = os.path.join(path_abx, "aad0917_Table S1.xls")
    if not os.path.isfile(filedest_abx):
        url_abx = "https://www.science.org/doi/suppl/10.1126/scitranslmed.aad0917/"
        "suppl_file/8-343ra81_table_s1.zip"
        _fetching_programmatically_not_allowed(url_abx, filedest_abx)
    paths_dict["abx"] = filedest_abx
    # get subcohort karelia
    path_karelia = os.path.join(dest_suppmat, "karelia_md")
    _make_dirs(path_karelia)
    filedest_karelia = os.path.join(path_karelia, "mmc2.xlsx")
    if not os.path.isfile(filedest_karelia):
        url_karelia = "https://www.cell.com/cms/10.1016/j.cell.2016.04.007/"
        "attachment/a61300b3-0fd7-43b1-acfc-4accd7e538de/mmc2.xlsx"
        _fetching_programmatically_not_allowed(url_karelia, filedest_karelia)
    paths_dict["karelia"] = filedest_karelia

    # get subcohort t1d
    path_t1d = os.path.join(dest_suppmat, "t1d_md")
    _make_dirs(path_t1d)
    filedest_t1d = os.path.join(path_t1d, "mmc2.xlsx")
    if not os.path.isfile(filedest_t1d):
        url_t1d = "https://www.cell.com/cms/10.1016/j.chom.2015.01.001/"
        "attachment/1f0883f8-1df7-447d-a47b-c1aa2bb2bbaf/mmc2.xlsx"
        _fetching_programmatically_not_allowed(url_t1d, filedest_t1d)
    paths_dict["t1d"] = filedest_t1d

    return paths_dict
