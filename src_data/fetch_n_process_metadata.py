import os

import pandas as pd
import qiime2 as q2
from srcd import utils_fetch as uf


def fetch_n_process_metadata(
    bioproject_id: str, study_map: dict, email: str = "my@mail.com", n_jobs: int = 6
):
    # setup
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
    uf._fetch_all_supp_material(path2data)

    print(meta.columns)
    return meta


if __name__ == "__main__":
    bioproject_id_vat19 = "PRJNA497734"
    study_map = {"vatanen19": ["PRJNA497734"]}
    fetch_n_process_metadata(bioproject_id_vat19, study_map)
