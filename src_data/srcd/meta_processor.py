"""Module to process metadata"""

import os
from typing import Tuple

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
