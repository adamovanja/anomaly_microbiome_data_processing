# General utility functions

import os
import pathlib
import subprocess
import tempfile

import matplotlib.pyplot as plt
import pandas as pd
import qiime2 as q2


def extract_color(style_name: str, color_index: int) -> str:
    """
    Extracts a specific color from a given style in matplotlib.

    Parameters:
    style_name (str): Name of the style in matplotlib
    color_index (int): Index of the color to be extracted

    Returns:
    str: The extracted color
    """
    plt.style.use(style_name)
    color_cycle = plt.style.library[style_name]["axes.prop_cycle"]

    # Extract the color at the specified index
    color = color_cycle.by_key()["color"][color_index]

    return color


def filter_md_by_ft(md: pd.DataFrame, ft: pd.DataFrame) -> pd.DataFrame:
    """
    Filters metadata by samples in feature table

    Parameters:
    md (pd.DataFrame): The metadata as a pandas DataFrame with a
    "host_id" column.
    ft (pd.DataFrame): The feature table.

    Returns: pd.DataFrame: Filtered metadata
    """
    ft_sample_ls = ft.index
    md_filt = md[md.index.isin(ft_sample_ls)]

    # assert that all samples in md are in ft_sample_ls
    assert len([x for x in ft_sample_ls if x not in md.index.tolist()]) == 0

    return md_filt


def transform_to_q2metadata(df_md):
    """Function that replaces column types not supported in Q2 with strings and
    returns respective Q2 metadata"""
    df_md = df_md.replace({True: "True", False: "True"}).copy()

    # convert boolean to string
    for col in df_md.columns:
        values_col = df_md[col].unique().tolist()
        if "True" in values_col or "False" in values_col:
            df_md[col] = df_md[col].astype("str")

    # convert date to string
    if "collection_date" in df_md.columns:
        df_md["collection_date"] = df_md["collection_date"].astype("str")
    if "sample_id" in df_md.columns:
        df_md["sample_id"] = df_md["sample_id"].astype("str")

    # transform to metadata artifact
    q2_md = q2.Metadata(df_md)

    return q2_md


def load_classifier(
    path_to_tax_classifier: str, file_tax_classifier: str
) -> q2.Artifact:
    path2save = os.path.join(path_to_tax_classifier, file_tax_classifier)
    if not os.path.isfile(path2save):
        command = f"srcd/get_silva_data.sh {path_to_tax_classifier} 6"
        subprocess.run(command, shell=True)

    return q2.Artifact.load(path2save)


def extract_file(viz, file):
    """
    This function reads specified in `file` from Q2 visualization artifact
    into a pandas dataframe. Only `tsv`, `csv` and `html` files are supported for now.
    """
    # code inspired by this Q2 forum post:
    # https://forum.qiime2.org/t/how-to-save-the-csv-create-a-table-from-the-barplot-visualisation-using-qiime2-api/17801/3
    with tempfile.TemporaryDirectory() as tmp:
        # export `data` directory from visualization into tmp
        viz.export_data(tmp)
        tmp_pathlib = pathlib.Path(tmp)
        extr = None
        for f in tmp_pathlib.iterdir():
            # print(f)
            if f.name == file:
                if f.name.endswith(".tsv"):
                    extr = pd.read_csv(f, sep="\t", index_col=0)
                elif f.name.endswith(".csv"):
                    extr = pd.read_csv(f, index_col=0)
                elif f.name.endswith(".html"):
                    extr = pd.read_html(f)
        if extr is None:
            raise ValueError(
                f"Requested file {file} does not exist or its format is not"
                " supported yet."
            )
        else:
            return extr
