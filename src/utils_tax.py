# Utility functions needed for taxonomic classification for neural ODE anomaly
# detection project

import os
from collections import Counter

import pandas as pd
import qiime2 as q2
from qiime2.plugins import feature_classifier
from src.utils import load_classifier

TAX_RANKS = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def perform_taxonomic_classification(
    path_to_tax_classifier, file_tax_classifier, seq_table, cohort_name="vatanen19"
):
    """Perform taxonomic classification of samples in seq_table with a pretrained
    classifier

    Args:
        path_to_tax_classifier (_type_): _description_
        file_tax_classifier (_type_): _description_
        seq_table (_type_): _description_
        cohort_name (_type_): _description_

    Returns:
        _type_: _description_
    """
    tax_classifier = load_classifier(path_to_tax_classifier, file_tax_classifier)
    # perform classification
    path_to_taxonomy = os.path.join(
        path_to_tax_classifier, f"otu_taxonomy_{cohort_name}.qza"
    )

    if not os.path.isfile(path_to_taxonomy):
        print("Reclassifying features with SILVA classifier...")
        (art_taxonomy,) = feature_classifier.methods.classify_sklearn(
            seq_table, tax_classifier
        )
        art_taxonomy.save(path_to_taxonomy)
        df_taxonomy = art_taxonomy.view(pd.DataFrame)
    else:
        df_taxonomy = q2.Artifact.load(path_to_taxonomy).view(pd.DataFrame)
    return df_taxonomy


def extract_taxonomic_entity(taxonomy_df: pd.DataFrame, tax_rank: str) -> dict:
    """
    Extract the taxonomic entity from a taxonomy DataFrame.

    Parameters: taxonomy_df (pd.DataFrame): A pandas DataFrame containing
    taxonomy information with a 'Taxon' column.

    Returns: dict: A dictionary with OTU ids as keys and respective taxonomic
    entities as values.
    """
    taxonomy_refined = taxonomy_df.copy()
    tax_ranks_dic = {rank: index for index, rank in enumerate(TAX_RANKS)}

    idx = tax_ranks_dic[tax_rank]
    prefix = tax_rank[0] + "__"

    taxonomy_refined[tax_rank] = (
        taxonomy_refined["Taxon"].str.split(";").str[idx].str.strip()
    )
    # remove unnecessary [] characters
    taxonomy_refined[tax_rank] = taxonomy_refined[tax_rank].str.replace(
        r"\[", "", regex=True
    )
    taxonomy_refined[tax_rank] = taxonomy_refined[tax_rank].str.replace(
        r"\]", "", regex=True
    )
    # flag unknown values at this level
    taxonomy_refined[tax_rank] = taxonomy_refined[tax_rank].replace(
        prefix, f"{prefix}unknown"
    )
    taxonomy_refined[tax_rank] = taxonomy_refined[tax_rank].fillna(f"{prefix}unknown")

    return taxonomy_refined[tax_rank].to_dict()


def get_counts_tax_dict(tax_dict: dict) -> dict:
    """
    Returns a dictionary containing the count of unique orders present in the
    input dictionary.

    Args:
        tax_dict: A dictionary containing taxonomic information.

    Returns:
        A dictionary containing the count of unique orders present in the input
        dictionary.
    """
    # get counts of orders present in the dictionary:
    tax_counts = Counter(tax_dict.values())
    sorted_tax_counts = dict(
        sorted(tax_counts.items(), key=lambda x: x[1], reverse=True)
    )
    print(f"Count of unique orders: {len(sorted_tax_counts)}")
    return sorted_tax_counts


def aggregate_ft_by_taxonomy(
    feature_table_df: pd.DataFrame, tax_dict: dict
) -> pd.DataFrame:
    """
    Aggregate a feature table DataFrame by taxonomic order in order_dict.

    Parameters: feature_table_df (pd.DataFrame): A pandas DataFrame containing
    the feature table with OTU ids as index.
    order_dict (dict): A dictionary
    with OTU ids as keys and taxonomic orders as values.

    Returns: pd.DataFrame: A pandas DataFrame with the aggregated feature table
    grouped by taxonomic order.
    """
    # Create a DataFrame with OTU and corresponding taxonomic entity
    tax_df = pd.DataFrame.from_dict(tax_dict, orient="index", columns=["tax_order"])
    tax_df.index.name = "OTU"

    # Merge the feature table with the order DataFrame
    # todo: treat unassigned features: raising warning
    merged_df = feature_table_df.T.merge(tax_df, left_index=True, right_index=True)

    # Group the merged DataFrame by the tax_order column and sum the rel abundances
    aggregated_df = merged_df.groupby("tax_order").sum()
    aggregated_df.index.name = "id"

    return aggregated_df.T


def transform_ft_by_tax_entity(
    ft: pd.DataFrame,
    tax_entity: str,
    df_taxonomy: pd.DataFrame,
    refine_unknown: str = None,
):
    # get taxonomic entities: dic[OTU] = tax_entity at level of interest only
    tax_dict = extract_taxonomic_entity(df_taxonomy, tax_entity)

    if refine_unknown:
        unknown_dict = extract_taxonomic_entity(df_taxonomy, refine_unknown)
        # replace __unknowns with tax. entity of different rank
        for key, value in tax_dict.items():
            tax_unknown = f"{tax_entity[0]}__unknown"
            if value == tax_unknown:
                tax_dict[key] = f"{tax_unknown[:-3]}_{unknown_dict[key]}"

    ft_aggregated = aggregate_ft_by_taxonomy(ft, tax_dict)

    return ft_aggregated


def explore_tax_grouping(df_taxonomy, freq_df) -> pd.DataFrame:
    """
    Retrieves feature count and fraction of unknown features for each taxonomic
    entity in freq_df following the taxonomic ranking in df_taxonomy.
    """
    dimred_tax_df = pd.DataFrame(index=TAX_RANKS)
    dimred_tax_df.index.name = "tax. entity"

    for tax_entity in TAX_RANKS:
        ft_aggregated = transform_ft_by_tax_entity(freq_df, tax_entity, df_taxonomy)

        # get relative abundances
        ft_aggregated_rel = ft_aggregated.div(ft_aggregated.sum(axis=1), axis=0)

        # avg fraction of unknown features
        unknown_col = tax_entity[0] + "__unknown"
        if unknown_col in ft_aggregated_rel.columns:
            frac_unknown = round(ft_aggregated_rel[unknown_col].mean() * 100, 2)
        else:
            frac_unknown = 0
        dimred_tax_df.loc[tax_entity, "nb features"] = ft_aggregated_rel.shape[1]
        dimred_tax_df.loc[tax_entity, "avg fraction of unknown [%]"] = frac_unknown

    return dimred_tax_df


def replace_special_tax_characters(x: str) -> str:
    """
    Function that replaces special characters in x
    """
    ls_spec = [" ", "  ", "-", ".", "(", ")"]
    for spec in ls_spec:
        x = x.replace(spec, "_")

    return x
