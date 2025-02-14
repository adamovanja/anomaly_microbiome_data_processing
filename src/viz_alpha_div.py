import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import AutoLocator, FixedLocator
from scipy.stats import mannwhitneyu, wilcoxon

warnings.simplefilter(action="ignore", category=FutureWarning)
plt.rcParams.update({"font.family": "DejaVu Sans"})
plt.style.use("tableau-colorblind10")


def read_and_prep_abx_exposure_data(path_to_abx):
    abx_df = pd.read_csv(path_to_abx, sep="\t")
    ls_abx_cols = ["host_id", "abx_start_age_months"]
    abx_df = abx_df[ls_abx_cols].sort_values(ls_abx_cols).drop_duplicates()

    # HOTFIX: one host namely has abx_start_date given with 7.6 instead of 7.5,
    # this was already wrongly written in the raw supp. data from original
    # authors
    abx_df.loc[abx_df["host_id"] == "E014403", "abx_start_age_months"] = 7.5

    return abx_df


def assign_columns_for_plots(md_df):
    md_df = md_df.assign(
        max_abx_w_microbiome=lambda df: df.groupby("host_id")[
            "abx_any_cumcount"
        ].transform("max"),
    )
    md_df["diet_milk_s"] = md_df["diet_milk"]
    md_df.loc[md_df["diet_milk"] == "mixed", "diet_milk_s"] = "bd"

    md_df.sort_values(
        [
            "abx_max_count_ever",
            "max_abx_w_microbiome",
            "host_id",
            "age_days",
        ],
        ascending=[True, True, True, True],
        inplace=True,
    )
    return md_df


def boxplot_all_div_metrics_over_time(
    noabx,
    div_metrics,
    x_axis,
    path2save,
    title="Development for infants w/o abx exposure",
):
    fig, axs = plt.subplots(
        4, 1, figsize=(14, 18), height_ratios=[1, 1, 1, 0.5], sharex=True, dpi=400
    )

    for i, metric in enumerate(div_metrics):
        sns.boxplot(x=x_axis, y=metric, data=noabx, ax=axs[i], color="lightblue")

        axs[i].set_xticklabels(axs[i].get_xticklabels(), rotation=90)
        metric_name = metric.replace("div_alpha_", "")
        if i == 0:
            axs[i].set_title(title)
        axs[i].set_ylabel(metric_name)
        if i < len(div_metrics) - 1:
            axs[i].set_xlabel("")
        else:
            # if it's the last plot
            axs[i].set_xlabel(x_axis)

    sns.countplot(x=x_axis, data=noabx, ax=axs[3], color="grey")
    axs[3].xaxis.set_major_locator(AutoLocator())
    # axs[3].xaxis.set_major_locator(FixedLocator(axs[3].get_xticks()))
    axs[3].set_xticklabels(axs[3].get_xticklabels(), rotation=90)
    axs[3].set_ylabel("Number of samples")
    fig.subplots_adjust()
    plt.tight_layout()
    filename = os.path.join(path2save, "alpha_noabx_boxplot_all.png")
    plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.show()


def lineplot_all_div_metrics_over_time(
    noabx, metric, x_axis, group_by_values, path2save
):
    # todo: fix to work independently of noabx properties
    # TODO: make bottom boxplot in matplotlib as above and sharex properly
    # TODO: remove x_label in top plots
    nb_covariates = len(group_by_values)
    fig, axs = plt.subplots(
        nb_covariates * 2,
        1,
        figsize=(12, 20),
        # figsize=(nb_covariates * 2, nb_covariates * 3.5),
        height_ratios=nb_covariates * [1, 0.3],
        dpi=400,
    )
    # min_age = noabx[x_axis].min()
    max_age = int(noabx[x_axis].max()) + 1
    index_df = pd.DataFrame({x_axis: np.arange(0, max_age, 1.0)})
    for i, group_by in enumerate(group_by_values):
        mean_df = noabx.groupby([x_axis, group_by])[metric].mean().reset_index()
        std_df = noabx.groupby([x_axis, group_by])[metric].std().reset_index()
        count_df = noabx.groupby([x_axis, group_by])[metric].count().reset_index()

        merged_df = pd.merge(
            mean_df, std_df, on=[x_axis, group_by], suffixes=("_mean", "_std")
        )
        merged_df = pd.merge(merged_df, count_df, on=[x_axis, group_by])
        merged_df.columns = list(merged_df.columns[:-1]) + ["count"]

        palette = sns.color_palette("husl", len(merged_df[group_by].unique()))
        line_plot = sns.lineplot(
            x=x_axis,
            y=f"{metric}_mean",
            hue=group_by,
            data=merged_df,
            ax=axs[2 * i],
            palette=palette,
        )
        axs[2 * i].legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        # ensure correct order of colors
        line_colors = [line.get_color() for line in line_plot.lines]
        color_dict = dict(zip(merged_df[group_by].unique(), line_colors))
        for name, group in merged_df.groupby(group_by):
            axs[2 * i].fill_between(
                group[x_axis],
                group[f"{metric}_mean"] - group[f"{metric}_std"],
                group[f"{metric}_mean"] + group[f"{metric}_std"],
                color=color_dict[name],
                alpha=0.3,
            )

        axs[2 * i].xaxis.set_major_locator(FixedLocator(np.arange(0, max_age, 1.0)))
        axs[2 * i].set_xticklabels(
            np.arange(0, max_age, 1.0), rotation=0, ha="center", fontsize=8
        )
        axs[2 * i].set_xticklabels(axs[2 * i].get_xticks().astype(int))

        metric_name = metric.replace("div_alpha_", "")
        group_by_name = group_by.replace("_", " ")
        axs[2 * i].set_title(
            f"Alpha diversity grouped by {group_by_name}", loc="center", fontsize=14
        )
        if metric_name == "faith_pd":
            metric_name = "Faith PD"
        axs[2 * i].set_ylabel(f"{metric_name}")
        axs[2 * i].set_xlim(-0.5, max_age - 0.5)
        axs[2 * i].set_xlabel("")

        count_df = pd.merge(index_df, count_df, on=x_axis, how="left")
        sns.barplot(
            x=x_axis,
            y=metric,
            data=count_df,
            ax=axs[2 * i + 1],
            hue=group_by,
            palette=palette,
        )
        axs[2 * i + 1].xaxis.set_major_locator(
            FixedLocator(axs[2 * i + 1].get_xticks())
        )
        axs[2 * i + 1].set_xticklabels(
            axs[2 * i + 1].get_xticklabels(), rotation=0, ha="center", fontsize=8
        )
        axs[2 * i + 1].set_xticklabels(axs[2 * i + 1].get_xticks().astype(int))
        axs[2 * i + 1].set_ylabel("# samples")
        axs[2 * i + 1].set_xlabel("Age [months]")
        axs[2 * i + 1].get_legend().remove()

    fig.subplots_adjust(hspace=-0.4, right=0.75, top=0.9)
    plt.tight_layout()
    filename = os.path.join(path2save, "alpha_all_lineplot_splitcov.pdf")
    print(filename)
    plt.savefig(filename, dpi=400, bbox_inches="tight", format="pdf")
    plt.show()


def calculate_mean_diversity_per_cov(noabx, metric, cov_groups):
    noabx_mean_metric = noabx.groupby(cov_groups)[metric].mean().reset_index()
    return noabx_mean_metric.rename(columns={metric: f"mean_{metric}"})


def perform_significance_tests(df, t0, t1_values, metric_to_evaluate="diff_metric"):
    results = []

    # Filter the DataFrame for t0
    df_t0 = df.loc[(df["diff_age_nth_abx"] == t0), ["host_id", metric_to_evaluate]]
    df_t0.rename(columns={metric_to_evaluate: "t0"}, inplace=True)

    for t1 in t1_values:
        if t1 == t0:
            (
                t1,
                h_stat_unpair,
                p_val_unpair,
                n_unpaired_samples,
                h_stat_pair,
                p_val_pair,
                n_paired_samples,
            ) = (t1, None, None, None, None, None, None)
        else:
            # Filter the DataFrame for each t1
            df_t1 = df.loc[
                (df["diff_age_nth_abx"] == t1), ["host_id", metric_to_evaluate]
            ]
            df_t1.rename(columns={metric_to_evaluate: "t1"}, inplace=True)

            # perform the mann-whitney u-test (unpaired/independent)
            t0_unpaired = df_t0["t0"].dropna()
            t1_unpaired = df_t1["t1"].dropna()
            # note: this includes counting potential doubles in case there are
            # multiple samples at t1 that can be compared to t0
            n_unpaired_samples = t1_unpaired.shape[0]
            if (t0_unpaired.shape[0] > 0) and (t1_unpaired.shape[0] > 0):
                h_stat_unpair, p_val_unpair = mannwhitneyu(
                    t0_unpaired, t1_unpaired, alternative="greater", method="exact"
                )
            else:
                h_stat_unpair, p_val_unpair = np.nan, np.nan

            # Perform the wilcoxon test (paired)
            df_t0_t1 = pd.merge(df_t0, df_t1, on="host_id", how="inner")
            df_t0_t1.dropna(inplace=True)
            n_paired_samples = df_t0_t1.shape[0]
            if n_paired_samples > 0:
                h_stat_pair, p_val_pair = wilcoxon(
                    df_t0_t1["t0"],
                    df_t0_t1["t1"],
                    alternative="greater",
                    method="exact",
                )
            else:
                h_stat_pair, p_val_pair = np.nan, np.nan
        results.append(
            (
                t1,
                h_stat_unpair,
                p_val_unpair,
                n_unpaired_samples,
                h_stat_pair,
                p_val_pair,
                n_paired_samples,
            )
        )

    # Create a DataFrame from the results
    result_df = pd.DataFrame(
        results,
        columns=[
            "t1",
            "H-statistic unpaired",
            "P-value unpaired",
            "# samples unpaired",
            "H-statistic paired",
            "P-value paired",
            "# samples paired",
        ],
    )
    result_df.set_index("t1", inplace=True)

    return result_df


def test_mean_zero(df, t_values, metric_to_evaluate="diff_metric"):
    """tests if distribution mean is significantly different from zero"""
    results = []
    for t in t_values:
        # get the samples at time t
        t_values = df.loc[(df["diff_age_nth_abx"] == t), metric_to_evaluate].dropna()
        if t_values.shape[0] > 0:
            h_stat_mzero, p_val_mzero = wilcoxon(t_values)
        else:
            h_stat_mzero, p_val_mzero = np.nan, np.nan
        results.append((t, h_stat_mzero, p_val_mzero))
    # create a dataframe
    result_df = pd.DataFrame(
        results,
        columns=["t1", "H-statistic zerom", "P-value zerom"],
    )
    result_df.set_index("t1", inplace=True)
    return result_df
