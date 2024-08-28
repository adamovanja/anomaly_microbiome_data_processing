import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoLocator, FixedLocator
from scipy.stats import mannwhitneyu, wilcoxon

warnings.simplefilter(action="ignore", category=FutureWarning)
plt.rcParams.update({"font.family": "DejaVu Sans"})
plt.style.use("tableau-colorblind10")


def read_and_prep_abx_exposure_data(path_to_abx):
    abx_df = pd.read_csv(path_to_abx, sep="\t")
    ls_abx_cols = ["host_id", "abx_start_age_months"]
    abx_df = abx_df[ls_abx_cols].sort_values(ls_abx_cols).drop_duplicates()
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
        figsize=(nb_covariates * 2, nb_covariates * 3.5),
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
        axs[2 * i].set_xticklabels(np.arange(0, max_age, 1.0), rotation=90)

        metric_name = metric.replace("div_alpha_", "")
        axs[2 * i].set_title(f"Grouped by {group_by}", loc="left")
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
        axs[2 * i + 1].set_xticklabels(axs[2 * i + 1].get_xticklabels(), rotation=90)
        axs[2 * i + 1].set_ylabel("Number of samples")
        axs[2 * i + 1].get_legend().remove()

    fig.subplots_adjust(hspace=-0.4, right=0.75, top=0.9)
    plt.tight_layout()
    filename = os.path.join(path2save, "alpha_noabx_lineplot_splitcov.png")
    plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.show()


def select_samples_around_nth_abx_exposure(md_df, abx_df, n=1):
    """
    get samples around n-th abx exposure (n=1 is first abx exposure, n=2 is
    second etc.)

    Args:
        md_df (_type_): _description_ abx_df (_type_): _description_ n (int,
        optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
    # indexing starts at zero
    n = n - 1
    # calculate age at n-th abx exposure for all hosts
    abx_nth_age = abx_df.groupby("host_id").nth(n)
    abx_nth_age = abx_nth_age.rename(columns={"abx_start_age_months": "age_nth_abx"})

    # add this column to all_samples
    all_samples = pd.merge(md_df, abx_nth_age, on="host_id", how="left")

    # calculate time of samples since n-th abx exposure
    all_samples = all_samples.assign(
        diff_age_nth_abx=all_samples["age_months_rounded05"]
        - all_samples["age_nth_abx"]
    )
    # round to full months for simplicity. note: added 0.01 since lots of 0.5
    # would otw be rounded down leading to uneven sample distribution
    all_samples["diff_age_nth_abx"] = all_samples["diff_age_nth_abx"] + 0.01
    all_samples["diff_age_nth_abx"] = all_samples["diff_age_nth_abx"].round(0)

    # select only samples before and after nth abx exposure
    abx_nth_samples = all_samples.loc[
        np.logical_and(
            ~all_samples["diff_age_nth_abx"].isna(),
            # really only samples around this n-th exposure
            np.logical_and(
                all_samples["abx_any_cumcount"] <= (n + 1),
                all_samples["abx_any_cumcount"] >= n,
            ),
        ),
        :,
    ]

    # only select samples that are up to 3 months prior to n-th abx exposure and
    # 12 months after
    abx_nth_samples = abx_nth_samples.loc[
        np.logical_and(
            abx_nth_samples["diff_age_nth_abx"] >= -3.0,
            abx_nth_samples["diff_age_nth_abx"] <= 12.0,
        ),
        :,
    ]
    return abx_nth_samples


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


def plot_diversity_difference(
    abx_nth_samples,
    n,
    path_to_output,
    div_metrics=[
        "div_alpha_faith_pd",
        "div_alpha_observed_features",
        "div_alpha_shannon",
    ],
    matched=True,
    test_zero_mean=False,
):
    dic_kw_results = {}
    for metric in div_metrics:
        if matched:
            # ! plot difference in diversity to matched samples
            # drop rows where mean metric is nan (no comparable group exists in noabx)
            nb_dropped = (abx_nth_samples[f"mean_{metric}"].isna()).sum()
            print(
                f"Number of samples disregarded because of lacking reference "
                f"in noabx: {nb_dropped}"
            )
            abx_nth_samples = abx_nth_samples[~abx_nth_samples[f"mean_{metric}"].isna()]
            abx_nth_samples = abx_nth_samples.assign(
                diff_metric=abx_nth_samples[metric] - abx_nth_samples[f"mean_{metric}"]
            )
            plotted_metric = "diff_metric"
            max_ylabel_ratio = 2
        else:
            # ! plot just diversity
            plotted_metric = metric
            max_ylabel_ratio = 1.8
        x_axis = "diff_age_nth_abx"
        abx_nth_samples["diff_age_nth_abx"].replace({-0.0: 0.0}, inplace=True)

        fig, axs = plt.subplots(
            2, 1, figsize=(12, 8), height_ratios=[1, 0.5], sharex=True, dpi=400
        )
        max_score = abx_nth_samples[plotted_metric].max()
        min_score = abx_nth_samples[plotted_metric].min()
        if min_score > 0:
            min_ylabel_ratio = 0.0
        else:
            min_ylabel_ratio = 1.1
        axs[0].set_ylim(min_ylabel_ratio * min_score, max_ylabel_ratio * max_score)
        # Define a custom color palette
        palette = ["grey"] * 3 + ["purple"] * 13
        # category needed to have consistent x-axis
        abx_nth_samples[f"{x_axis}_cat"] = pd.Categorical(
            abx_nth_samples[x_axis], categories=np.arange(-3.0, 13.0)
        )
        sns.boxplot(
            x=f"{x_axis}_cat",
            y=plotted_metric,
            data=abx_nth_samples,
            ax=axs[0],
            palette=palette,
            # legend=False,
        )
        axs[0].legend([], [], frameon=False)
        axs[0].set_xlim([-3, 13])
        axs[0].axvline(3, color="darkred")
        axs[0].axhline(0, color="grey")
        if matched:
            axs[0].set_ylabel("Difference to matched alpha diversity")
        else:
            axs[0].set_ylabel("Alpha diversity")
        axs[0].set_title(metric.replace("div_alpha_", ""))

        t1_values = [x for x in range(-3, 13)]
        result_df = perform_significance_tests(
            abx_nth_samples, -1.0, t1_values, plotted_metric
        )
        dic_kw_results[metric] = result_df

        # Add a star above the boxplots if the p-value < 0.10
        unpaired_color = "sandybrown"
        paired_color = "darkgreen"
        dic_tests = {"unpaired": [1.2, unpaired_color], "paired": [1.1, paired_color]}
        if matched and test_zero_mean:
            # test if mean is significantly different from zero
            t1_val_zerom = [x for x in range(-3, 13)]
            result_df_zerom = test_mean_zero(
                abx_nth_samples, t1_val_zerom, plotted_metric
            )
            result_df = pd.merge(result_df, result_df_zerom, on="t1", how="outer")
            dic_tests["zerom"] = [1.0, "darkblue"]

        for test, (y_shift, color) in dic_tests.items():
            for t1, p_val in zip(result_df.index, result_df[f"P-value {test}"]):
                if p_val < 0.05:
                    sign = "**"
                elif p_val < 0.1:
                    sign = "*"

                if p_val < 0.1:
                    max_y = abx_nth_samples[plotted_metric].max()
                    axs[0].text(
                        t1 + 3,
                        y_shift * max_y,
                        sign,
                        color=color,
                        ha="center",
                        fontsize=25,
                    )

        # add a count barplot for paired and unpaired
        # axs[1] is the barplot
        grouped_counts = (
            abx_nth_samples.groupby(x_axis)[plotted_metric]
            .count()
            .reset_index(name="counts")
        )
        sns.barplot(
            x=x_axis, y="counts", data=grouped_counts, color="peachpuff", ax=axs[1]
        )

        df = dic_kw_results[metric].reset_index()
        for k, v in {
            "unpaired": ["none", unpaired_color],
            "paired": ["none", paired_color],
        }.items():
            df_p = df[["t1", f"# samples {k}"]].copy()
            sns.barplot(
                x=df_p["t1"].astype(float),
                y=df_p[f"# samples {k}"],
                ax=axs[1],
                facecolor=v[0],
                edgecolor=v[1],
            )
        axs[1].axvline(3, color="darkred")
        axs[1].set_ylabel("Number of samples")
        axs[1].set_xlabel(f"Months since {n}. abx exposure")

        # Create a custom legend
        custom_lines = [
            Line2D([0], [0], color=unpaired_color, lw=3),
            Line2D([0], [0], color=paired_color, lw=3),
        ]
        custom_cross = [
            Line2D(
                [0],
                [0],
                color=unpaired_color,
                marker="*",
                markersize=12,
                linestyle="None",
            ),
            Line2D(
                [0],
                [0],
                color=paired_color,
                marker="*",
                markersize=12,
                linestyle="None",
            ),
        ]
        legend_txt = ["unpaired to -1.0", "paired to -1.0"]
        if matched and test_zero_mean:
            custom_cross.append(
                Line2D(
                    [0],
                    [0],
                    color="darkblue",
                    marker="*",
                    markersize=12,
                    linestyle="None",
                ),
            )
            legend_txt.append("diff. to zero mean")
        axs[0].legend(custom_cross, legend_txt, loc="upper right")
        axs[1].legend(custom_lines, legend_txt)

        plt.tight_layout()
        filename = os.path.join(
            path_to_output, f"alpha_abx_boxplot_diff_{n}th_abx_{metric}.png"
        )
        plt.savefig(filename, dpi=400, bbox_inches="tight")
        plt.show()
    return fig, dic_kw_results


def calculate_nth_abx_effect_on_diversity(
    md_df,
    abx_df,
    n,
    path_to_output,
    cov_groups=["age_months", "delivery_mode", "diet_milk_s"],
    div_metrics=[
        "div_alpha_faith_pd",
        "div_alpha_observed_features",
        "div_alpha_shannon",
    ],
    test_zero_mean=False,
):
    # only select samples that are around n-th abx exposure
    abx_nth_samples = select_samples_around_nth_abx_exposure(md_df, abx_df, n=n)

    # also replace nan in cov_group columns with "unknown" -> better matching
    abx_nth_samples[cov_groups] = abx_nth_samples[cov_groups].fillna("unknown")

    # get average diversity metric for noabx infants grouped by covariates
    noabx = md_df.loc[md_df["max_abx_w_microbiome"] == 0.0, :].copy()
    # replace nan in cov_group columns with "unknown" to not loose these
    noabx.loc[:, cov_groups] = noabx[cov_groups].fillna("unknown")

    dic_mean_div_noabx = {}
    for metric in div_metrics:
        dic_mean_div_noabx[metric] = calculate_mean_diversity_per_cov(
            noabx, metric, cov_groups
        )

    # match to mean_div from no_abx & calculate difference
    for div_metric in div_metrics:
        noabx_mean_metric = dic_mean_div_noabx[div_metric]

        abx_nth_samples = pd.merge(
            abx_nth_samples, noabx_mean_metric, on=cov_groups, how="left"
        )

    fig, dic_kw_results = plot_diversity_difference(
        abx_nth_samples,
        n,
        path_to_output,
        div_metrics,
        matched=True,
        test_zero_mean=test_zero_mean,
    )
    return fig, dic_kw_results


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
