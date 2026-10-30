import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import AutoLocator
from scipy.stats import mannwhitneyu, wilcoxon

from src.alpha_div_stats import format_p

warnings.simplefilter(action="ignore", category=FutureWarning)
plt.rcParams.update(
    {"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42}
)
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


FIGSIZE = (5.17, 5.5)
FIGURE_STEM = "alpha_all_lineplot_splitcov"
LETTERS = "abcdefgh"


def group_trajectory_stats(df, metric, x_axis, group):
    """Mean, standard deviation and count per age bin and group; the standard
    deviation is NaN for bins with a single sample."""
    stats = df.groupby([x_axis, group])[metric].agg(["mean", "std", "count"])
    return stats.reset_index()


def level_colors(levels, reference):
    """Colors of the seaborn husl palette assigned in sorted level order; the
    returned dictionary lists the reference level first (legend order)."""
    ordered = sorted(set(levels))
    palette = dict(zip(ordered, sns.color_palette("husl", len(ordered)).as_hex()))
    legend_order = [reference] + [lvl for lvl in ordered if lvl != reference]
    return {lvl: palette[lvl] for lvl in legend_order}


def _style_axis(ax, base_font):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(labelsize=base_font, length=2, width=0.6, pad=1.5)


def _draw_trajectories(ax, panel, metric, x_axis, base_font):
    traj = group_trajectory_stats(panel.data, metric, x_axis, panel.group)
    colors = level_colors(traj[panel.group], panel.reference)
    for level, color in colors.items():
        g = traj[traj[panel.group] == level].sort_values(x_axis)
        ax.plot(g[x_axis], g["mean"], color=color, lw=1.1, label=panel.label(level))
        ax.fill_between(
            g[x_axis],
            g["mean"] - g["std"],
            g["mean"] + g["std"],
            color=color,
            alpha=0.2,
            lw=0,
        )
    ax.legend(
        loc="lower right",
        fontsize=base_font - 0.5,
        frameon=False,
        handlelength=1.2,
        labelspacing=0.2,
        borderpad=0.2,
        borderaxespad=0.3,
    )
    return traj, colors


def _draw_counts(ax, traj, colors, x_axis, group, max_age):
    counts = (
        traj.pivot(index=x_axis, columns=group, values="count")
        .reindex(np.arange(1, max_age + 1))
        .fillna(0)
    )
    bottom = np.zeros(len(counts))
    for level, color in colors.items():
        if level not in counts:
            continue
        ax.bar(counts.index, counts[level], bottom=bottom, width=0.8, color=color, lw=0)
        bottom += counts[level].values


def panels_by_p(panels, stats):
    """Panels ordered by increasing joint p-value of their contrasts."""
    p_joint = stats.groupby("panel")["p_joint"].first()
    return sorted(panels, key=lambda panel: p_joint[panel.title])


def latex_pct(estimate, low, high):
    """E.g. '$-9$~\\% [$-13$, $-5$]' (rounded, no negative zero)."""
    est, lo, hi = (0.0 + round(v) for v in (estimate, low, high))
    return f"${est:+.0f}$~\\% [${lo:.0f}$, ${hi:.0f}$]"


def latex_p(p):
    text = format_p(p)
    return f"${text}$" if text.startswith("<") else text


def latex_table(panels, stats):
    """LaTeX tabular of the age-adjusted differences, variables ordered by
    increasing joint p-value; the variable name and the p-value only on the
    first row of each variable."""
    lines = [
        r"\begin{tabular}{l l l r}",
        r"\toprule",
        r"Variable & Contrast & Difference [95\% CI] & $p$ \\",
        r"\midrule",
    ]
    for panel in panels_by_p(panels, stats):
        rows = stats[stats["panel"] == panel.title].reset_index(drop=True)
        for j, row in rows.iterrows():
            cells = [
                panel.title if j == 0 else "",
                f"{panel.label(row['level'])} vs {panel.label(row['reference'])}",
                latex_pct(row["pct_diff"], row["ci_low"], row["ci_high"]),
                latex_p(row["p_joint"]) if j == 0 else "",
            ]
            lines.append(" & ".join(cells).rstrip() + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def write_latex_table(panels, stats, path):
    """Write the LaTeX table of age-adjusted differences and return its path."""
    path = Path(path)
    path.write_text(latex_table(panels, stats))
    return path


def build_group_figure(panels, metric, x_axis, base_font=7.0, ncols=2, figsize=FIGSIZE):
    """Mean trajectories with standard deviation bands and monthly sample counts
    per panel."""
    nrows = -(-len(panels) // ncols)
    fig = plt.figure(figsize=figsize, dpi=300)
    outer = fig.add_gridspec(
        nrows,
        ncols,
        left=0.085,
        right=0.995,
        top=0.965,
        bottom=0.055,
        hspace=0.5,
        wspace=0.3,
    )
    max_age = int(max(p.data[x_axis].max() for p in panels))
    ymax = 0.0
    line_axes, count_axes = [], []
    for i, panel in enumerate(panels):
        inner = outer[divmod(i, ncols)].subgridspec(
            2, 1, height_ratios=[1.0, 0.33], hspace=0.1
        )
        ax = fig.add_subplot(inner[0])
        axc = fig.add_subplot(inner[1], sharex=ax)
        traj, colors = _draw_trajectories(ax, panel, metric, x_axis, base_font)
        _draw_counts(axc, traj, colors, x_axis, panel.group, max_age)
        ymax = max(ymax, float((traj["mean"] + traj["std"].fillna(0)).max()))
        ax.set_title(panel.title, fontsize=base_font + 1, loc="left", pad=3)
        ax.text(
            -0.19,
            1.03,
            LETTERS[i],
            transform=ax.transAxes,
            fontsize=base_font + 2,
            fontweight="bold",
            va="bottom",
            ha="left",
        )
        ax.set_ylabel("Faith PD", fontsize=base_font)
        axc.set_ylabel("# samples", fontsize=base_font)
        axc.set_xlabel("Age [months]", fontsize=base_font, labelpad=1)
        ax.tick_params(labelbottom=False)
        for a in (ax, axc):
            _style_axis(a, base_font)
        line_axes.append(ax)
        count_axes.append(axc)
    for ax in line_axes:
        ax.set_xlim(0.5, max_age + 0.5)
        ax.set_xticks(np.arange(0, max_age + 1, 6))
        ax.set_ylim(0, ymax * 1.04)
    for col in range(ncols):
        fig.align_ylabels(
            [a for i, a in enumerate(line_axes) if i % ncols == col]
            + [a for i, a in enumerate(count_axes) if i % ncols == col]
        )
    return fig


def lineplot_alpha_by_groups(panels, metric, x_axis, path2save, stem=FIGURE_STEM):
    """Save the grouped trajectory figure as PDF (and PNG) and return the PDF path."""
    fig = build_group_figure(panels, metric, x_axis)
    path2save = Path(path2save)
    pdf = path2save / f"{stem}.pdf"
    fig.savefig(pdf, format="pdf", bbox_inches="tight")
    fig.savefig(path2save / f"{stem}.png", format="png", dpi=200, bbox_inches="tight")
    print(pdf)
    return pdf


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
