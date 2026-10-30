"""
Cohort description figure of the manuscript (Fig. 3): microbiome samples and
antibiotic exposures per infant over age (a), diet covariates of the samples
over age (b) and sample characteristics (c). The figure is drawn at its
printed width with one font size; see scripts/describe_dataset.ipynb.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

from src.color_maps import all_color_maps

# manuscript text width (5.17 in) and one size for all axis text
FIGSIZE = (5.17, 4.3)
BASE_FONT = 7
LETTER_FONT = 10
FIGURE_STEM = "dataset"
SAMPLE_PALETTE = {True: "#f4a461", False: "#014487"}
EVENT_COLOR = "darkred"
UNKNOWN_COLOR = "#d9d9d8"
DIET_ORDER = {
    "diet_weaning": ["no", "yes", "unknown"],
    "diet_milk": ["bd", "mixed", "fd", "unknown"],
}
CHARACTERISTIC_COLORS = {
    "# samples per host": "#ababab",
    "Duration between samples [days]": "#844f7f",
    "Age at first sample [days]": "#6c9cc2",
}
AGE_TICKS = [5, 10, 15, 20, 25, 30, 35]
SUPTITLE = "Distribution of microbial samples over time"
SORT_COLUMNS = [
    "abx_max_count_ever",
    "max_abx_w_microbiome",
    "host_id",
    "age_months_rounded05",
]
# layout in figure fractions: row a (two scatter panels, legend in the right
# margin), row b (left) and row c (right)
SCATTER_GRID = dict(left=0.06, right=0.68, top=0.905, bottom=0.60, wspace=0.12)
DIET_GRID = dict(left=0.075, right=0.53, top=0.47, bottom=0.075, hspace=0.15)
CHARACTERISTIC_GRID = dict(left=0.60, right=0.985, top=0.47, bottom=0.075, hspace=0.95)
LETTER_POSITIONS = {"a": (0.005, 0.995), "b": (0.005, 0.515), "c": (0.545, 0.515)}
SUPTITLE_POSITION = (0.37, 0.975)
LEGEND_KWS = dict(frameon=False, handletextpad=0.4, labelspacing=0.3, borderaxespad=0.0)


def rc_params(base_font=BASE_FONT):
    """Matplotlib rc parameters of the figure at its printed size."""
    return {
        "font.family": "DejaVu Sans",
        "font.size": base_font,
        "axes.titlesize": base_font,
        "axes.labelsize": base_font,
        "xtick.labelsize": base_font,
        "ytick.labelsize": base_font,
        "legend.fontsize": base_font,
        "legend.title_fontsize": base_font,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.pad": 2.0,
        "ytick.major.pad": 2.0,
        "axes.titlepad": 3.0,
        "axes.labelpad": 2.0,
    }


def prepare_samples(md_df, th_months=3):
    """
    Samples with the flag of being collected up to `th_months` after an
    antibiotic exposure and the number of exposures per host before its
    samples, sorted as in the notebook.
    """
    md_df = md_df.assign(
        sample_lt_xm_after_abx=lambda df: df["abx_any_last_t_dmonths"] <= th_months,
        max_abx_w_microbiome=lambda df: df.groupby("host_id")[
            "abx_any_cumcount"
        ].transform("max"),
    )
    return md_df.sort_values(SORT_COLUMNS).reset_index(drop=True)


def split_by_exposure(md_df, abx_events):
    """Samples and exposure events of the hosts without and with exposure."""
    groups = {}
    for name, exposed in (("w/o", False), ("with", True)):
        mask = md_df["max_abx_w_microbiome"] > 0
        samples = md_df[mask if exposed else ~mask].copy()
        events = abx_events[abx_events["host_id"].isin(samples["host_id"].unique())]
        groups[name] = (samples, events.copy())
    return groups


def diet_counts_by_age(md_df, diet_var, age_var="age_months_rounded1"):
    """Sample counts per age bin and diet category (missing = unknown)."""
    values = md_df[diet_var].astype(object).where(md_df[diet_var].notna(), "unknown")
    counts = pd.crosstab(md_df[age_var], values)
    return counts.reindex(columns=DIET_ORDER[diet_var], fill_value=0).sort_index()


def sample_characteristics(md_df):
    """Samples per host, days between consecutive samples and age at first sample."""
    by_host = md_df.sort_values(["host_id", "age_days"]).groupby("host_id")["age_days"]
    return {
        "# samples per host": by_host.count(),
        "Duration between samples [days]": by_host.diff().dropna(),
        "Age at first sample [days]": by_host.min(),
    }


def diet_colors(diet_var):
    """Colours of the diet categories in their plotting order."""
    colors = {**all_color_maps[diet_var], "unknown": UNKNOWN_COLOR}
    return [colors[category] for category in DIET_ORDER[diet_var]]


def sample_legend_handles(th_months):
    """Handles and labels of the legend of the sample panels."""
    handles = [
        Line2D([0], [0], color=EVENT_COLOR, marker="x", ls="None", markersize=4),
        Line2D(
            [0], [0], color=SAMPLE_PALETTE[True], marker="o", ls="None", markersize=3
        ),
        Line2D(
            [0], [0], color=SAMPLE_PALETTE[False], marker="o", ls="None", markersize=3
        ),
    ]
    labels = [
        "abx event",
        f"sample <= {th_months}m after abx",
        f"sample > {th_months}m after abx",
    ]
    return handles, labels


def _draw_samples(ax, samples, events, title, ylabel, markersize=3):
    sns.scatterplot(
        x="age_months_rounded05",
        y="host_id",
        hue="sample_lt_xm_after_abx",
        palette=SAMPLE_PALETTE,
        data=samples,
        ax=ax,
        s=markersize,
        linewidth=0,
        legend=False,
    )
    sns.scatterplot(
        x="abx_start_age_months",
        y="host_id",
        data=events,
        ax=ax,
        s=markersize * 1.5,
        marker="x",
        color=EVENT_COLOR,
        linewidth=0.5,
        legend=False,
    )
    ax.set_title(title)
    ax.set_xlabel("Age [months]")
    ax.set_ylabel(ylabel)
    ax.margins(y=0.005)
    ax.tick_params(left=False, labelleft=False)


def _draw_diet_bars(ax, counts, diet_var):
    bottom = np.zeros(len(counts))
    for category, color in zip(counts.columns, diet_colors(diet_var)):
        ax.bar(
            counts.index,
            counts[category],
            bottom=bottom,
            color=color,
            width=0.8,
            label=category,
            lw=0,
        )
        bottom += counts[category].values
    ax.set_xticks(AGE_TICKS)
    ax.set_ylabel("# samples")
    ax.legend(
        title=diet_var.replace("_", " "),
        loc="upper right",
        frameon=True,
        framealpha=0.9,
        edgecolor="none",
        handlelength=1.0,
        handleheight=0.8,
        labelspacing=0.2,
        borderpad=0.3,
    )


def _draw_violin_box(ax, values, color, title):
    sns.violinplot(x=values, inner=None, color=color, ax=ax, orient="h", linewidth=0.6)
    sns.boxplot(
        x=values,
        width=0.1,
        boxprops={"facecolor": "white", "edgecolor": "black", "zorder": 2},
        flierprops={
            "marker": "o",
            "markerfacecolor": "none",
            "markeredgecolor": "black",
            "markersize": 2,
        },
        linewidth=0.6,
        ax=ax,
        orient="h",
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xlim(left=0)
    ax.tick_params(left=False, labelleft=False)


def build_dataset_figure(
    md_df, abx_events, th_months=3, base_font=BASE_FONT, figsize=FIGSIZE
):
    """
    Composes the figure: samples and exposures per host (a), diet covariates
    over age (b) and sample characteristics (c).

    Args:
        md_df (pd.DataFrame): Sample metadata of the feature table.
        abx_events (pd.DataFrame): Start age (months) of each exposure per host.
        th_months (float): Samples up to this many months after an exposure
        are highlighted.

    Returns:
        matplotlib.figure.Figure
    """
    md_df = prepare_samples(md_df, th_months)
    groups = split_by_exposure(md_df, abx_events)
    with plt.rc_context(rc_params(base_font)):
        fig = plt.figure(figsize=figsize, dpi=300)

        scatter_grid = fig.add_gridspec(1, 2, **SCATTER_GRID)
        ylabels = {
            "w/o": r"Host ID$_{\mathrm{no\ abx}}$",
            "with": r"Host ID$_{\mathrm{abx}}$",
        }
        scatter_axes = []
        for i, (name, (samples, events)) in enumerate(groups.items()):
            ax = fig.add_subplot(
                scatter_grid[0, i], sharex=scatter_axes[0] if scatter_axes else None
            )
            title = f"Hosts {name} abx exposure ({samples['host_id'].nunique()})"
            _draw_samples(ax, samples, events, title, ylabels[name])
            scatter_axes.append(ax)
        handles, labels = sample_legend_handles(th_months)
        scatter_axes[-1].legend(
            handles, labels, loc="upper left", bbox_to_anchor=(1.04, 1.0), **LEGEND_KWS
        )
        fig.text(
            *SUPTITLE_POSITION, SUPTITLE, ha="center", va="top", fontsize=base_font
        )

        diet_grid = fig.add_gridspec(2, 1, **DIET_GRID)
        diet_axes = []
        for i, diet_var in enumerate(DIET_ORDER):
            ax = fig.add_subplot(
                diet_grid[i, 0], sharex=diet_axes[0] if diet_axes else None
            )
            _draw_diet_bars(ax, diet_counts_by_age(md_df, diet_var), diet_var)
            diet_axes.append(ax)
        diet_axes[0].tick_params(labelbottom=False)
        diet_axes[-1].set_xlabel("Age [months]")
        fig.align_ylabels(diet_axes)

        characteristic_grid = fig.add_gridspec(3, 1, **CHARACTERISTIC_GRID)
        for i, (title, values) in enumerate(sample_characteristics(md_df).items()):
            ax = fig.add_subplot(characteristic_grid[i, 0])
            _draw_violin_box(ax, values, CHARACTERISTIC_COLORS[title], title)

        for letter, (x, y) in LETTER_POSITIONS.items():
            fig.text(
                x,
                y,
                letter,
                fontsize=LETTER_FONT,
                fontweight="bold",
                ha="left",
                va="top",
            )
    return fig


def save_dataset_figure(fig, path2save, stem=FIGURE_STEM):
    """Saves the figure as PDF (and PNG preview) at its page size;
    returns the PDF path."""
    path2save = Path(path2save)
    pdf = path2save / f"{stem}.pdf"
    with plt.rc_context(rc_params()):
        fig.savefig(pdf, format="pdf")
        fig.savefig(path2save / f"{stem}.png", format="png", dpi=200)
    return pdf
