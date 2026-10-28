"""Sensitivity of the post-antibiotic diversity deficit to the alpha diversity metric.

Repeats the covariate-matched comparison of Fig. 4c of the manuscript for Faith PD and
Shannon diversity on the final feature table (``data/final``). For every sample of an
antibiotic-exposed infant the deficit relative to unexposed infants is

    deficit = mean(unexposed, matched stratum) - observed value,

with strata defined by monthly age bin, delivery mode, milk diet and weaning status
(missing covariates form their own stratum). To compare the metrics on one axis, each
deficit is expressed in percent of the matched mean of its metric (relative deficit).
Samples are then placed in monthly bins
relative to the first, second and third antibiotic exposure exactly as in Fig. 4c:
window of -3 to +6 months, only samples with the exposure count of that window, last
sample per infant and bin, and all pre-exposure samples collapsed to the last sample
before the exposure (bin -1). Exposures are limited to the first 24 months.

Outputs (``--out-dir``): the per-sample windows, the per-bin tests (deficit vs. the
pre-exposure bin per metric, Shannon vs. Faith PD within bins, Spearman concordance of
post-exposure deficits) and the figure.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
    }
)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "data" / "final"
FT_DEFAULT = DATA / "ft_vat19_anomaly_v20240806_entero_genus.tsv"
ABX_DEFAULT = DATA / "ts_vat19_abx_v20240806.tsv"

DAYS_PER_MONTH = 30.437
LIMIT_MONTHS = 24.0
MIN_MONTHS, MAX_MONTHS = -3.0, 6.0
PRE_BIN = -1.0
MATCHING_VARS = ["month_bin", "delivery_mode", "diet_milk", "diet_weaning"]
METRICS = {"div_alpha_faith_pd": "Faith PD", "div_alpha_shannon": "Shannon"}
COLORS = {"div_alpha_faith_pd": "#0072B2", "div_alpha_shannon": "#D55E00"}
METRIC, REFERENCE = "div_alpha_shannon", "div_alpha_faith_pd"
COMPARISON = "shannon_vs_faith_pd"
EXPOSURES = (1, 2, 3)
SUFFIX = {1: "st", 2: "nd", 3: "rd"}


def add_month_bins(df: pd.DataFrame) -> pd.DataFrame:
    """Monthly and half-monthly age bins as used by the modelling repository."""
    df = df.copy()
    df["day"] = df["age_days"].astype(int)
    df["month_bin"] = (df["day"] / DAYS_PER_MONTH).round().astype(int)
    df["month5_bin"] = (df["day"] / DAYS_PER_MONTH * 2).round() / 2
    return df


def load_feature_table(path: Path = FT_DEFAULT) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "sample_id"
    df = add_month_bins(df.reset_index())
    df["unexposed"] = df.groupby("host_id")["abx_any_cumcount"].transform("max") == 0
    df[MATCHING_VARS[1:]] = df[MATCHING_VARS[1:]].fillna("unknown")
    return (
        df[df["month5_bin"] <= LIMIT_MONTHS]
        .sort_values(["host_id", "day"])
        .reset_index(drop=True)
    )


def load_abx(path: Path = ABX_DEFAULT) -> pd.DataFrame:
    abx = pd.read_csv(path, sep="\t")[["host_id", "abx_start_age_months"]]
    # one start age is given as 7.6 instead of 7.5 in the original supplementary data
    abx.loc[abx["host_id"] == "E014403", "abx_start_age_months"] = 7.5
    abx = abx[abx["abx_start_age_months"] <= LIMIT_MONTHS]
    return abx.sort_values(["host_id", "abx_start_age_months"]).reset_index(drop=True)


def add_matched_deficits(ft: pd.DataFrame, metrics=tuple(METRICS)) -> pd.DataFrame:
    """Absolute and percent deficit of exposed samples vs. matched unexposed infants."""
    ref = ft[ft["unexposed"]]
    means = ref.groupby(MATCHING_VARS)[list(metrics)].mean()
    means.columns = [f"mean_{m}" for m in metrics]
    exposed = ft[~ft["unexposed"]].merge(means, on=MATCHING_VARS, how="left")
    for m in metrics:
        exposed[f"delta_{m}"] = exposed[f"mean_{m}"] - exposed[m]
        exposed[f"rel_{m}"] = 100 * exposed[f"delta_{m}"] / exposed[f"mean_{m}"]
    return exposed


def select_around_exposure(
    exposed: pd.DataFrame, abx: pd.DataFrame, n: int, value_cols: list[str]
) -> pd.DataFrame:
    """Samples in monthly bins around the n-th exposure, as in Fig. 4c."""
    nth = abx.groupby("host_id").nth(n - 1)[["host_id", "abx_start_age_months"]]
    nth = nth.rename(columns={"abx_start_age_months": "age_nth_abx"})
    df = exposed.merge(nth, on="host_id", how="inner")
    df["months_since"] = df["month5_bin"] - df["age_nth_abx"]
    df = df[df["months_since"].between(MIN_MONTHS, MAX_MONTHS)]
    df = df[df["abx_any_cumcount"].between(n - 1, n)].dropna(subset=value_cols)
    edges = np.arange(MIN_MONTHS, MAX_MONTHS + 1.0, 1.0)
    df = df.assign(
        bin=pd.cut(
            df["months_since"], bins=edges, labels=edges[:-1], right=False
        ).astype(float)
    )
    df = df.sort_values(["host_id", "day"]).groupby(["host_id", "bin"]).tail(1)
    pre = (
        df[df["bin"] < 0]
        .sort_values("bin")
        .groupby("host_id")
        .tail(1)
        .assign(bin=PRE_BIN)
    )
    return (
        pd.concat([pre, df[df["bin"] >= 0]])
        .sort_values(["host_id", "day"])
        .assign(exposure=n)
    )


def tests_vs_pre(window: pd.DataFrame, col: str) -> pd.DataFrame:
    """One-sided tests of a larger deficit in each post-exposure bin than in the
    pre-exposure bin."""
    pre = window.loc[window["bin"] == PRE_BIN, ["host_id", col]].rename(
        columns={col: "pre"}
    )
    rows = []
    for b, g in window[window["bin"] >= 0].groupby("bin"):
        p_mwu = (
            mannwhitneyu(pre["pre"], g[col], alternative="less").pvalue
            if len(pre) and len(g)
            else np.nan
        )
        paired = pre.merge(g[["host_id", col]], on="host_id")
        p_wil = (
            wilcoxon(paired["pre"], paired[col], alternative="less").pvalue
            if len(paired)
            else np.nan
        )
        rows.append(
            dict(
                bin=b,
                n=len(g),
                median=g[col].median(),
                mean=g[col].mean(),
                p_mwu_vs_pre=p_mwu,
                n_paired=len(paired),
                p_wilcoxon_vs_pre=p_wil,
            )
        )
    return pd.DataFrame(rows)


def cross_metric_tests(window: pd.DataFrame) -> pd.DataFrame:
    """Two-sided paired test of the relative deficit, Shannon vs. Faith PD, per bin."""
    rows = []
    for b, g in window.groupby("bin"):
        diff = g[f"rel_{METRIC}"] - g[f"rel_{REFERENCE}"]
        p = wilcoxon(diff).pvalue if len(diff) and (diff != 0).any() else np.nan
        rows.append(
            dict(
                bin=b,
                comparison=COMPARISON,
                n=len(g),
                median_difference=diff.median(),
                p_wilcoxon=p,
            )
        )
    return pd.DataFrame(rows)


def post_exposure_samples(windows: dict[int, pd.DataFrame]) -> pd.DataFrame:
    """Unique post-exposure samples (bins >= 0) pooled over exposures."""
    post = pd.concat([w[w["bin"] >= 0] for w in windows.values()])
    return post.drop_duplicates("sample_id")


def concordance(windows: dict[int, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    groups = {str(n): w[w["bin"] >= 0] for n, w in windows.items()}
    groups["pooled"] = post_exposure_samples(windows)
    for label, g in groups.items():
        rho = (
            spearmanr(g[f"rel_{METRIC}"], g[f"rel_{REFERENCE}"]).statistic
            if len(g) > 2
            else np.nan
        )
        rows.append(
            dict(exposure=label, comparison=COMPARISON, n=len(g), spearman_rho=rho)
        )
    return pd.DataFrame(rows)


def run(ft_path: Path = FT_DEFAULT, abx_path: Path = ABX_DEFAULT) -> dict:
    ft = load_feature_table(ft_path)
    exposed = add_matched_deficits(ft)
    rel_cols = [f"rel_{m}" for m in METRICS]
    windows = {
        n: select_around_exposure(exposed, load_abx(abx_path), n, rel_cols)
        for n in EXPOSURES
    }
    vs_pre = pd.concat(
        [
            tests_vs_pre(w, f"rel_{m}").assign(exposure=n, metric=METRICS[m])
            for n, w in windows.items()
            for m in METRICS
        ]
    )
    cross = pd.concat(
        [cross_metric_tests(w).assign(exposure=n) for n, w in windows.items()]
    )
    return dict(
        ft=ft,
        windows=windows,
        vs_pre=vs_pre,
        cross=cross,
        concordance=concordance(windows),
    )


def _stars(p: float) -> str:
    return "**" if p < 0.05 else "*" if p < 0.1 else ""


def _draw_boxes(
    ax, window: pd.DataFrame, vs_pre: pd.DataFrame, n: int, star_rows: dict[str, float]
) -> None:
    bins = [PRE_BIN] + list(np.arange(0.0, MAX_MONTHS, 1.0))
    width, offsets = 0.34, (-0.2, 0.2)
    for (m, label), off in zip(METRICS.items(), offsets):
        data = [window.loc[window["bin"] == b, f"rel_{m}"].values for b in bins]
        pos = np.arange(len(bins)) + off
        ax.boxplot(
            data,
            positions=pos,
            widths=width,
            patch_artist=True,
            showfliers=True,
            boxprops=dict(
                facecolor=COLORS[m], alpha=0.55, edgecolor=COLORS[m], linewidth=0.9
            ),
            medianprops=dict(color="black", linewidth=1.0),
            whiskerprops=dict(color=COLORS[m], linewidth=0.9),
            capprops=dict(color=COLORS[m], linewidth=0.9),
            flierprops=dict(
                marker="o",
                markersize=2.2,
                markerfacecolor="none",
                markeredgecolor=COLORS[m],
                alpha=0.7,
            ),
        )
        sub = vs_pre[(vs_pre["exposure"] == n) & (vs_pre["metric"] == label)].set_index(
            "bin"
        )
        for i, b in enumerate(bins):
            if b in sub.index and _stars(sub.loc[b, "p_mwu_vs_pre"]):
                ax.text(
                    pos[i],
                    star_rows[m],
                    _stars(sub.loc[b, "p_mwu_vs_pre"]),
                    color=COLORS[m],
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                )
    ax.set_xticks(np.arange(len(bins)))
    ax.set_xticklabels([f"{int(b)}" for b in bins], fontsize=7.5)
    ax.axhline(0, color="grey", linewidth=0.8, zorder=0)
    ax.axvline(0.5, color="darkred", linewidth=1.0)
    ax.set_xlabel(f"Months since {n}{SUFFIX[n]} antibiotic exposure", fontsize=8.5)
    ax.set_title(f"{n}{SUFFIX[n]} exposure", fontsize=9.5)
    ax.tick_params(axis="y", labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _draw_scatter(ax, post: pd.DataFrame, conc: pd.DataFrame) -> None:
    x, y = post[f"rel_{REFERENCE}"], post[f"rel_{METRIC}"]
    ax.scatter(x, y, s=9, color=COLORS[METRIC], alpha=0.55, linewidths=0)
    lim = [min(x.min(), y.min()) - 5, max(x.max(), y.max()) + 5]
    ax.plot(lim, lim, color="grey", linewidth=0.8, linestyle="--", zorder=0)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    row = conc[conc["exposure"] == "pooled"].iloc[0]
    ax.text(
        0.03,
        0.96,
        f"Spearman $\\rho$ = {row['spearman_rho']:.2f}\nn = {int(row['n'])} samples",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.set_xlabel(f"Relative deficit, {METRICS[REFERENCE]} [%]", fontsize=8.5)
    ax.set_ylabel(f"Relative deficit, {METRICS[METRIC]} [%]", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def build_figure(res: dict, out_pdf: Path, out_png: Path | None = None) -> None:
    fig = plt.figure(figsize=(9.5, 7.0), dpi=400)
    grid = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.95], hspace=0.5, wspace=1.1)
    top = [fig.add_subplot(grid[0, 2 * i : 2 * i + 2]) for i in range(3)]
    all_rel = pd.concat(res["windows"].values())[[f"rel_{m}" for m in METRICS]]
    rel_min, rel_max = all_rel.min().min(), all_rel.max().max()
    span = rel_max - rel_min
    star_rows = {m: rel_max + span * (0.05 + 0.06 * i) for i, m in enumerate(METRICS)}
    for ax, n in zip(top, EXPOSURES):
        _draw_boxes(ax, res["windows"][n], res["vs_pre"], n, star_rows)
        ax.set_ylim(rel_min - 0.04 * span, rel_max + 0.17 * span)
    top[0].set_ylabel(
        "Deficit vs. matched unexposed infants\n[% of matched mean]", fontsize=8.5
    )
    handles = [
        Patch(facecolor=COLORS[m], alpha=0.55, edgecolor=COLORS[m], label=label)
        for m, label in METRICS.items()
    ]
    handles += [
        Line2D([0], [0], color="darkred", linewidth=1.0, label="antibiotic exposure"),
        Line2D(
            [0],
            [0],
            marker="$*$",
            color="black",
            linestyle="None",
            markersize=7,
            label="deficit larger than in bin -1 (* p<0.1, ** p<0.05)",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        fontsize=7.5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.955),
    )
    bottom = fig.add_subplot(grid[1, 1:5])
    _draw_scatter(bottom, post_exposure_samples(res["windows"]), res["concordance"])
    for ax, letter in zip(top + [bottom], "abcd"):
        ax.text(
            -0.14,
            1.08,
            letter,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="bottom",
        )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.08)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    if out_png is not None:
        fig.savefig(out_png, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_outputs(res: dict, out_dir: Path, out_png: Path | None = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    keep = [
        "sample_id",
        "host_id",
        "exposure",
        "bin",
        "day",
        "month_bin",
        "abx_any_cumcount",
    ] + MATCHING_VARS[1:]
    keep += [c for m in METRICS for c in (m, f"mean_{m}", f"delta_{m}", f"rel_{m}")]
    pd.concat(res["windows"].values())[keep].to_csv(
        out_dir / "metric_sensitivity_windows.csv", index=False
    )
    res["vs_pre"].to_csv(out_dir / "metric_sensitivity_tests_vs_pre.csv", index=False)
    res["cross"].to_csv(
        out_dir / "metric_sensitivity_tests_between_metrics.csv", index=False
    )
    res["concordance"].to_csv(
        out_dir / "metric_sensitivity_concordance.csv", index=False
    )
    build_figure(res, out_dir / "supp_metric_sensitivity.pdf", out_png)


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ft", type=Path, default=FT_DEFAULT)
    parser.add_argument("--abx", type=Path, default=ABX_DEFAULT)
    parser.add_argument("--out-dir", type=Path, default=HERE)
    parser.add_argument("--out-png", type=Path, default=None)
    args = parser.parse_args(argv)
    res = run(args.ft, args.abx)
    write_outputs(res, args.out_dir, args.out_png)
    n_unexposed = res["ft"].loc[res["ft"]["unexposed"], "host_id"].nunique()
    print(f"unexposed infants: {n_unexposed}")
    print(res["concordance"].to_string(index=False))
    sig = res["cross"][res["cross"]["p_wilcoxon"] < 0.1]
    print(
        sig.to_string(index=False)
        if len(sig)
        else "no Shannon vs. Faith PD difference at p<0.1 in any bin"
    )
    return res


if __name__ == "__main__":
    main()
