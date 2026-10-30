"""Tests for the metric sensitivity analysis: binning, matching, windows, outputs."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import metric_sensitivity as ms

HERE = Path(__file__).resolve().parent
REAL_DATA = ms.FT_DEFAULT.exists() and ms.ABX_DEFAULT.exists()


def _ft(rows):
    cols = [
        "sample_id",
        "host_id",
        "age_days",
        "abx_any_cumcount",
        "delivery_mode",
        "diet_milk",
        "diet_weaning",
        "div_alpha_faith_pd",
        "div_alpha_shannon",
        "div_alpha_observed_features",
    ]
    return pd.DataFrame(rows, columns=cols)


def test_month_bins_match_modelling_convention():
    df = ms.add_month_bins(pd.DataFrame({"age_days": [45, 0, 730]}))
    assert df["month_bin"].tolist() == [1, 0, 24]
    assert df["month5_bin"].tolist() == [1.5, 0.0, 24.0]


def test_unexposed_flag_and_unknown_covariates(tmp_path):
    ft = _ft(
        [
            ["s1", "A", 30, 0, "vaginal", None, "no", 1, 1, 1],
            ["s2", "A", 60, 1, "vaginal", "bd", None, 1, 1, 1],
            ["s3", "B", 30, 0, "vaginal", "bd", "no", 1, 1, 1],
            ["s4", "C", 800, 0, "vaginal", "bd", "no", 1, 1, 1],
        ]
    )
    path = tmp_path / "ft.tsv"
    ft.set_index("sample_id").to_csv(path, sep="\t")
    loaded = ms.load_feature_table(path)
    assert loaded.set_index("sample_id")["unexposed"].to_dict() == {
        "s1": False,
        "s2": False,
        "s3": True,
    }
    assert loaded.loc[loaded["sample_id"] == "s1", "diet_milk"].item() == "unknown"
    assert loaded.loc[loaded["sample_id"] == "s2", "diet_weaning"].item() == "unknown"


def test_matched_deficit_relative_to_stratum_mean():
    rows = [
        [f"u{i}", f"U{i}", 30, 0, "vaginal", "bd", "no", v, v, 2 * v]
        for i, v in enumerate([1.0, 2.0, 3.0, 4.0, 5.0])
    ]
    rows += [
        ["u5", "U5", 30, 0, "cesarean", "bd", "no", 10.0, 10.0, 10.0]
    ]  # other stratum, same age bin
    rows += [
        ["e1", "E1", 30, 1, "vaginal", "bd", "no", 1.0, 1.5, 1.0],
        [
            "e2",
            "E2",
            30,
            1,
            "cesarean",
            "fd",
            "no",
            1.0,
            1.0,
            1.0,
        ],  # no matching stratum
        ["e3", "E3", 300, 1, "vaginal", "bd", "no", 1.0, 1.0, 1.0],
    ]  # no unexposed sample at that age
    ft = ms.add_month_bins(_ft(rows))
    ft["unexposed"] = ft["abx_any_cumcount"] == 0
    out = ms.add_matched_deficits(ft).set_index("sample_id")
    assert out.loc["e1", "delta_div_alpha_faith_pd"] == pytest.approx(3.0 - 1.0)
    # deficit in percent of the matched mean of the same metric
    assert out.loc["e1", "rel_div_alpha_faith_pd"] == pytest.approx(100 * 2.0 / 3.0)
    assert out.loc["e1", "rel_div_alpha_shannon"] == pytest.approx(100 * 1.5 / 3.0)
    assert np.isnan(out.loc["e2", "rel_div_alpha_faith_pd"])
    assert np.isnan(out.loc["e3", "rel_div_alpha_faith_pd"])


def test_window_selection_mirrors_fig4c():
    months = [7.5, 9.0, 9.5, 10.0, 10.5, 12.0, 17.0]
    cum = [0, 0, 0, 1, 1, 1, 2]
    rows = [
        [
            f"s{i}",
            "H",
            int(round(m * ms.DAYS_PER_MONTH)),
            c,
            "vaginal",
            "bd",
            "no",
            1.0,
            1.0,
            1.0,
        ]
        for i, (m, c) in enumerate(zip(months, cum))
    ]
    exposed = ms.add_month_bins(_ft(rows))
    exposed["rel_x"] = 1.0
    abx = pd.DataFrame({"host_id": ["H", "H"], "abx_start_age_months": [10.0, 15.0]})
    win = ms.select_around_exposure(exposed, abx, 1, ["rel_x"]).set_index("sample_id")
    # pre-exposure samples collapse to the last one before the exposure; the two
    # samples in bin 0 keep the later one; the sample at 17 months lies after the
    # 2nd exposure
    assert win["bin"].to_dict() == {"s2": -1.0, "s4": 0.0, "s5": 2.0}
    assert (win["exposure"] == 1).all()
    # 2nd exposure at 15 months: the sample 3 months before enters bin -3 and
    # collapses to -1, the samples around the 1st exposure lie outside the window
    win2 = ms.select_around_exposure(exposed, abx, 2, ["rel_x"]).set_index("sample_id")
    assert win2["bin"].to_dict() == {"s5": -1.0, "s6": 2.0}


def test_tests_vs_pre_direction():
    win = pd.DataFrame(
        {
            "host_id": list("abcdefgh") * 2,
            "bin": [-1.0] * 8 + [0.0] * 8,
            "rel": list(range(8)) + [v + 5 for v in range(8)],
        }
    )
    res = ms.tests_vs_pre(win, "rel").set_index("bin")
    assert res.loc[0.0, "n"] == 8 and res.loc[0.0, "n_paired"] == 8
    assert (
        res.loc[0.0, "p_mwu_vs_pre"] < 0.05 and res.loc[0.0, "p_wilcoxon_vs_pre"] < 0.05
    )


def test_cross_metric_tests_direction_and_identical_metrics():
    win = pd.DataFrame(
        {
            "bin": [0.0] * 6,
            "rel_div_alpha_faith_pd": [1, 2, 3, 4, 5, 6.0],
            "rel_div_alpha_shannon": [2, 3, 4, 5, 6, 7.0],
        }
    )
    res = ms.cross_metric_tests(win).set_index("bin")
    assert (
        res.loc[0.0, "median_difference"] == 1.0 and res.loc[0.0, "p_wilcoxon"] < 0.05
    )
    same = ms.cross_metric_tests(
        win.assign(rel_div_alpha_shannon=win["rel_div_alpha_faith_pd"])
    )
    assert same["p_wilcoxon"].isna().all()


@pytest.mark.skipif(not REAL_DATA, reason="final feature table not available")
def test_real_data_outputs(tmp_path):
    res = ms.run()
    assert res["ft"].loc[res["ft"]["unexposed"], "host_id"].nunique() == 140
    for n, win in res["windows"].items():
        counts = win.groupby("bin").size()
        assert counts.loc[-1.0] > 0 and counts.loc[0.0] > 0
        assert not win[[f"rel_{m}" for m in ms.METRICS]].isna().any().any()
        assert win.groupby(["host_id", "bin"]).size().max() == 1
    ms.write_outputs(res, tmp_path, tmp_path / "preview.png")
    for name in [
        "supp_metric_sensitivity.pdf",
        "preview.png",
        "metric_sensitivity_windows.csv",
        "metric_sensitivity_tests_vs_pre.csv",
        "metric_sensitivity_tests_between_metrics.csv",
        "metric_sensitivity_concordance.csv",
    ]:
        assert (tmp_path / name).stat().st_size > 0
    pooled = res["concordance"].set_index("exposure")
    assert pooled.loc["pooled", "spearman_rho"] > 0.5
    assert set(res["ft"].columns) >= set(ms.METRICS) and len(ms.METRICS) == 2


@pytest.mark.skipif(not REAL_DATA, reason="final feature table not available")
def test_legend_below_the_box_plot_row_above_the_scatter():
    fig = ms.compose_figure(ms.run())
    fig.canvas.draw()
    (legend,) = fig.legends
    box = legend.get_window_extent()
    top, scatter = fig.axes[:3], fig.axes[3]
    # below the x labels of the box plots, which carry the stars
    assert box.y1 <= min(ax.xaxis.label.get_window_extent().y0 for ax in top)
    # above the scatter panel and its letter
    letter = [t for t in scatter.texts if t.get_text() == "d"][0]
    assert box.y0 >= letter.get_window_extent().y1
    centre = (top[0].get_window_extent().x0 + top[-1].get_window_extent().x1) / 2
    assert abs((box.x0 + box.x1) / 2 - centre) < 1.0
