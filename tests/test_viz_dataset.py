"""Tests for the cohort description figure (Fig. 3 of the manuscript)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src import viz_dataset as viz  # noqa: E402


def _cohort(rng, n_unexposed=6, n_exposed=5, drop_fd=False):
    """Samples of unexposed and exposed infants with the notebook's columns."""
    rows, events = [], []
    milk = ["bd", "mixed", "fd", None] if not drop_fd else ["bd", "mixed", None]
    for i in range(n_unexposed + n_exposed):
        host = f"h{i:02d}"
        exposed = i >= n_unexposed
        exposure = 8.0 + i if exposed else np.nan
        if exposed:
            events.append(dict(host_id=host, abx_start_age_months=exposure))
        ages = list(rng.choice(np.arange(40, 1100, 20), size=4, replace=False))
        if exposed:
            # one sample within the first month after the exposure
            ages.append(int(round((exposure + 1.0) * 30.437)))
        for age in sorted(set(ages)):
            months = age / 30.437
            since = months - exposure if exposed and months >= exposure else np.nan
            rows.append(
                dict(
                    host_id=host,
                    age_days=int(age),
                    age_months_rounded05=round(months * 2) / 2,
                    age_months_rounded1=float(round(months)),
                    diet_milk=milk[i % len(milk)],
                    diet_weaning=["no", "yes", None][i % 3],
                    abx_any_last_t_dmonths=since,
                    abx_any_cumcount=float(not np.isnan(since)),
                    abx_max_count_ever=float(exposed),
                )
            )
    return pd.DataFrame(rows), pd.DataFrame(events)


def test_prepare_samples_marks_samples_within_the_threshold_and_exposed_hosts():
    md_df, _ = _cohort(np.random.default_rng(0))
    prepared = viz.prepare_samples(md_df, th_months=3).set_index(
        ["host_id", "age_days"]
    )
    expected = (md_df["abx_any_last_t_dmonths"] <= 3).set_axis(
        pd.MultiIndex.from_frame(md_df[["host_id", "age_days"]])
    )
    assert prepared["sample_lt_xm_after_abx"].sort_index().equals(expected.sort_index())
    assert prepared["sample_lt_xm_after_abx"].sum() > 0
    per_host = prepared.groupby("host_id")["max_abx_w_microbiome"].first()
    assert (per_host.loc[["h00", "h06"]] == [0.0, 1.0]).all()


def test_split_by_exposure_separates_hosts_and_their_events():
    md_df, events = _cohort(np.random.default_rng(1))
    groups = viz.split_by_exposure(viz.prepare_samples(md_df, 3), events)
    assert list(groups) == ["w/o", "with"]
    assert groups["w/o"][0]["host_id"].nunique() == 6
    assert groups["with"][0]["host_id"].nunique() == 5
    assert groups["w/o"][1].empty
    assert set(groups["with"][1]["host_id"]) == set(groups["with"][0]["host_id"])


def test_diet_counts_by_age_keeps_category_order_and_fills_missing_with_zero():
    md_df, _ = _cohort(np.random.default_rng(2), drop_fd=True)
    counts = viz.diet_counts_by_age(md_df, "diet_milk")
    assert list(counts.columns) == ["bd", "mixed", "fd", "unknown"]
    assert (counts["fd"] == 0).all()
    assert counts["unknown"].sum() == md_df["diet_milk"].isna().sum()
    assert counts.values.sum() == len(md_df)
    weaning = viz.diet_counts_by_age(md_df, "diet_weaning")
    assert list(weaning.columns) == ["no", "yes", "unknown"]


def test_sample_characteristics_are_per_host_series():
    md_df, _ = _cohort(np.random.default_rng(3))
    chars = viz.sample_characteristics(md_df)
    assert list(chars) == [
        "# samples per host",
        "Duration between samples [days]",
        "Age at first sample [days]",
    ]
    n_hosts = md_df["host_id"].nunique()
    assert len(chars["# samples per host"]) == n_hosts
    assert len(chars["Duration between samples [days]"]) == len(md_df) - n_hosts
    assert (chars["Duration between samples [days]"] > 0).all()
    assert len(chars["Age at first sample [days]"]) == n_hosts


def _all_texts(ax):
    texts = [ax.title, ax.xaxis.label, ax.yaxis.label]
    texts += ax.get_xticklabels() + ax.get_yticklabels()
    if ax.get_legend() is not None:
        texts += ax.get_legend().get_texts()
    return [t for t in texts if t.get_text() != ""]


def test_build_dataset_figure_layout_fonts_and_legend():
    md_df, events = _cohort(np.random.default_rng(4))
    fig = viz.build_dataset_figure(md_df, events, th_months=3)
    np.testing.assert_allclose(fig.get_size_inches(), viz.FIGSIZE)
    assert len(fig.axes) == 7
    scatter, bars, violins = fig.axes[:2], fig.axes[2:4], fig.axes[4:]
    assert [ax.get_title() for ax in scatter] == [
        "Hosts w/o abx exposure (6)",
        "Hosts with abx exposure (5)",
    ]
    assert scatter[0].get_ylabel() == r"Host ID$_{\mathrm{no\ abx}}$"
    assert scatter[1].get_ylabel() == r"Host ID$_{\mathrm{abx}}$"
    assert scatter[0].get_legend() is None
    assert [t.get_text() for t in scatter[1].get_legend().get_texts()] == [
        "abx event",
        "sample <= 3m after abx",
        "sample > 3m after abx",
    ]
    assert [ax.get_legend().get_title().get_text() for ax in bars] == [
        "diet weaning",
        "diet milk",
    ]
    assert bars[1].get_xlabel() == "Age [months]"
    assert [ax.get_title() for ax in violins] == list(viz.CHARACTERISTIC_COLORS)
    for ax in violins:
        assert ax.get_xlim()[0] == 0
    letters = {t.get_text(): t for t in fig.texts if t.get_text() in "abc"}
    assert sorted(letters) == ["a", "b", "c"]
    for text in letters.values():
        assert text.get_fontweight() == "bold"
        assert text.get_fontsize() == viz.LETTER_FONT
    sizes = {t.get_fontsize() for ax in fig.axes for t in _all_texts(ax)}
    assert sizes == {viz.BASE_FONT}
    plt.close(fig)


def test_save_dataset_figure_writes_pdf(tmp_path):
    md_df, events = _cohort(np.random.default_rng(5))
    fig = viz.build_dataset_figure(md_df, events)
    pdf = viz.save_dataset_figure(fig, tmp_path)
    assert pdf == tmp_path / "dataset.pdf"
    assert pdf.read_bytes()[:4] == b"%PDF"
    plt.close(fig)
