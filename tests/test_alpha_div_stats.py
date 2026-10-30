"""Tests for the age-adjusted group differences in alpha diversity."""

import numpy as np
import pandas as pd
import pytest
import statsmodels.formula.api as smf

from src import alpha_div_stats as ads


def _cohort(
    rng,
    n_infants=150,
    levels=("A", "B"),
    effects=None,
    infant_sd=0.2,
    noise_sd=0.1,
    n_samples=8,
    max_age=24,
):
    """Infants with a constant group label and log-normal Faith PD
    around an age trend."""
    effects = effects or {}
    rows = []
    for i in range(n_infants):
        level = levels[i % len(levels)]
        offset = rng.normal(0, infant_sd)
        ages = rng.choice(np.arange(1, max_age + 1), size=n_samples, replace=False)
        for age in ages:
            effect = effects.get(level, 0.0)
            log_y = (
                np.log(5 + age) + np.log1p(effect) + offset + rng.normal(0, noise_sd)
            )
            rows.append(
                dict(
                    host_id=f"h{i}",
                    age_months_rounded1=float(age),
                    group=level,
                    div_alpha_faith_pd=float(np.exp(log_y)),
                    abx_any_cumcount=0.0,
                )
            )
    return pd.DataFrame(rows)


def test_flag_unexposed_marks_infants_without_exposure_before_any_sample():
    df = pd.DataFrame(
        {
            "host_id": ["A", "A", "B", "B", "C", "C"],
            "abx_any_cumcount": [0, 0, 0, 1, 0, 0],
            "abx_max_count_ever": [0, 0, 1, 1, 2, 2],
        }
    )
    flags = ads.flag_unexposed(df)
    # C was exposed only after its last sample and counts as unexposed
    assert flags.tolist() == [True, True, False, False, True, True]


def test_contrasts_recover_multiplicative_offset():
    df = _cohort(np.random.default_rng(1), effects={"B": 0.25})
    res = ads.age_adjusted_contrasts(df, "group", "A")
    assert res["level"].tolist() == ["B"] and res["reference"].tolist() == ["A"]
    row = res.iloc[0]
    assert row["pct_diff"] == pytest.approx(25.0, abs=5.0)
    assert row["ci_low"] < 25.0 < row["ci_high"]
    assert row["p_joint"] < 0.01
    assert row["n_level"] + row["n_reference"] == len(df)
    assert row["infants_level"] == 75 and row["infants_reference"] == 75


def test_contrasts_without_effect_cover_zero():
    df = _cohort(np.random.default_rng(2))
    row = ads.age_adjusted_contrasts(df, "group", "A").iloc[0]
    assert abs(row["pct_diff"]) < 5.0
    assert row["ci_low"] < 0.0 < row["ci_high"]
    assert row["p_joint"] > 0.05


def test_contrasts_drop_samples_with_missing_group():
    df = _cohort(np.random.default_rng(3))
    df.loc[df.index[:20], "group"] = np.nan
    row = ads.age_adjusted_contrasts(df, "group", "A").iloc[0]
    assert row["n_level"] + row["n_reference"] == len(df) - 20


def test_contrasts_with_three_levels_share_the_joint_p():
    df = _cohort(
        np.random.default_rng(4),
        levels=("Finland", "Estonia", "Russia"),
        effects={"Estonia": 0.1, "Russia": -0.05},
    )
    res = ads.age_adjusted_contrasts(df, "group", "Finland")
    assert res["level"].tolist() == ["Estonia", "Russia"]
    assert res["p_joint"].nunique() == 1
    assert res.set_index("level").loc["Estonia", "pct_diff"] == pytest.approx(10, abs=5)
    assert res.set_index("level").loc["Russia", "pct_diff"] == pytest.approx(-5, abs=5)


def test_confidence_intervals_account_for_between_infant_variation():
    df = _cohort(np.random.default_rng(5), infant_sd=0.4, noise_sd=0.05)
    row = ads.age_adjusted_contrasts(df, "group", "A").iloc[0]
    naive = smf.ols(
        "np.log(div_alpha_faith_pd) ~ C(age_months_rounded1) + C(group)", data=df
    ).fit()
    naive_width = 100 * (
        np.exp(naive.conf_int().loc["C(group)[T.B]", 1])
        - np.exp(naive.conf_int().loc["C(group)[T.B]", 0])
    )
    assert row["ci_high"] - row["ci_low"] > 1.5 * naive_width


def test_variance_components_are_recovered():
    df = _cohort(np.random.default_rng(6), n_infants=80, infant_sd=0.4, noise_sd=0.1)
    row = ads.age_adjusted_contrasts(df, "group", "A").iloc[0]
    assert row["sd_infant"] == pytest.approx(0.4, abs=0.1)
    assert row["sd_residual"] == pytest.approx(0.1, abs=0.03)


def test_summarise_panels_returns_one_row_per_contrast():
    rng = np.random.default_rng(8)
    panels = [
        ads.Panel(_cohort(rng), "group", "A", "Two levels"),
        ads.Panel(
            _cohort(rng, levels=("x", "y", "z")), "group", "x", "Three levels", "all"
        ),
    ]
    table = ads.summarise_panels(panels)
    assert table["panel"].tolist() == ["Two levels", "Three levels", "Three levels"]
    assert table["population"].tolist() == ["unexposed infants", "all", "all"]
    assert table["level"].tolist() == ["B", "y", "z"]
    assert table.columns.tolist() == [
        "panel",
        "population",
        "group",
        "level",
        "reference",
        "n_level",
        "n_reference",
        "infants_level",
        "infants_reference",
        "pct_diff",
        "ci_low",
        "ci_high",
        "p_joint",
        "sd_infant",
        "sd_residual",
    ]
    assert table.drop(columns=["panel", "population", "group"]).notna().all().all()


def test_formatting_of_percentages_and_p_values():
    assert ads.format_pct(4.2, -4.4, 13.6) == "+4% [−4, 14]"
    assert ads.format_pct(-10.1, -19.9, 1.0) == "−10% [−20, 1]"
    assert ads.format_p(0.0004) == "<0.001"
    assert ads.format_p(0.009) == "0.009"
    assert ads.format_p(0.07) == "0.07"
    assert ads.format_p(0.643) == "0.64"
