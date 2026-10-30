"""Tests for the grouped alpha diversity trajectory figure."""

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src import alpha_div_stats as ads  # noqa: E402
from src import viz_alpha_div as viz  # noqa: E402

METRIC, X = "div_alpha_faith_pd", "age_months_rounded1"


def _cohort(rng, levels=("A", "B"), n_infants=30):
    rows = []
    for i in range(n_infants):
        offset = rng.normal(0, 0.2)
        for age in rng.choice(np.arange(1, 25), size=6, replace=False):
            rows.append(
                dict(
                    host_id=f"h{i}",
                    age_months_rounded1=float(age),
                    group=levels[i % len(levels)],
                    div_alpha_faith_pd=float(
                        (5 + age) * np.exp(offset + rng.normal(0, 0.1))
                    ),
                    abx_any_cumcount=0.0,
                )
            )
    return pd.DataFrame(rows)


def test_group_trajectory_stats_omit_std_for_single_samples():
    df = pd.DataFrame(
        {X: [1.0, 1.0, 2.0], "group": ["A", "A", "A"], METRIC: [1.0, 3.0, 5.0]}
    )
    out = viz.group_trajectory_stats(df, METRIC, X, "group").set_index([X, "group"])
    assert out.loc[(1.0, "A"), "mean"] == 2.0
    assert out.loc[(1.0, "A"), "std"] == np.sqrt(2.0)
    assert out.loc[(1.0, "A"), "count"] == 2
    assert out.loc[(2.0, "A"), "count"] == 1
    assert np.isnan(out.loc[(2.0, "A"), "std"])


def test_level_colors_follow_the_husl_palette_in_sorted_level_order():
    three = viz.level_colors(["Russia", "Finland", "Estonia"], "Finland")
    assert three == {"Estonia": "#f77189", "Finland": "#50b131", "Russia": "#3ba3ec"}
    two = viz.level_colors(["yes", "no"], "no")
    assert two == {"no": "#f77189", "yes": "#36ada4"}


def test_build_group_figure_has_line_and_count_axes_only():
    rng = np.random.default_rng(0)
    panels = [
        ads.Panel(_cohort(rng), "group", "A", "Two levels"),
        ads.Panel(_cohort(rng, levels=("x", "y", "z")), "group", "x", "Three", "all"),
    ]
    fig = viz.build_group_figure(panels, METRIC, X)
    assert len(fig.axes) == 2 * len(panels)
    line_axes = fig.axes[0::2]
    for ax, panel in zip(line_axes, panels):
        legend = ax.get_legend()
        assert legend is not None
        assert len(legend.get_texts()) == panel.data["group"].nunique()
    matplotlib.pyplot.close(fig)


def test_lineplot_alpha_by_groups_writes_pdf_and_png(tmp_path):
    rng = np.random.default_rng(1)
    panels = [ads.Panel(_cohort(rng), "group", "A", "Two levels")]
    pdf = viz.lineplot_alpha_by_groups(panels, METRIC, X, tmp_path)
    assert pdf == tmp_path / "alpha_all_lineplot_splitcov.pdf"
    assert pdf.stat().st_size > 0
    assert (tmp_path / "alpha_all_lineplot_splitcov.png").stat().st_size > 0


def test_panels_by_p_orders_panels_by_increasing_joint_p():
    rng = np.random.default_rng(2)
    panels = [ads.Panel(_cohort(rng), "group", "A", title) for title in ["a", "b", "c"]]
    stats = pd.DataFrame(
        {"panel": ["a", "b", "b", "c"], "p_joint": [0.5, 0.01, 0.01, 0.2]}
    )
    assert [p.title for p in viz.panels_by_p(panels, stats)] == ["b", "c", "a"]


def _stats_two_panels():
    return pd.DataFrame(
        {
            "panel": ["Two levels", "Three levels", "Three levels"],
            "level": ["B", "y", "z"],
            "reference": ["A", "x", "x"],
            "pct_diff": [1.0, 2.0, 3.0],
            "ci_low": [-1.0, 0.0, 1.0],
            "ci_high": [3.0, 4.0, 5.0],
            "p_joint": [0.5, 0.01, 0.01],
        }
    )


def test_latex_pct_and_p_use_math_mode():
    assert viz.latex_pct(-9.4, -12.9, -4.6) == r"$-9$~\% [$-13$, $-5$]"
    assert viz.latex_pct(3.2, -7.9, 15.7) == r"$+3$~\% [$-8$, $16$]"
    assert viz.latex_p(0.0004) == r"$<0.001$"
    assert viz.latex_p(0.005) == "0.005"
    assert viz.latex_p(0.94) == "0.94"


def test_latex_table_orders_variables_by_p_without_panel_letters():
    rng = np.random.default_rng(3)
    panels = [
        ads.Panel(_cohort(rng), "group", "A", "Two levels"),
        ads.Panel(_cohort(rng, levels=("x", "y", "z")), "group", "x", "Three levels"),
    ]
    table = viz.latex_table(panels, _stats_two_panels())
    rows = [line for line in table.splitlines() if line.endswith(r"\\")]
    assert table.startswith(r"\begin{tabular}") and r"\bottomrule" in table
    assert rows[0] == r"Variable & Contrast & Difference [95\% CI] & $p$ \\"
    assert rows[1] == r"Three levels & y vs x & $+2$~\% [$0$, $4$] & 0.01 \\"
    assert rows[2] == r" & z vs x & $+3$~\% [$1$, $5$] & \\"
    assert rows[3] == r"Two levels & B vs A & $+1$~\% [$-1$, $3$] & 0.50 \\"
    assert len(rows) == 4


def test_write_latex_table_writes_the_file(tmp_path):
    rng = np.random.default_rng(4)
    panels = [
        ads.Panel(_cohort(rng), "group", "A", "Two levels"),
        ads.Panel(_cohort(rng, levels=("x", "y", "z")), "group", "x", "Three levels"),
    ]
    path = viz.write_latex_table(panels, _stats_two_panels(), tmp_path / "t.tex")
    assert path == tmp_path / "t.tex"
    assert path.read_text() == viz.latex_table(panels, _stats_two_panels())
