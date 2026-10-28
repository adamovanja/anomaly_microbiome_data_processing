"""Age-adjusted differences in alpha diversity between groups of infants.

For a grouping variable (e.g. delivery mode or country of residence), the
difference in Faith PD between each level and a reference level is estimated
with a linear mixed-effects model of log Faith PD on monthly age bins and the
grouping variable, with a random intercept per infant for the repeated samples
(fitted by restricted maximum likelihood). Coefficients are reported as percent
differences at a given age with 95% confidence intervals and a joint Wald test
of the grouping variable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

AGE_COL = "age_months_rounded1"
METRIC_COL = "div_alpha_faith_pd"
HOST_COL = "host_id"
ABX_COUNT_COL = "abx_any_cumcount"
LOG_COL = "_log_metric"
MINUS = "−"


@dataclass
class Panel:
    """Data and grouping of one panel of the grouped alpha diversity figure."""

    data: pd.DataFrame
    group: str
    reference: str
    title: str
    population: str = "unexposed infants"
    labels: dict = field(default_factory=dict)

    def label(self, level) -> str:
        return self.labels.get(level, str(level))

    @property
    def levels(self) -> list:
        """Reference level first, the other levels in sorted order."""
        others = sorted(set(self.data[self.group].dropna()) - {self.reference})
        return [self.reference] + others


def flag_unexposed(df: pd.DataFrame) -> pd.Series:
    """True for samples of infants without antibiotic exposure before any sample."""
    return df.groupby(HOST_COL)[ABX_COUNT_COL].transform("max") == 0


def _prepare(df: pd.DataFrame, group: str, reference: str, metric: str) -> pd.DataFrame:
    d = df.dropna(subset=[group, metric]).copy()
    d[group] = d[group].astype(str)
    if reference not in set(d[group]):
        raise ValueError(f"reference level {reference!r} not found in {group!r}")
    d[LOG_COL] = np.log(d[metric])
    return d


def _design(d: pd.DataFrame, age: str, group: str, reference: str) -> pd.DataFrame:
    """Design matrix: one indicator per age bin (no intercept) and one per
    non-reference level."""
    design = pd.get_dummies(d[age].astype(int), prefix="age", dtype=float)
    for lvl in sorted(d[group].unique()):
        if lvl != reference:
            design[f"level:{lvl}"] = (d[group] == lvl).astype(float)
    return design


def _fit(d: pd.DataFrame, design: pd.DataFrame):
    result = sm.MixedLM(d[LOG_COL], design, groups=d[HOST_COL]).fit(reml=True)
    if not result.converged:
        raise RuntimeError("the mixed model did not converge")
    return result


def _wald_p(result, names: list) -> float:
    restrictions = np.zeros((len(names), len(result.params)))
    for i, name in enumerate(names):
        restrictions[i, result.params.index.get_loc(name)] = 1.0
    return float(result.wald_test(restrictions, scalar=True).pvalue)


def _pct(value: float) -> float:
    return 100.0 * (np.exp(value) - 1.0)


def age_adjusted_contrasts(
    df: pd.DataFrame,
    group: str,
    reference: str,
    metric: str = METRIC_COL,
    age: str = AGE_COL,
) -> pd.DataFrame:
    """Percent difference of each level to the reference level within age bins,
    with the standard deviations of the random intercept and of the residuals."""
    d = _prepare(df, group, reference, metric)
    result = _fit(d, _design(d, age, group, reference))
    names = [n for n in result.params.index if n.startswith("level:")]
    p_joint = _wald_p(result, names)
    ci = result.conf_int()
    counts = d.groupby(group)[HOST_COL].agg(["size", "nunique"])
    rows = []
    for name in names:
        level = name[len("level:") :]
        rows.append(
            dict(
                level=level,
                reference=reference,
                n_level=int(counts.loc[level, "size"]),
                n_reference=int(counts.loc[reference, "size"]),
                infants_level=int(counts.loc[level, "nunique"]),
                infants_reference=int(counts.loc[reference, "nunique"]),
                pct_diff=_pct(result.params[name]),
                ci_low=_pct(ci.loc[name, 0]),
                ci_high=_pct(ci.loc[name, 1]),
                p_joint=p_joint,
                sd_infant=float(np.sqrt(result.cov_re.iloc[0, 0])),
                sd_residual=float(np.sqrt(result.scale)),
            )
        )
    return pd.DataFrame(rows)


def summarise_panels(
    panels: list, metric: str = METRIC_COL, age: str = AGE_COL
) -> pd.DataFrame:
    """One row per contrast of each panel."""
    tables = []
    for panel in panels:
        table = age_adjusted_contrasts(
            panel.data, panel.group, panel.reference, metric, age
        )
        table.insert(0, "group", panel.group)
        table.insert(0, "population", panel.population)
        table.insert(0, "panel", panel.title)
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def format_pct(estimate: float, low: float, high: float) -> str:
    """E.g. '+4% [−4, 14]' with a typographic minus sign."""
    text = f"{estimate:+.0f}% [{low:.0f}, {high:.0f}]"
    return text.replace("-", MINUS)


def format_p(p: float) -> str:
    if np.isnan(p):
        return "n/a"
    if p < 0.001:
        return "<0.001"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.2f}"


def write_summary(table: pd.DataFrame, path: Path) -> Path:
    path = Path(path)
    table.to_csv(path, index=False)
    return path
