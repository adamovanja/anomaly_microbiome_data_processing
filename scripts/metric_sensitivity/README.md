# Sensitivity of the post-antibiotic diversity deficit to the alpha diversity metric

Repeats the covariate-matched comparison of Fig. 4c of the manuscript for Faith PD and
Shannon diversity (manuscript revision, reviewer 1, point 4). Deficits relative to matched
unexposed infants are expressed in percent of the matched mean so that both metrics share one
axis; binning around the first to third exposure follows
Fig. 4c (see the module docstring).

- `metric_sensitivity.py`: analysis and figure; reads `data/final`.
- `test_metric_sensitivity.py`: binning, matching, window selection, test directions and outputs.
- Outputs written next to the script: `metric_sensitivity_windows.csv` (per sample),
  `metric_sensitivity_tests_vs_pre.csv`, `metric_sensitivity_tests_between_metrics.csv`,
  `metric_sensitivity_concordance.csv` and `supp_metric_sensitivity.pdf`.

Run from this directory (needs matplotlib, pandas, scipy and pytest, see `requirements.txt`):

```bash
pytest -q
python metric_sensitivity.py --out-png preview.png
```

The figure is copied into the manuscript repository as `figures/supp_metric_sensitivity.pdf`.
