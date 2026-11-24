# Manuscript figures

Base R, no packages beyond those of the toolchain. Each script reads only
results files of the canonical run (`sections/*/results/`) and of Section 8,
checks them where it can (the power curve against the tabulated values, the
I-squared fractions against the summary table), and writes to `FIG_OUT`
(default `../manuscript/figures/`, outside this repository):

- `<name>.pdf` (cairo, fonts embedded) and `<name>.png` (600 dpi), 174 mm wide;
- `<name>_sources.tsv`: every file the figure read, with its MD5.

```bash
for f in scripts/dev/figures/fig[0-9]*.R; do Rscript "$f"; done
```

| Script | Figure |
|---|---|
| `fig2_power.R` | Predicted against realised power over the replicates; smallest detectable OR against MAF |
| `fig3_calibration.R` | lambda of the logistic scan by number of components; permutation null of lambda in the Cox scan |
| `fig4_replicates.R` | I-squared at the causal variants over the replicates; ABF credible-set coverage and rank of the causal variant |
