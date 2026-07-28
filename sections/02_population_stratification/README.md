# Section 2 — Population stratification

**Lead:** Elif Öz

Ancestry structure in the QC-passed data, definition of the analysis set, and the
principal components used as covariates in Section 4A.

## Contents

- `02_population_stratification.qmd` — the website tutorial page.
- `scripts/` — one script per step, PLINK 2 and base R.
- `results/` — outputs of one locked-in run on `pdac_demo`, displayed on the page.

## Manuscript code tags

| Tag | Script |
|---|---|
| CODE 2-01 | `scripts/01_ld_pruning.sh` |
| CODE 2-02 | `scripts/02_pca_all.sh` |
| CODE 2-03 | `scripts/03_pca_plots.R` |
| CODE 2-04 | `scripts/04_define_analysis_set.sh` |
| CODE 2-05 | `scripts/05_pca_within_eur.sh` |

## Key numbers from the demonstration run

| | |
|---|---|
| Variants after LD pruning | 222,538 of 414,695 |
| Variance explained, PC1 / PC2 / PC3 | 47.0% / 38.1% / 3.7% |
| QC-passed by group | 625 EUR, 319 EAS, 295 AFR |
| Analysis set | 611 European: 224 cases, 387 controls |
| Effective sample size | 567.5 |
| Variance explained within EUR, PC1 / PC2 | 23.9% / 11.8% |

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- Heavy comments: the reader has no bioinformatics background.
- References as `[DOI/PMID]` in square brackets, highlighted yellow (Vancouver at the end).
- Open a PR to `main` when ready; a maintainer reviews and merges.
