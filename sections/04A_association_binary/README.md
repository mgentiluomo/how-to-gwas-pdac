# Section 4A — Association testing (binary outcome)

**Lead:** Elif Öz
**Status:** ☐ not started ☐ in draft ☑ ready for review ☐ harmonised

Single-variant case-control association testing on the European analysis set defined
in Section 2, with the diagnostics needed before any result is interpreted.

## Contents

- `04A_association_binary.qmd` — the website tutorial page.
- `scripts/` — one script per step, PLINK 2 and base R.
- `results/` — outputs of one locked-in run on `pdac_demo`, displayed on the page.

## Manuscript code tags

| Tag | Script |
|---|---|
| CODE 4A-01 | `scripts/01_make_covariates.R` |
| CODE 4A-02 | `scripts/02_association.sh` |
| CODE 4A-03 | `scripts/03_qq_lambda.R` |
| CODE 4A-04 | `scripts/04_manhattan.R` |
| CODE 4A-05 | `scripts/05_summary.R` |

## Key numbers from the demonstration run

| | |
|---|---|
| Analysis set | 224 cases, 387 controls |
| Variants tested | 396,396 autosomal; 392,331 with a valid P |
| Firth fallback invoked | 2,221 variants |
| Genomic inflation | λ = 1.028, λ₁₀₀₀ = 1.098 |
| P < 5 × 10⁻⁸ | 1 variant |
| P < 1 × 10⁻⁵ | 8 variants |
| Lead variant | 9:133249045:A:G (*ABO*), OR 3.90 (95% CI 2.47–6.15), P = 5.5 × 10⁻⁹ |
| Simulated truth at that variant | OR 2.40 — the recovered estimate is inflated by ~58%, the winner's curse |

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- Association testing is restricted to the autosomes: chrX in this dataset is simulated
  for the sex check only.
- Open a PR to `main` when ready; a maintainer reviews and merges.
