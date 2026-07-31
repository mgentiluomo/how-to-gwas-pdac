# Section 4A — Association testing (binary outcome)

**Lead:** Elif Öz

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
| Analysis set | 219 cases, 392 controls |
| Variants tested | 409,105 autosomal; 405,000 with a valid P |
| Firth fallback invoked | 2,460 variants |
| Genomic inflation | λ = 1.016, λ₁₀₀₀ = 1.058 |
| P < 5 × 10⁻⁸ | 11 variants, all at *ABO* |
| P < 1 × 10⁻⁵ | 13 variants |
| Lead variant | 9:133266804:G:T (*ABO*), **not the causal one**; OR 2.41 per the risk allele (95% CI 1.81 to 3.20), P = 1.4 × 10⁻⁹ |
| Simulated causal variant | 9:133273682:A:T, third by P; OR 2.34, P = 4.4 × 10⁻⁹ |
| Scored against the truth | generative OR 2.60, so the recovered estimate sits slightly **below** it. At 99% power there is no selection to inflate it |
| Allele orientation | the risk allele is the reference and the major allele, so PLINK prints 0.415 rather than 2.41 |

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- Association testing is restricted to the autosomes: chrX in this dataset is simulated
  for the sex check only.
- Open a PR to `main` when ready; a maintainer reviews and merges.
