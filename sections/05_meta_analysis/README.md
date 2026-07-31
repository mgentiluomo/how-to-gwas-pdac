# Section 5 — Meta-analysis

**Lead:** Pelin Ünal

Ancestry-stratified association testing, harmonisation of summary statistics, and
fixed- and random-effects meta-analysis with heterogeneity.

## Contents

- `05_meta_analysis.qmd` — the website tutorial page.
- `scripts/` — four steps, PLINK 2 and base R.
- `results/` — outputs of one locked-in run on `pdac_demo`.

## Scope and honesty note

The three strata are ancestry groups within a single genotype file, not independently
recruited cohorts. The statistical machinery demonstrated is real; the harmonisation
problems of a genuine multi-cohort meta-analysis cannot arise here, and the page says so
explicitly. The harmonisation checks are nonetheless run in full, because a script that
assumes agreement will never detect disagreement.

## Manuscript code tags

| Tag | Script |
|---|---|
| CODE 5-01 | `scripts/01_stratified_gwas.sh` |
| CODE 5-02 | `scripts/02_harmonise.R` |
| CODE 5-03 | `scripts/03_meta_analysis.R` |
| CODE 5-04 | `scripts/04_meta_plots.R` |

## Key numbers from the demonstration run

| | |
|---|---|
| Strata | EUR 219/392, AFR 138/149, EAS 141/176 |
| Effective N, combined | 1,161.7 against 562.0 for Europeans alone |
| Effect-allele flips required | 37,738 AFR, 35,843 EAS, within one genotype file |
| Strand-ambiguous variants | 1,916 (0.6%), 354 in the unresolvable frequency band |
| Variants in all three strata | 311,478 of 405,000 |
| λ per stratum / meta | 1.015, 1.022, 0.997 / 0.994 |
| Causal variant P | 1.3 × 10⁻⁸ European alone → **2.9 × 10⁻¹⁸ meta** |
| Causal variant OR | 2.36 per risk allele; I² = 41.5%, P het = 0.18 |
| Simulated truth | OR 2.60; the meta estimate sits slightly below it, not above |
| TERT/CLPTM1L | 0.16 in Europeans → 0.021 across the strata: a real effect, still not recovered |
| Heterogeneity | the **causal** variant has the highest I² in the region while three of its proxies have none |

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- The meta-analysis arithmetic is written out rather than delegated to METAL or GWAMA,
  so that the reader can see what those tools compute.
- Open a PR to `main` when ready; a maintainer reviews and merges.
