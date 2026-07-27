# Section 5 — Meta-analysis

**Lead:** Pelin Ünal
**Status:** ☐ not started ☐ in draft ☑ ready for review ☐ harmonised

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
| Strata | EUR 229/396, AFR 142/153, EAS 141/178 |
| Effective N, combined | 1,189.7 against 580.4 for Europeans alone |
| Effect-allele flips required | 32,829 AFR, 30,878 EAS — within one genotype file |
| Strand-ambiguous variants | 1,833 (0.6%), 339 in the unresolvable frequency band |
| Variants in all three strata | 296,678 of 392,331 |
| λ per stratum / meta | 1.022, 0.991, 1.013 / 1.003 |
| Lead variant P | 5.4 × 10⁻⁹ European alone → **2.4 × 10⁻²² meta** |
| Lead variant OR | 3.74 (2.87–4.87); I² = 0%, P het = 0.44 |
| Simulated truth | OR 2.40 — meta-analysis does not cure the winner's curse |
| Heterogeneity at LD proxies | up to I² = 79%; one variant drops from P 1.3 × 10⁻⁹ to 0.0035 under random effects |

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- The meta-analysis arithmetic is written out rather than delegated to METAL or GWAMA,
  so that the reader can see what those tools compute.
- Open a PR to `main` when ready; a maintainer reviews and merges.
