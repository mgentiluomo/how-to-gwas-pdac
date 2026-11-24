# Section 5: Meta-analysis

**Lead:** Pelin Ünal

Ancestry-stratified association testing, harmonisation of summary statistics, and
fixed- and random-effects meta-analysis with heterogeneity.

## Contents

- `05_meta_analysis.qmd`: the website tutorial page.
- `scripts/`: four steps, PLINK 2 and base R.
- `results/`: outputs of one locked-in run on `pdac_demo`.

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
| Strata | EUR 219/392, AFR 138/148, EAS 141/175 |
| Effective N, combined | 1,160.0 against 562.0 for Europeans alone |
| Effect-allele flips required | 37,766 AFR, 35,819 EAS, within one genotype file |
| Strand-ambiguous variants | 1,918 (0.6%), 354 in the unresolvable frequency band |
| Variants in all three strata | 311,843 of the 405,004 of the European stratum |
| λ per stratum (EUR, AFR, EAS) / meta | 1.007, 1.003, 1.025 / 0.993 |
| Causal variant P | 5.1 × 10⁻⁹ European alone → **1.2 × 10⁻¹⁸ meta** |
| Causal variant OR | 2.40 per risk allele; I² = 27.5%, P het = 0.25 |
| Simulated truth | OR 2.60; the meta estimate sits slightly below it, not above |
| TERT/CLPTM1L | P 0.20 in Europeans → 0.037 across the strata: a real effect, still not recovered |
| Heterogeneity | of the 13 variants at P < 5 × 10⁻⁸, eleven have I² = 0; the **causal** variant has 27.5% and one neighbour 31.7% |

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- The meta-analysis arithmetic is written out rather than delegated to METAL or GWAMA,
  so that the reader can see what those tools compute.
- Open a PR to `main` when ready; a maintainer reviews and merges.
