# Section 6 — Fine mapping

**Lead:** Burçak Otlu, Erdi Kılıç

From an association peak to a credible set of candidate causal variants, with the
resolution scored against a known truth.

## Contents

- `06_finemapping_annotation.qmd` — the website tutorial page.
- `scripts/` — three steps, PLINK 2 and base R.
- `results/` — outputs of one locked-in run.

## Method note

Fine mapping is implemented as Wakefield's approximate Bayes factor, written out in
full in `02_finemap_abf.R`, which assumes a single causal variant per region. This is
a deliberate choice for a teaching page: it is short enough to read and understand,
and at the sample size of the demonstration a second signal would not be detectable.
The page states the assumption explicitly and directs the reader to SuSiE, FINEMAP or
CAVIAR whenever multiplicity is plausible, citing *TERT-CLPTM1L* as the PDAC locus
where it demonstrably is.

## Manuscript code tags

| Tag | Script |
|---|---|
| CODE 6-01 | `scripts/01_region_and_ld.sh` |
| CODE 6-02 | `scripts/02_finemap_abf.R` |
| CODE 6-03 | `scripts/03_finemap_plots.R` |

## Key numbers from the demonstration run

| | |
|---|---|
| Region | chr9:132.9–133.6 Mb, 204 variants |
| Prior variance on the log odds ratio | W = 0.04; coverage 95% |
| European set, credible set size | 11 variants spanning 12.2 kb |
| European set, top variant | `9:133266804:G:T`, posterior probability 0.24, **not the causal one** |
| European set, causal variant | ranked **third**, posterior probability 0.13 |
| Three strata combined, credible set size | 9 variants |
| Three strata combined, top variant | `9:133273682:A:T`, posterior probability 0.43, **the causal one** |
| Scored against the simulated truth | present in both sets; promoted from third to first by combining ancestries |

The section's central result is the change in **ranking**, not in size: the set
shrinks only from eleven variants to nine, while the causal variant moves from third
place to first and its posterior nearly doubles. Linkage disequilibrium differs
between populations, so variants that cannot be distinguished within one ancestry are
distinguished by disagreement between ancestries. More samples of the same ancestry
would have sharpened the P value without necessarily changing the order.

The mirror-image argument, that the causal variant should show no heterogeneity while
its proxies do, is **not** visible in this dataset. Section 5 reports the opposite,
because I² is barely estimable with fewer than 320 individuals per stratum. The
resolution gain above is real and is visible in the posteriors; the heterogeneity
signature that supposedly explains it is not.

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- In-sample LD throughout, since individual-level data are available.
- Open a PR to `main` when ready; a maintainer reviews and merges.
