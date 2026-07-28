# Section 6 — Fine mapping

**Lead:** Burçak Otlu, Erdi Kılıç
**Status:** ☐ not started ☐ in draft ☑ ready for review ☐ harmonised

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
| European set, credible set size | 77 variants |
| European set, top variant | `9:133249045:A:G`, posterior probability 0.51, ranked first |
| Three strata combined, credible set size | **1 variant** |
| Three strata combined, top variant | `9:133249045:A:G`, posterior probability 0.99 |
| Scored against the simulated truth | the causal variant is recovered in both, and uniquely identified only after combining ancestries |

The collapse from 77 variants to one is the section's central result. It is not a
power effect: linkage disequilibrium differs between populations, so variants that
cannot be distinguished within one ancestry are distinguished by disagreement between
ancestries. It is the mirror image of the heterogeneity pattern reported in Section 5,
where the causal variant showed I² = 0 while its correlated neighbours reached 79%.

## Conventions

- Build GRCh38; PLINK 2 and base R, no packages to install.
- In-sample LD throughout, since individual-level data are available.
- Open a PR to `main` when ready; a maintainer reviews and merges.
