# Statistical power

**Lead:** Riccardo Farinella

What the analysis set can and cannot detect, computed before the association test.
Placed between population stratification and association testing because that is
where the manuscript places it: after the analysis set is fixed, before any result
exists to be interpreted.

## Contents

- `03B_statistical_power.qmd` — the website tutorial page.
- `scripts/` — two steps, base R only.
- `results/` — outputs of one locked-in run.

## Manuscript code tags

| Tag | Script |
|---|---|
| CODE 7-01 | `scripts/01_power_curves.R` |
| CODE 7-02 | `scripts/02_power_events.R` |

## Key numbers from the demonstration run

| | |
|---|---|
| Analysis set | 219 cases, 392 controls; 611 people behaving like 562 |
| All strata combined | effective N 1,161.7 |
| Smallest detectable OR at MAF 0.30 | 2.27 in the analysis set, 1.77 combined |
| Power at OR 1.3, MAF 0.30 | 0.0% |
| Power at OR 1.5, MAF 0.30 | 1.0% |
| *ABO*, generative OR 2.60 at MAF 0.39 | 99.1% power in the analysis set, 100% combined |
| TERT/CLPTM1L, generative OR 1.50 at MAF 0.49 | 2.0% power in the analysis set, 28.6% combined |
| Survival: cases with follow up | 219, the European analysis set |
| Survival: observed events | 167, the number that governs power |
| Smallest detectable HR at MAF 0.30 | 2.12 |
| Events needed for HR 1.3 at MAF 0.30 | 1,370 |

## Why this section carries weight

It shows that the shape of the result was fixed before the analysis ran: a typical
cancer susceptibility variant, odds ratio 1.3, had no chance of detection at this
size. The two causal variants make the point from both ends: *ABO* at 99.1% power
was always going to be found, TERT/CLPTM1L at 2.0% was not, and combining the
strata lifts the second only to 28.6%, which is still not enough.

It also settles what the winner's curse is and is not. At 99% power there is no
selection to inflate anything, and the recovered odds ratio, 2.34, sits slightly
below the generative 2.60. The curse is a consequence of low power rather than a
property of discovery, so the question to ask of a published effect size is what
power produced it.

## Conventions

- Base R only, no packages to install.
- Additive model, genome-wide threshold 5 × 10⁻⁸, target power 80%.
- Covariate adjustment deflates information by a factor of (1 − R²) and is discussed
  in the manuscript; the calculations here are unadjusted and therefore optimistic,
  which is stated rather than hidden.
