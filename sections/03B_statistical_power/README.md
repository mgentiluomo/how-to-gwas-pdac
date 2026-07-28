# Statistical power

**Lead:** Riccardo Farinella
**Status:** ☐ not started ☐ in draft ☑ ready for review ☐ harmonised

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
| Analysis set | 224 cases, 387 controls; 611 people behaving like 568 |
| All strata combined | effective N 1,189.7 |
| Smallest detectable OR at MAF 0.30 | 2.24 in the analysis set, 1.76 combined |
| Power at OR 1.3, MAF 0.30 | 0.0% |
| Power at OR 1.5, MAF 0.30 | 1.1% |
| Simulated causal variant, OR 2.40 at MAF 0.085 | 9.8% power in the analysis set, 69.3% combined |
| Survival: individuals with follow up | 625 |
| Survival: observed events | 144, the number that governs power |
| Smallest detectable HR at MAF 0.30 | 2.25 |
| Events needed for HR 1.3 at MAF 0.30 | 1,370 |

## Why this section carries weight

It shows that the shape of the result was fixed before the analysis ran: a typical
cancer susceptibility variant, odds ratio 1.3, had no chance of detection at this
size. It also makes the mechanical link between low power and the winner's curse:
the causal variant had 9.8% power, was detected anyway, and its odds ratio came out
at 3.90. Detection in an underpowered scan is a selection
event, and the inflation follows from it.

## Conventions

- Base R only, no packages to install.
- Additive model, genome-wide threshold 5 × 10⁻⁸, target power 80%.
- Covariate adjustment deflates information by a factor of (1 − R²) and is discussed
  in the manuscript; the calculations here are unadjusted and therefore optimistic,
  which is stated rather than hidden.
