# Section 4B: Association testing, survival outcome

**Lead:** Pelin Ünal


## What goes in this folder
- `*.md`: the narrative text for this section (Markdown).
- annotated `*.R` / `*.sh` scripts demonstrating the step on the demo dataset
  (`../../demo_dataset/`).
- any small inputs needed to reproduce the example (no large genotype files).

## Conventions
- Build GRCh38; PLINK 2 + R; heavy comments for a non-bioinformatics audience.
- References as `[DOI/PMID]` in square brackets, highlighted yellow (Vancouver at the end).
- Open a PR to `main` when ready; a maintainer reviews and merges.

## Outside the canonical run
- `06_random_covariate_null.R`: the permutation null of `04_permutation_null.R`
  repeated with ten random N(0, 1) covariates in place of the ten within-set
  components (same permutations, same number of parameters), and the count of
  tied survival times. `scripts/dev/run_all.sh` does not call it; it reads the
  outputs of the canonical run in `results/`. `scripts/dev/run_experiments.sh`
  runs it with Section 8 and records it in `EXPERIMENTS_MANIFEST.txt`.
  Output: `results/pdac_demo_04B_randcov_*` and `pdac_demo_04B_ties.tsv`.
