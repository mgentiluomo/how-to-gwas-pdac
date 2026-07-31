# Dataset construction scripts

## What is here

| Script | Purpose |
|---|---|
| `04_simulate_phenotype.R` | Generates the case/control phenotype, the covariates and the survival outcome from the released genotypes |

This is the script that produced the phenotype files in the current data
release. It is committed here so that the released data and the code that
generated it live in the same place and move together.

You do not need to run it to follow the tutorial. Download the released data
with `scripts/dev/download_demo_data.sh`.

```bash
Rscript 04_simulate_phenotype.R plink2 pdac_demo sample_ancestry.tsv .
```

Seed 2026, fixed in advance and never selected on the outcome. The script
re-estimates every generative parameter after simulating and stops with an error
rather than writing a dataset whose effects do not match their targets.

## What is not here

The genotype layer of the dataset, that is the filtering of the HGDP and
1000 Genomes callset to GSA positions, the anonymisation of identifiers, the
injected missingness and the simulated chromosome X, was produced by separate
scripts that are not currently in this repository. Those steps have not been
re-run since the first release and the genotypes are unchanged.

This is a gap in the provenance chain and it is recorded here rather than left
for a reader to discover. It should be closed before the resource is frozen and
deposited.

## What is deliberately withheld

The identifier mapping between the anonymised `MG####` labels and the source
sample identifiers is not distributed, here or in the release.

The causal-variant truth table **is** distributed, as `truth.tsv` in the release.
It is not a secret: the generative effect sizes are stated in the manuscript, in
the dataset README and in this script's header, and releasing the table is what
makes the exercise of scoring an analysis against a known truth reproducible.
