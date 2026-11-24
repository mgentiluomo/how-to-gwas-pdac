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

## What is elsewhere, and what is not here

The first step of the genotype layer is in
[`scripts/data_release/build_gwas_tutorial_dataset.sh`](../../scripts/data_release/build_gwas_tutorial_dataset.sh):
the selection of the individuals from the HGDP and 1000 Genomes panels, their
ancestry labels and the filtering of the callset to GSA positions. It is the
script as it was run when the dataset was built, with the server paths turned
into environment variables; it has not been run again. Its header names the
public input files, with versions and addresses, and states what was checked
afterwards and what was not recorded.

The ancestry label in `sample_ancestry.tsv` is the genetic-ancestry group
that gnomAD v3.1 inferred for the sample (`population_inference.pop`, with
`nfe` and `fin` written as `eur`), not the population the sample was
collected under. The African and East Asian individuals are a random draw of
350 from those eligible, made with an awk shuffle whose result is not
guaranteed to be the same with another awk. One of the 350 African samples
drawn is not in the phased panel, so 349 were released.

The later steps, that is the pseudonymisation of identifiers, the injected
missingness and the simulated chromosome X, were produced by separate scripts
that are not in this repository. Those steps have not been re-run since the
first release and the genotypes are unchanged.

This is a gap in the provenance chain and it is recorded here rather than left
for a reader to discover. It should be closed before the resource is frozen and
deposited.

## What is not published, and what that does and does not mean

The identifier mapping between the pseudonymised `MG####` labels and the
source sample identifiers is not distributed, here or in the release. The
`MG####` identifiers are pseudonyms. Their correspondence with the source
identifiers can be reconstructed from public files (the order of the gnomAD
metadata file and the script in `scripts/data_release/`), and in any case
from the genotypes, which are public. It reveals nothing about a real person
beyond what the public panels already state, because case status, age,
survival and the X chromosome are simulated. The mapping is not published,
but withholding it is not a protective measure.

The causal-variant truth table **is** distributed, as `truth.tsv` in the release.
It is not a secret: the generative effect sizes are stated in the manuscript, in
the dataset README and in this script's header, and releasing the table is what
makes the exercise of scoring an analysis against a known truth reproducible.
