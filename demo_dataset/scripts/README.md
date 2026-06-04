# Dataset construction scripts

These scripts document how the demonstration dataset was built on real
HGDP+1KG reference genotypes. You do NOT need to run them to follow the
tutorial (just use `download_data.sh`); they are here for provenance and
reproducibility (FAIR, Section 7).

Run order:

1. `01_build_dataset.sh`        — filter HGDP+1KG VCFs to GSA positions,
                                  subset ancestries, make PLINK binaries
2. `02_degrade_and_anonymize.R` — inject realistic missingness; anonymise IDs
3. `03_make_chrX.R`             — add simulated chrX consistent with sex
4. `04_simulate_phenotype.R`    — simulate ABO-driven PDAC phenotype

All random steps use fixed seeds (2026/2027/2028/2029) for reproducibility.

NOTE: the private keys produced during construction (ID mapping, injected sex
discordances, causal-locus truth table) are intentionally NOT included here or
in the released data.
