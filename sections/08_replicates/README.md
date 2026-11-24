# Section 8: Replicates of the simulated phenotype

The released dataset is one draw of the simulated phenotype. One draw shows
what can happen; it cannot show how often. This section draws new
case-control phenotypes on the same genotypes, with the same generator as the
release and a third causal locus, and analyses each one as the canonical run
does. The released dataset does not change, and this section is not part of
the canonical run: it reads the canonical run's outputs in `results/`.

## Design

- **Generator:** the logic of `demo_dataset/scripts/04_simulate_phenotype.R`.
  Step 01 checks that, with seed 2026 and the two release loci, it reproduces
  the released `phenotype.txt` and ages exactly, and stops otherwise.
- **Causal loci:**

  | Locus | Variant | Generative OR |
  |---|---|---|
  | ABO | 9:133273682:A:T | 2.60 |
  | TERT/CLPTM1L | 5:1286401:C:A | 1.50 |
  | 13q22.1 | 13:73362027:C:T | 1.97 |

  The 13q22.1 variant was chosen by a fixed rule: the array variant nearest
  to rs9543325 that is non-palindromic, has MAF 0.25 to 0.5 in the European
  set and at least 0.10 in the other strata. Its OR gives a predicted power
  of 50.7% in the European analysis set with the formula of Section 3B
  (N_eff 562 of the canonical set). Step 04 shows that this formula is
  optimistic under the generator: see `pdac_demo_08_expected_se.txt`.
- **Seeds:** replicate *i* uses seed 20261124 + *i*.
- **Number of replicates:** `GWAS_08_N_REP`, default 200.
- **Per replicate:** the European scan of Section 4A and the three-stratum
  meta-analysis of Section 5, on the LD-pruned European set (for lambda) plus
  every variant within 350 kb of each causal variant; ABF fine-mapping in
  those windows.
- **Components:** the European scan adjusts for SEX, AGE and PC1 to PCk,
  k from `pdac_demo_04A_npc_chosen.txt` (2 in the final canonical run); each
  stratum of the meta-analysis uses the components of its own `covar.txt`.
  Results committed in 01e7f98 used all 20 components in the European scan
  by mistake; corrected on 24 November 2026.
- **Held fixed:** samples, QC, strata, components and the number of
  components are those of the canonical run. Relatedness pruning and the
  Hardy-Weinberg test in controls depend on the phenotype and are not
  repeated per replicate.

## Run

`bash scripts/dev/run_experiments.sh` runs every step below and the 4B step
with random covariates, after the canonical run, and writes
`EXPERIMENTS_MANIFEST.txt` (commit, canonical run read, seeds, software
versions, MD5 of every result file). The steps one by one:

```bash
Rscript sections/08_replicates/scripts/01_simulate_replicates.R
Rscript sections/08_replicates/scripts/02_analyse_replicates.R
Rscript sections/08_replicates/scripts/03_summarise.R
Rscript sections/08_replicates/scripts/04_expected_information.R
Rscript sections/08_replicates/scripts/05_test_statistics.R
Rscript sections/08_replicates/scripts/06_components.R          # about 1,000 genome-wide scans
Rscript sections/08_replicates/scripts/06b_random_covariates.R
Rscript sections/08_replicates/scripts/07_causal_variants_and_structure.R
Rscript sections/08_replicates/scripts/08_expected_gain.R
Rscript sections/08_replicates/scripts/09_set_composition.R     # needs the unpublished identifier mapping
Rscript sections/08_replicates/scripts/10_simulated_detectable_or.R
Rscript sections/08_replicates/scripts/11_plink2_variant_set.R
```

Steps 05 to 09 were added on 24 November 2026 in answer to a second internal
review of the manuscript:

- **05, test statistics.** Wald (what PLINK 2 reports), score,
  likelihood-ratio and Firth tests at the three causal variants in every
  replicate, European scan and meta-analysis, and the calibration of the
  first three on 2,000 null variants. It decides whether realised power below
  the Section 3B formula belongs to the formula or to the Wald test.
- **06, components.** The European scan of every replicate with 0, 2, 5, 10
  and 20 components: lambda, power, and the number the rule of 04A step 01b
  would choose in each replicate.
- **06b, random covariates.** The scans with 5 and 20 components against
  scans with 5 and 20 covariates drawn from N(0, 1), on the pruned set: is
  the rise of lambda with the components the number of parameters alone?
- **07, causal variants and structure.** Whether the causal genotypes follow
  the within-European components, and whether the components explain the
  correlation between the ABO and TERT/CLPTM1L genotypes.
- **08, expected gain.** For the comparison with relatives retained (4C) and
  for the two relatedness rules (4A sensitivity): the change in the ABO
  chi-square expected from the change in effective sample size, beside the
  observed one. Reads committed tables only.
- **09, set composition.** Counts of the source populations in each analysis
  set. It needs the identifier mapping, which is published neither in the repository
  nor in the release; without it the step writes nothing. The committed table
  holds counts per population only.

- **10, simulated detectable OR.** The OR at 80% power by simulation at one
  frequency, with the risk on the minor or on the major allele, against the
  value of the Section 3B formula.
- **11, PLINK 2 and the variants loaded together.** One documented case in
  which the standard error of a variant changes with the other variants in
  the command.

How the ancestry labels were made: `demo_data/sample_ancestry.tsv` holds the
genetic-ancestry group gnomAD v3.1 inferred for each sample (column
`population_inference.pop` of the HGDP+1KG metadata), with `nfe` and `fin`
merged into `eur`. Step 09 checks this on every released individual. The
label is inferred from the genotypes and is not the region a population was
sampled in, so two individuals of the African stratum come from populations
of other regions (table `pdac_demo_08_label_rule.tsv`). The script that
wrote the file when the dataset was built is recorded in
`scripts/data_release/build_gwas_tutorial_dataset.sh`; it was not run again.

PLINK 2 and exact agreement: step 02 scans an extracted subset of variants.
The same model on all variants (step 06), or refitted in R (step 05), agrees
with it to about four significant digits and not to the last one, because
the result PLINK 2 returns for a variant depends slightly on the variants
loaded beside it (our reading: an iterative fit in single precision started
from the previous variant; not verified in its source). The headers of steps
05 and 06 give the sizes.

## Outputs

Written to `results/`. The replicate phenotypes in `results/sim/` are
regenerated by step 01 and are not committed.

| File | Content |
|---|---|
| `pdac_demo_08_power.tsv` | Realised against predicted power |
| `pdac_demo_08_winners_curse.tsv` | Estimated / generative OR, overall and given P < 5e-8; the same on the log OR scale with bootstrap intervals |
| `pdac_demo_08_finemap.tsv` | Credible-set coverage and rank of the causal variant |
| `pdac_demo_08_i2.tsv` | I-squared at the causal variants, where the true effect is identical by construction |
| `pdac_demo_08_lambda_summary.tsv` | Distribution of lambda |
| `pdac_demo_08_per_replicate.tsv` | Results for every replicate |
| `pdac_demo_08_lambda.tsv` | Lambda for every replicate |
| `pdac_demo_08_realised.tsv` | Realised ORs for every replicate |
| `pdac_demo_08_summary.txt` | Plain-text report |
| `pdac_demo_08_expected_se.tsv`, `.txt` | Standard error and power of the causal log OR: Section 3B formula, expected information of the joint generative model and of the single-locus model the scan fits (also with the effects of sex and age, or of everything but the locus, removed from the generator), realised |
| `pdac_demo_08_genotype_r.tsv` | Correlation of the risk-allele counts of the causal variants within each stratum |
| `pdac_demo_08_tests_power.tsv`, `_tests_per_replicate.tsv`, `_tests_null.tsv`, `_tests_summary.txt` | Step 05: realised power of the Wald, score, likelihood-ratio and Firth tests; null calibration |
| `pdac_demo_08_npc_summary.tsv`, `_npc_per_replicate.tsv`, `_npc_summary.txt` | Step 06: lambda and power by number of components |
| `pdac_demo_08_npc_rule.tsv`, `_npc_rule_summary.tsv`, `_npc_rule_counts.tsv` | Step 06: the number of components the rule of 04A step 01b chooses in each replicate |
| `pdac_demo_08_randcov_summary.tsv`, `_randcov_per_replicate.tsv`, `_randcov_summary.txt` | Step 06b: lambda with k components against k random covariates |
| `pdac_demo_08_causal_structure.tsv`, `_causal_structure_by_pc.tsv`, `_causal_partial_r.tsv`, `_causal_structure.txt` | Step 07: causal genotypes against the within-European components; partial correlations |
| `pdac_demo_08_expected_gain.tsv`, `.txt` | Step 08: observed against expected gain at the ABO variant |
| `pdac_demo_08_set_composition.tsv`, `_by_project.tsv`, `.txt` | Step 09: source populations and panels of the analysis sets, counts |
| `pdac_demo_08_label_rule.tsv` | Step 09: released label against the gnomAD inferred group and the region of the population, counts |
| `pdac_demo_08_simulated_detectable_or.tsv`, `_simulated_power_grid.tsv`, `.txt` | Step 10: OR at 80% power by simulation; power by OR |
| `pdac_demo_08_plink2_variant_set.tsv`, `.txt` | Step 11: one variant tested alone, within a window, and in R |
| `pdac_demo_08_tert_adjusted.tsv` | European TERT/CLPTM1L estimate, unadjusted and adjusted for ABO (and 13q22.1), replicates and released phenotype |
