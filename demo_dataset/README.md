# PDAC GWAS Tutorial — Demonstration Dataset

The running example for the manuscript

> *"How to carry out a Genome-Wide Association Study: a step-by-step annotated
> guide using pancreatic cancer as a case study"*
> (TRANSPAN COST Action CA21116, Working Group 1).

## ⚠️ Important: this is a teaching dataset

- **Genotypes are REAL** (public HGDP + 1000 Genomes reference data) but
  **sample identifiers are anonymised** (`MG0001`–`MG1461`) and bear no
  relation to any real individual's identity.
- **The phenotype is ENTIRELY SIMULATED.** No real pancreatic cancer cases are
  included. Case/control status, age, sex effects and survival outcomes were
  generated *in silico*.
- **The X chromosome is SIMULATED** for the sole purpose of demonstrating the
  sex-check QC step. It does not represent real chrX genotype data.
- This dataset must **not** be used for any biological or clinical inference
  about pancreatic cancer. Its only purpose is to teach GWAS methodology.

---

## Getting the data

The genotype/phenotype files are **not** stored in git (they are large
binaries). Download them from the GitHub Release with:

```bash
bash demo_dataset/download_data.sh
```

This fetches the files into `demo_dataset/data/` (git-ignored). The same files
are mirrored on Zenodo (DOI to be added on publication).

---

## Files

| File | Description |
|------|-------------|
| `pdac_demo.bed/.bim/.fam` | PLINK binary genotypes (430,000 variants, 1,461 individuals) |
| `phenotype.txt`     | `FID IID PHENO` (1 = control, 2 = case) |
| `covariates.txt`    | `FID IID SEX AGE` |
| `survival.txt`      | `FID IID TIME` (months) `EVENT` (1 = event, 0 = censored) |
| `sample_ancestry.tsv` | `IID`, genetic-ancestry group (`eur`/`afr`/`eas`) |

## Dataset composition

- **1,461 individuals**: 762 EUR, 349 AFR, 350 EAS — a deliberately balanced
  multi-ancestry design so that population stratification and PCA (Section 2)
  can be taught on a clearly structured dataset.
- **430,000 variants**: 424,000 autosomal + 6,000 simulated chrX. Variants are
  restricted to Illumina GSA-24v3 positions (GRCh38).
- Realistic genotype missingness was intentionally introduced (8,480
  low-call-rate variants; 25 low-call-rate samples) so that QC filters
  (`--geno`, `--mind`) have a visible effect.
- Variant IDs in `chr:pos:ref:alt` format (GRCh38).

## Simulated phenotype design

- **Single causal locus in the ABO region** (9q34.2), the most robust and
  most-replicated PDAC risk locus (Amundadottir et al. 2009; Rizzato et al.
  2013), with no polygenic background.
- The effect is deliberately larger than the true ABO effect on pancreatic
  cancer, which is near 1.2. Without that amplification no locus would be
  detectable at this sample size and the guide would have nothing to
  demonstrate downstream. The direction matches the literature.
- **Important caveat on the effect size.** The phenotype was generated under a
  liability-threshold model, while the analysis estimates an odds ratio on the
  log-odds scale. The two scales differ by a factor of roughly 1.3 to 1.8, so
  the odds ratio observed in the released data is not the parameter specified
  during construction. Estimated from the two ancestry strata that were not
  used for discovery, and therefore free of selection bias, the per-allele
  odds ratio in the released realisation is approximately 3.65 (95% CI 2.63 to
  5.09). Use that figure, not any generative value, when scoring an analysis
  of this dataset against a known truth. A replacement simulation on the
  log-odds scale, which validates itself by re-estimating the parameters it
  generated, is in preparation.
- Covariates: sex (slightly higher male risk) and age (cases older).
- The intended workflow analyses the **EUR subset only**. Before quality
  control the European group contains 254 cases and 508 controls; after it,
  the analysis set is **224 cases and 387 controls**, an effective sample
  size of 568.

## Expected results from the canonical run

| Quantity | Value |
|---|---|
| Samples after quality control | 1,217 |
| Variants after quality control | 414,695 |
| European analysis set | 224 cases, 387 controls |
| Genomic inflation factor | 1.028 (lambda1000 = 1.10) |
| Genome-wide significant loci | 1, at ABO |
| Lead variant | 9:133249045:A:G |
| P at the lead variant | 5.5e-9 |
| Odds ratio at the lead variant | 3.90 (95% CI 2.47 to 6.15) |

If your numbers differ, check that you are using the release tag the guide
links to rather than the most recent one.

## How the dataset was built (reproducibility)

The dataset was generated on real reference genotypes by the scripts in
[`scripts/`](scripts/), in this order:

| Script | Purpose |
|--------|---------|
| `01_build_dataset.sh`      | Filter HGDP+1KG VCFs to GSA positions, subset ancestries, make PLINK binaries |
| `02_degrade_and_anonymize.R` | Inject realistic missingness; anonymise IDs to MG#### |
| `03_make_chrX.R`           | Add a simulated chrX consistent with karyotypic sex (for the sex-check) |
| `04_simulate_phenotype.R`  | Simulate the ABO-driven PDAC phenotype + covariates + survival |

These scripts document provenance; you do **not** need to run them to follow
the tutorial — just download the released data. The private keys used during
construction (ID mapping, injected sex discordances, causal-locus truth table)
are **not** distributed.

## Provenance of genotypes

Real genotypes derived from the **HGDP + 1000 Genomes** joint callset
(gnomAD v3.1), phased reference panel:

- gnomAD v3.1 HGDP+1KG (Koenig et al., *Genome Research* 2024)
- HGDP (Bergström et al., *Science* 2020)
- 1000 Genomes Project (*Nature* 2015)
- Source phased panel: Zenodo record `10.5281/zenodo.18156285`

Only autosomal biallelic SNPs at GSA positions were retained; the chrX was
simulated separately (see above).

## License

CC-BY 4.0. Please cite the tutorial paper and the sources above.

---

*Maintained centrally by M. Gentiluomo & R. Farinella (University of Pisa).*
