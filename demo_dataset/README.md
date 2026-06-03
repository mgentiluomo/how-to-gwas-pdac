# PDAC GWAS Tutorial — Demo Dataset Preparation

Pipeline to build a small, fully reproducible, **didactic** GWAS dataset for the
manuscript:

> *"How to carry out a Genome-Wide Association Study: a step-by-step annotated
> guide using pancreatic cancer as a case study"*
> (TRANSPAN COST Action CA21116, Working Group 1).

The dataset is derived from **HAPNEST**, a publicly available, fully synthetic
collection of genotypes generated within the INTERVENE H2020 project. Because
HAPNEST contains **no real participants**, it can be redistributed, cited, and
will remain available for the long term — solving the three problems that make
real cohorts unsuitable for a teaching paper (data hosting, long-term
persistence, and patient anonymity), while still preserving realistic patterns
of **linkage disequilibrium (LD)** and **allele frequency**.

> Source: HAPNEST pre-generated data, EBI BioStudies accession **S-BSST936**.
> Citing reference for the manuscript: [PMID/DOI placeholder — HAPNEST,
> Wharrie et al., *Bioinformatics* 2023] <!-- highlight yellow -->

---

## Design choices

| Parameter        | Value                                          | Rationale                                                                 |
|------------------|------------------------------------------------|---------------------------------------------------------------------------|
| Source dataset   | HAPNEST pre-generated (EBI S-BSST936)          | Synthetic, public, citable; preserves real LD / MAF                       |
| Final sample     | **1,000 cases + 5,000 controls** (n = 6,000)   | Realistic scale of an early PDAC GWAS (e.g. PanScan I)                     |
| Ancestry mix     | EUR 80% / AFR 10% / EAS 5% / SAS 5%            | Mimics a European-led cohort; enables teaching of population stratification |
| Chromosomes used | 1, 5, 9, 13, 16, 17                            | Each carries one simulated causal locus (keeps the toy dataset small)     |
| Causal loci      | 6 **real, well-replicated** PDAC GWAS regions  | Students "rediscover" loci they can look up in the literature             |
| Effect sizes     | per-allele OR 1.15–1.25                         | Honest, realistic expectations about GWAS effect sizes                     |
| Heritability     | h² ≈ 0.15 (liability scale)                     | Consistent with PDAC family/twin estimates                                |
| Build            | GRCh38                                          | Matches HAPNEST                                                           |

### The six simulated causal loci

These are **genuine** pancreatic ductal adenocarcinoma (PDAC) susceptibility
regions reported by the PanScan / PanC4 consortia. The pipeline does **not**
hard-code real rsIDs (HAPNEST uses its own synthetic variant IDs); instead it
picks, within each region, the common variant nearest the target position and
designates it as causal.

| Chr | Region   | Nearest gene  | Approx. GRCh38 position | Reference (placeholder)        |
|-----|----------|---------------|-------------------------|--------------------------------|
| 1   | 1q32.1   | NR5A2         | chr1:200,000,000        | [PMID placeholder] <!-- y -->  |
| 5   | 5p15.33  | TERT / CLPTM1L| chr5:1,320,000          | [PMID placeholder] <!-- y -->  |
| 9   | 9q34.2   | ABO           | chr9:136,149,000        | [PMID placeholder] <!-- y -->  |
| 13  | 13q22.1  | KLF5          | chr13:73,330,000        | [PMID placeholder] <!-- y -->  |
| 16  | 16q23.1  | BCAR1         | chr16:75,180,000        | [PMID placeholder] <!-- y -->  |
| 17  | 17q25.1  | LINC00673     | chr17:73,261,000        | [PMID placeholder] <!-- y -->  |

> Verify the exact coordinates/genes against the current literature before
> final submission, and replace the placeholders with Vancouver-formatted
> references.

---

## A note on statistical power (read this before you panic)

With **n = 6,000** and **realistic** per-allele odds ratios (1.15–1.25), most
of the six causal loci will **not** reach the conventional genome-wide
significance threshold of *p* < 5 × 10⁻⁸. **This is intentional.** It is one of
the most important lessons of the tutorial: real GWAS need tens of thousands of
samples precisely because individual common-variant effects are small. The
strongest one or two loci (e.g. ABO) will usually pop out; the rest illustrate
why consortia and meta-analysis (Section 5) exist.

If, for a particular teaching session, you want a cleaner Manhattan plot with
all six loci visible, you can inflate the effect sizes via the `--or-scale`
option in `03_simulate_phenotype.R` — but please keep the realistic default for
the published example.

---

## Pipeline overview

```
01_download_hapnest.sh     Download the HAPNEST subset (selected chromosomes)
02_subset_individuals.R    Choose 6,000 individuals with the target ancestry mix
02b_apply_subset.sh        Extract those individuals with PLINK2 and merge chromosomes
03_simulate_phenotype.R    Simulate a PDAC-like phenotype from the 6 causal loci
04_run_gwas.sh             Reference run: PCA + logistic-regression GWAS (PLINK2)
05_plots.R                 Manhattan plot, QQ plot and genomic inflation (λ_GC)
```

Run them in order. Each script prints, at the end, the exact command for the
next step.

```bash
bash    scripts/01_download_hapnest.sh
Rscript scripts/02_subset_individuals.R
bash    scripts/02b_apply_subset.sh
Rscript scripts/03_simulate_phenotype.R
bash    scripts/04_run_gwas.sh
Rscript scripts/05_plots.R
```

---

## Software requirements

| Tool     | Version (tested) | Notes                                             |
|----------|------------------|---------------------------------------------------|
| PLINK    | **2.0** (a6+)    | `plink2` must be on your `PATH`                   |
| R        | ≥ 4.2            | packages: `data.table`, `optparse`, `ggplot2`     |
| bash     | ≥ 4              | GNU coreutils, `wget` or `curl`                   |

Install the R packages once:

```r
install.packages(c("data.table", "optparse", "ggplot2"))
```

---

## Directory layout

```
gwas_tutorial/
├── README.md
├── scripts/      # the six pipeline scripts
├── data/         # downloaded + intermediate genotype data (git-ignored)
├── results/      # GWAS summary statistics
└── docs/         # figures for the manuscript
```

> The `data/` folder can become large. Add it to `.gitignore`; only the
> scripts, the final summary statistics, and the figures need to live in the
> University of Pisa institutional repository.

---

## Reproducibility

Every random step uses a fixed seed (`--seed 2026`, configurable). Re-running
the pipeline on the same HAPNEST download reproduces an identical dataset,
phenotype, and set of results — a FAIR-data requirement discussed in Section 7.

---

*Maintained centrally by M. Gentiluomo & R. Farinella (University of Pisa).
Code style and final harmonisation are coordinated by the editors; please do
not change file names or the numbering scheme without prior agreement.*
