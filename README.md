# How to carry out a genome-wide association study of a rare and complex trait

**A decision-focused, reproducible guide using pancreatic ductal adenocarcinoma (PDAC) as a worked example**

[![Companion website](https://img.shields.io/badge/companion%20site-live-blue)](https://mgentiluomo.github.io/how-to-gwas-pdac/)
[![Demo dataset](https://img.shields.io/badge/demo%20dataset-v0.1--data-green)](https://github.com/mgentiluomo/how-to-gwas-pdac/releases/tag/v0.1-data)
[![License: MIT](https://img.shields.io/badge/code-MIT-lightgrey)](LICENSE)

This repository is the companion resource for the methodology manuscript developed within the **TRANSPAN COST Action (CA21116), Working Group 1**; target journal to be defined. It pairs a narrative on how to design and run a GWAS of a rare, genetically complex trait such as pancreatic cancer with an annotated PLINK 2 and R pipeline, demonstrated on a purpose-built demonstration dataset.

The manuscript explains why each decision is made and when an analysis is or is not legitimate given the data in hand. This repository shows how each step is actually run. Neither replaces the other; they are intended to be used together.

**Companion website (start here to run the pipeline):** [mgentiluomo.github.io/how-to-gwas-pdac](https://mgentiluomo.github.io/how-to-gwas-pdac/)

---

## Status

Nine stages are published as executable worked analyses, each with annotated
scripts and committed results from a single canonical run. The remaining stages
of the workflow are covered in the manuscript narrative but are not yet worked
on the demonstration data, and are listed here as not built rather than as
forthcoming, so that what is available is never overstated.

| Stage | Status |
|---|---|
| Setup and first PLINK command | Published, with automated project and tool setup |
| Study design and genotyping technologies | Published, narrative |
| Genotyping quality control | Published: nine steps, annotated scripts, committed results |
| Ancestry and population structure | Published, with scripts and results |
| Statistical power | Published, with scripts and results |
| Association testing, binary outcome | Published, with scripts and results |
| Meta-analysis across ancestry strata | Published, with scripts and results |
| Fine mapping | Published, with scripts and results |
| Linkage disequilibrium and what it decides | Published, with scripts and results |
| Imputation | Not built |
| Association testing, time-to-event outcome | Not built |
| Association software comparison | Not built |
| Functional annotation | Not built |
| Reporting, replication and FAIR data | Not built |

All figures on the site and in the manuscript come from one canonical run of
the published pipeline. Where a number appears in more than one place it is the
same number; if you find one that is not, it is a defect and an issue is
welcome.

## The demonstration dataset

`pdac_demo` uses real, anonymised genotypes from the public HGDP and 1000 Genomes reference panel, restricted to Illumina GSA-24v3 positions on GRCh38, with a fully simulated PDAC phenotype and deliberately injected quality control artefacts. No controlled-access application is required.

| | |
|---|---|
| Individuals | 1,461, 762 EUR, 349 AFR, 350 EAS |
| Variants | 430,000, 424,000 autosomal, 6,000 simulated chrX |
| Injected artefacts | 8,480 low-call-rate variants, 25 low-call-rate samples, 4 sex discordances |
| Causal locus | *ABO* (9q34.2), single locus, no polygenic background |
| Pre-QC EUR case:control | 254 cases, 508 controls, approximately 1:2 |
| After QC | 1,217 individuals retained (83.3%) and 414,695 variants retained (96.4%) |

Genotypes and phenotypes are not stored in git, see `.gitignore`. Download them with `bash scripts/dev/download_demo_data.sh`, or from the [release page](https://github.com/mgentiluomo/how-to-gwas-pdac/releases/tag/v0.1-data). Full provenance, construction scripts and the didactic rationale for the dataset design are in [`demo_dataset/README.md`](demo_dataset/README.md).

This dataset teaches GWAS under realistic, difficult conditions on purpose, it is not an idealised example in which every locus succeeds.

## Repository structure

```
how-to-gwas-pdac/
├── README.md                you are here
├── CONTRIBUTING.md          workflow for contributors, branches, PRs, review
├── LICENSE                  MIT for code, manuscript text CC-BY-4.0
├── index.qmd                companion website home page
├── getting_started.qmd      setup, project structure, tool install, first command
├── _quarto.yml              website navigation and build configuration
├── styles.css
├── demo_dataset/            the pdac_demo dataset: README, download script, construction scripts
├── sections/                one folder per manuscript section
│   └── 01B_genotyping_qc/   the QC walkthrough, its scripts and its committed results
├── scripts/dev/             project bootstrap and tool installation helpers
├── docs/                    figures, the pipeline diagram, setup helper documents
└── env/                     software versions and R package list, for reproducibility
```

## Section leads

| Section | Manuscript topic | Lead | Folder |
|---|---|---|---|
| 0 | Introduction and overview | Daniele Campa | `sections/00_introduction/` |
| 1A | Study design and consortium | Bogdan Ungureanu, Anca Riza Costache | `sections/01A_study_design/` |
| 1B | Genotyping quality control | Murat Güler | `sections/01B_genotyping_qc/` |
| 2 | Population stratification | Elif Öz | `sections/02_population_stratification/` |
| 3 | Imputation | Pelin Ünal | `sections/03_imputation/` |
| 4A | Association testing, binary outcome | Elif Öz | `sections/04A_association_binary/` |
| 4B | Association testing, survival outcome | Pelin Ünal | `sections/04B_association_survival/` |
| 4C | Association testing, pipeline and software | Murat Güler | `sections/04C_association_other/` |
| 5 | Meta-analysis | Pelin Ünal | `sections/05_meta_analysis/` |
| 6 | Fine mapping | Burçak Otlu, Erdi Kılıç | `sections/06_finemapping_annotation/` |
| 6.5 | Functional annotation | Chiara Corradi | `sections/06.5_gwas_annotations/` |
| 7 | Reporting, replication, FAIR data | Manuel Gentiluomo | `sections/07_reporting_fair/` |
| — | Shared helper scripts, not a manuscript section | Murat Güler | `sections/08_helper_scripts/` |

Code production: Riccardo Farinella and Murat Güler. Final code review and standardisation: Manuel Gentiluomo.

## Quick start

```bash
git clone https://github.com/mgentiluomo/how-to-gwas-pdac.git
cd how-to-gwas-pdac
bash scripts/dev/download_demo_data.sh        # fetches pdac_demo into demo_data/
quarto preview                            # live preview of the companion website
```

Readers following the website do not need to clone anything: `getting_started.qmd` sets up the project folders and downloads the data and the tools automatically.

## Conventions

* Tools: PLINK 1.9 and PLINK 2.0, plus R. One step, one script, heavily commented for a non-bioinformatics audience.
* Build: GRCh38 throughout, variant IDs in `chr:pos:ref:alt` format.
* Reproducibility: all random steps use a fixed seed, default 2026.
* References: provisional placeholders as `[DOI/PMID]`, highlighted yellow in the manuscript, final formatting is Vancouver style.
* Data: raw or intermediate genotype files, `.bed`, `.bim`, `.fam`, `.vcf`, are never committed, see `.gitignore`. Only scripts, narrative, summary statistics and figures live in git.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the branch and pull request workflow. Each contributor works inside their own `sections/<name>/` folder and opens a pull request to `main` when ready.

## Citation

If you use this guide, the demonstration dataset or the accompanying scripts, please cite:

> Gentiluomo M, Campa D, Farinella R, et al. How to carry out a genome-wide association study of a rare and complex trait: a decision-focused, reproducible guide using pancreatic ductal adenocarcinoma as a worked example. Target journal to be defined. In preparation.

A citable Zenodo archive of the code and demonstration dataset will be linked here on publication.

## License

Code is released under the [MIT License](LICENSE). Manuscript text and figures are released under CC-BY-4.0.

## Acknowledgement

This work is supported by the TRANSPAN COST Action (CA21116). Maintainers: M. Gentiluomo and R. Farinella, University of Pisa.
