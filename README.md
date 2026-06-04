# How to carry out a GWAS — a step-by-step annotated guide (PDAC case study)

Working repository for the manuscript:

> **"How to carry out a Genome-Wide Association Study: a step-by-step annotated
> guide using pancreatic cancer as a case study"**

Developed within the **TRANSPAN COST Action (CA21116)**, Working Group 1.
The guide targets Master's students and clinicians **without** a bioinformatics
background, and pairs a hybrid review/methods paper with fully annotated
**PLINK / R** code and a small, reproducible demo dataset.

This repository is the single shared base of files for the writing group: text
drafts, annotated code, and figures all live here so everyone works against the
same version.

---

## Repository layout

```
how-to-gwas-pdac/
├── README.md                ← you are here
├── CONTRIBUTING.md          ← how to work in this repo (READ FIRST)
├── LICENSE                  ← MIT (code); see note below for text/figures
├── .gitignore
├── demo_dataset/            ← demo dataset (real HGDP+1KG genotypes, simulated phenotype)
│   ├── README.md
│   └── scripts/             ← 01 … 05 pipeline (download → simulate → GWAS → plots)
├── sections/                ← one folder per manuscript section
│   ├── 00_introduction/
│   ├── 01A_study_design/
│   ├── 01B_genotyping_qc/
│   ├── 02_population_stratification/
│   ├── 03_imputation/
│   ├── 04A_association_binary/
│   ├── 04B_association_continuous/
│   ├── 04C_association_other/
│   ├── 05_meta_analysis/
│   ├── 06_finemapping_annotation/
│   ├── 06.5_gwas_annotations/
│   └── 07_reporting_fair/
├── docs/                    ← figures, the pipeline diagram, shared assets
└── env/                     ← software versions & R package list (reproducibility)
```

Each `sections/<n>/` folder contains a short `README.md` stating the section
lead, the current status, and what belongs there. Put the narrative text as
Markdown and the annotated code alongside it.

---

## Section assignments

| Section | Topic | Lead(s) |
|---|---|---|
| 0 | Introduction / overview | Daniele Campa |
| 1A | Study design & consortium | Bogdan Ungureanu, Anca Riza Costache |
| 1B | Genotyping QC | Murat Güler |
| 2 | Population stratification | Elif Öz |
| 3 | Imputation | Pelin Ünal |
| 4A | Association testing (binary) | Elif Öz |
| 4B | Association testing (continuous) | Pelin Ünal |
| 4C | Association testing (other) | Murat Güler |
| 5 | Meta-analysis | Pelin Ünal |
| 6 | Fine mapping & functional annotation | Burçak Otlu |
| 6.5 | GWAS annotations | Chiara Corradi |
| 7 | Reporting, replication & FAIR data | Manuel Gentiluomo |
| — | Code production & data generation | Riccardo Farinella, Murat Güler |

**Central harmonisation and final review of all code and text: M. Gentiluomo
(with R. Farinella).** Nothing is merged to `main` without maintainer review —
see `CONTRIBUTING.md`.

---

## The demo dataset

Everything you need to build the running example is in `demo_dataset/`. It
uses **real HGDP + 1000 Genomes genotypes** (anonymised) with a fully simulated
PDAC phenotype driven by a single well-replicated causal locus (**ABO**, 9q34.2).
The data are downloaded from the GitHub Release via `bash demo_dataset/download_data.sh`. See `demo_dataset/README.md` for the full design, teaching-dataset caveats, and provenance.

---

## Conventions

- **Code:** PLINK 1.9 (`plink`) and PLINK 2.0 (`plink2`), plus R; one self-contained, heavily commented script per step.
- **Documentation:** Markdown.
- **References:** provisional placeholders as `[DOI/PMID]` in square brackets,
  highlighted yellow; final formatting in **Vancouver** style.
- **Build:** GRCh38 throughout.
- **No data in git:** raw/intermediate genotype files stay out of the repo
  (see `.gitignore`); only scripts, final summary statistics, and figures are
  committed.

---

## License

Code is released under the **MIT License** (see `LICENSE`). Manuscript text and
figures, if hosted here, are intended for **CC-BY-4.0**; confirm with the
journal (Human Genomics, open access) before final release.

---

## Acknowledgement

This work is supported by the **TRANSPAN COST Action (CA21116)**.
Maintainers: M. Gentiluomo & R. Farinella, University of Pisa.
