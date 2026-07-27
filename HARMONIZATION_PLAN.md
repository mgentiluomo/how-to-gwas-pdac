# Harmonisation plan: one published, manuscript-aligned repository

**Goal.** Everything the group has produced ends up consolidated, harmonised and published in a
single place, `https://github.com/mgentiluomo/how-to-gwas-pdac`, structured to match the
manuscript's Methods sections. This document is the inventory, the record of decisions already
taken, and the checklist of what remains.

---

## 1. Decisions taken

**Quality control (settled).** The nine-step pipeline developed on `murat_v1` is the canonical QC
of the guide. The five-filter `qc.qmd` at the repository root is superseded and removed. The
pipeline is brought into `main` through a normal, reviewed pull request rather than a wholesale
merge, so that authorship and the review process described in `CONTRIBUTING.md` are preserved.
The manuscript's Table 2, Table 4 and Results section already report the nine-step numbers, and
they match the committed run exactly: 1,461 samples and 430,000 variants at entry, 1,239 samples
and 401,909 variants retained.

**Setup and tooling (settled).** The automated project setup in `scripts/dev/` and the expanded
`getting_started.qmd` replace the manual PLINK download instructions previously on `main`.

**`mgentiluomo/GWAS_training_session` (settled).** It stays a separate, standalone repository. It
is generic, built on a public 1000 Genomes subset rather than on `pdac_demo`, and it does not
carry the didactic design that is the point of this guide. It is not merged in. One piece of it is
worth reusing: `00_Intro/Introduction_from_rad_data_to_binary.md`, the explanation of how raw array
output becomes a PLINK binary dataset, which is what Section 1A needs as a starting point.

---

## 2. Where each remaining section's content comes from

| Manuscript section | Folder | Source | Action |
|---|---|---|---|
| 0, Study design and consortium | `sections/00_introduction/` | Daniele Campa's draft in the manuscript | Transcribe the narrative, no code |
| 1A, Genotyping technologies | `sections/01A_study_design/` | `GWAS_training_session/00_Intro/`, plus the section leads' draft | Adapt the raw-to-PLINK explainer for a non-bioinformatics audience |
| 1B, Quality control | `sections/01B_genotyping_qc/` | Done | Merged from `murat_v1` |
| 2, Population stratification | `sections/02_population_stratification/` | Elif Öz's manuscript draft | Write the walkthrough on `pdac_demo`: LD pruning, PCA, outliers, ancestry inference |
| 3, Imputation | `sections/03_imputation/` | Pelin Ünal's manuscript draft | Write the walkthrough: VCF preparation, TOPMed server submission, post-imputation filtering |
| 4A, Association, binary | `sections/04A_association_binary/` | Elif Öz's manuscript draft | Logistic regression with Firth fallback, QQ and Manhattan plots, lambda |
| 4B, Association, survival | `sections/04B_association_continuous/` | Pelin Ünal's manuscript draft | REGENIE `--t2e` walkthrough |
| 4C, Software comparison | `sections/04C_association_other/` | Murat Güler, writing plan only | Same analysis in two tools, side by side |
| 5, Meta-analysis | `sections/05_meta_analysis/` | Not yet drafted | Described rather than worked: the demo dataset is a single cohort |
| 6, Fine mapping | `sections/06_finemapping_annotation/` | Burçak Otlu's manuscript draft | SuSiE-RSS on the *ABO* locus, credible sets, LD diagnostics |
| 6.5, Functional annotation | `sections/06.5_gwas_annotations/` | Chiara Corradi's manuscript draft | VEP, regulatory overlap, colocalisation |
| 7, Reporting and FAIR | `sections/07_reporting_fair/` | Manuel Gentiluomo | STREGA checklist, summary statistics deposition, Zenodo archiving |

Two chains have to stay consistent as these are written, because the manuscript quotes numbers
from both: the QC-passed dataset feeds Section 2, and the European analysis set defined there,
229 cases and 396 controls, feeds Sections 4A, 6 and 6.5.

---

## 3. Order of operations

1. **Merge the QC pull request** into `main` with the blocking review items resolved, then delete
   the `murat_v1` branch. See `PR_REVIEW_murat_v1.md`.
2. **Apply the alignment commit**: `README.md`, `index.qmd`, `_quarto.yml` and this file.
3. **Run the association step on the QC-passed dataset.** This is now the critical path: it
   produces the genomic inflation factor and the lead variant P value that are the last two
   numeric placeholders in the manuscript's Results section, and it is the input to Sections 2,
   4A and 6 of the tutorial.
4. **Draft Section 1A** from the training-session explainer, with the section leads.
5. **Draft Section 2** on `pdac_demo`, since the ancestry step is what defines the European
   analysis set every later section depends on.
6. **Uncomment each sidebar block in `_quarto.yml`** as its section acquires real content, and
   confirm `quarto preview` renders cleanly before each merge to `main`.

## 4. Repository hygiene, outstanding

* Committed results are heavy in places: one illustration video of 5.9 MB and two full pairwise
  relatedness tables of 4.5 MB and 3.6 MB. The tables are regenerable from the scripts.
* `env/software_versions.md` is still a template. Section 7 of the manuscript commits us to
  reporting exact versions, so it should be completed for the locked-in QC run.
* Repository visibility. The repository is currently readable publicly. That is the right end
  state, but it should be a deliberate choice made now rather than discovered later, since it
  contains results that are not yet published.
* A mapping between the manuscript's `[CODE 1B-nn]` tags and the script filenames should live in
  each section's `README.md`, so a reader of the paper can find the corresponding script.
