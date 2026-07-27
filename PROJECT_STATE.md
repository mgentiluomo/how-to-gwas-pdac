# Project state

The working record of this project: what has been decided, what the numbers are and
which script produces each of them, what has been verified, and what remains. It is
kept in the repository so that it survives interruptions, so that a co-author can pick
up the thread without reconstructing it, and because a methodology paper about
reproducibility should be able to account for its own.

Last updated: 27 July 2026.

---

## 1. What this project is

A methodology manuscript for *Human Genomics* and its companion tutorial website,
developed within the TRANSPAN COST Action (CA21116), Working Group 1.

**Title.** How to carry out a genome-wide association study of a rare and complex
trait: a decision-focused, reproducible guide using pancreatic ductal adenocarcinoma
as a worked example.

**The argument.** GWAS methodology for rare, genetically complex traits requires
explicit decisions at every step, not the defaults designed for large, well-powered
common-disease studies. The manuscript supplies the reasoning; the tutorial supplies
the commands; the demonstration dataset makes both verifiable.

**Editorial principle.** The manuscript narrative carries the considerations specific
to rare and complex traits, meaning why a choice is appropriate, when an analysis is
legitimate and whether a signal can be interpreted. Generic instruction belongs in the
tutorial. These are not interchangeable and should not be mixed.

**Writing conventions.** UK English. Punctuation limited to the full stop, comma,
semicolon and colon. The word "framework" is not used. References are provisional as
`[DOI:...]` or `[PMID:...]`, highlighted yellow, and converted to Vancouver at
submission.

---

## 2. Decisions taken

| Decision | Resolution |
|---|---|
| Canonical quality control | The nine-step pipeline, merged from `murat_v1`. The earlier five-filter `qc.qmd` is superseded and removed. |
| Demonstration data | Everything runs on `pdac_demo` or on material generated from it. No real participant genotypes are distributed, anywhere. |
| Array manifest | Not redistributed. Section 1A requires the reader to supply the manifest matching their own product, and verifies the match before converting. |
| Repository | `mgentiluomo/how-to-gwas-pdac`, public, single source of truth. `GWAS_training_session` remains separate; two pieces of it were adapted. |
| Site build | Quarto, published to `gh-pages` with `quarto publish`. `_quarto.yml` restricts rendering to `.qmd` and `.md`, so section scripts are never executed as documents. |
| Association testing | Autosomes only. The simulated chromosome X exists for the sex check and carries no phenotype signal. |
| GenCall filtering | Applied at conversion, threshold as a parameter, failures written as missing rather than deleted, counts logged. |
| Indels | Dropped, deliberately and counted, not filtered away by a condition written for SNPs. |
| Fine mapping method | Single causal variant approximate Bayes factor, written out in full, with the assumption stated and multi-signal methods signposted. |

---

## 3. Key numbers, and the script that produces each

Every figure below is reproducible from the released data and the committed scripts.

### Quality control, `sections/01B_genotyping_qc/scripts/`

| | |
|---|---|
| Entering | 1,461 samples, 430,000 variants |
| Sample call rate `--mind 0.02` | 25 removed |
| Sex check | 4 removed |
| Heterozygosity, F ± 3 SD | 2 removed |
| Variant call rate `--geno 0.05` | 8,480 removed |
| Hardy-Weinberg in controls, `1e-6` | 14,573 removed |
| Relatedness, KING > 0.1875, phenotype-aware | 191 removed, 154 controls and 37 cases |
| MAF `0.01` | 5,038 removed |
| Retained | 1,239 samples, 401,909 variants; 84.8% and 93.5% |

Independently reproduced on a separate machine with a separately compiled PLINK 2,
returning identical counts at every step.

### Ancestry and analysis set, `sections/02_population_stratification/scripts/`

| | |
|---|---|
| After LD pruning | 221,352 of 401,909 variants |
| Variance explained, PC1 / PC2 / PC3 | 47.0% / 38.1% / 3.7% |
| Composition after QC | 625 EUR, 319 EAS, 295 AFR |
| European analysis set | 229 cases, 396 controls; effective N 580.4 |
| Variance explained within EUR, PC1 / PC2 | 23.9% / 11.8% |
| Stability check | lead variant P 5.31 × 10⁻⁹ with ten components, 5.41 × 10⁻⁹ with five |

### Association, `sections/04A_association_binary/scripts/`

| | |
|---|---|
| Autosomal variants tested | 396,396; 392,331 with a valid P |
| Firth fallback invoked | 2,221 variants |
| Genomic inflation | λ = 1.028, λ₁₀₀₀ = 1.098 |
| P < 5 × 10⁻⁸ | 1 variant |
| P < 1 × 10⁻⁵ | 8 variants |
| Lead variant | `9:133249045:A:G`, OR 3.80 (95% CI 2.43–5.95), P 5.3 × 10⁻⁹ |
| Simulated truth at that variant | OR 2.40; the estimate is inflated by about 58%, the winner's curse |
| Power at that variant, European set | approximately 10%; detection in an underpowered scan is itself a selection event |

### Meta-analysis, `sections/05_meta_analysis/scripts/`

| | |
|---|---|
| Strata | EUR 229/396, AFR 142/153, EAS 141/178 |
| Effective N, combined | 1,189.7 against 580.4 |
| Effect-allele flips required | 63,707, within a single genotype file |
| Strand-ambiguous variants | 1,833 (0.6%), 339 unresolvable by frequency |
| Variants present in all three strata | 296,678 of 392,331 |
| λ per stratum / meta | 1.022, 0.991, 1.013 / 1.003 |
| Lead variant | P 2.4 × 10⁻²², OR 3.74 (2.87–4.87), I² = 0%, P het 0.44 |
| Heterogeneity at LD proxies | up to I² = 79%; one variant falls from P 1.3 × 10⁻⁹ to 0.0035 under random effects |
| Power at the causal variant, combined | approximately 69% |

### Raw data handling, `sections/01A_study_design/scripts/`

Measured on synthetic data generated from `pdac_demo`:

| | |
|---|---|
| Manifest match | 7,520 of 8,000 (94.0%) |
| Positional names agreeing with the manifest coordinate | 0 of 1,860; 97 also differ in chromosome |
| Minus-strand design | 3,194 of 7,520 (42.5%) |
| Merge | DEMO002 clean; DEMO003 conflicts on 142 variants, resolved by flipping |
| Unobserved alleles repaired from the manifest | 1,546 of 1,546 |

Measured on a real Illumina GSAMD-24v3 report, quoted but not distributed:

| | |
|---|---|
| GenCall by marker type | rs-named 0.720 median, 1.2% failures; position-named custom content 0.441 median, 4.5% failures |
| Variants with positional names | 12,192 of 730,059 (1.7%) |
| Positional names agreeing with the manifest | 776 of 5,745; the rest differ by 384–543 kb |
| Minus-strand design | 43% |
| Indel probes | 7.7% |

### Fine mapping, `sections/06_finemapping_annotation/scripts/`

| | |
|---|---|
| Region | chr9:132.9–133.6 Mb, 198 variants |
| European set, credible set | 55 variants; top variant PP 0.51, the true causal variant, ranked first |
| Three strata combined, credible set | **1 variant**; PP 0.99, the true causal variant |

---

## 4. Verified references

Checked against the publisher record. Others in the manuscript remain to be verified.

| Reference | Identifier |
|---|---|
| Verlouw JAM, Clemens E, et al. A comparison of genotyping arrays. *Eur J Hum Genet* 2021;29:1611–1624 | doi:10.1038/s41431-021-00917-7 |
| Johnson EO, Hancock DB, et al. Imputation across genotyping arrays for GWAS. *Hum Genet* 2013;132:509–522 | doi:10.1007/s00439-013-1266-7, PMID:23334152 |
| Klein AP, Wolpin BM, et al. Genome-wide meta-analysis identifies five new susceptibility loci for pancreatic cancer. *Nat Commun* 2018 | PMID:29422604, GWAS Catalog GCST005434 |
| Martin AR, Kanai M, et al. Clinical use of current polygenic risk scores may exacerbate health disparities. *Nat Genet* 2019;51:584–591 | doi:10.1038/s41588-019-0379-x |

**Corrections made.** A DOI supplied for Martin 2019 pointed to the 2021 publisher
correction rather than the article. Two references in an earlier AI-generated
bibliography were found not to exist and were removed. Verify before citing.

**Linked but not cited.** Illumina DRAGEN Array and the Broad Institute Birdsuite FAQ
are product and project pages, linked in the tutorial and not entered in the reference
list.

---

## 5. State of the deliverables

### Tutorial: six sections of twelve

| Section | State |
|---|---|
| 1A Genotyping technologies and file formats | complete |
| 1B Quality control | complete |
| 2 Population stratification | complete |
| 4A Association testing | complete |
| 5 Meta-analysis | complete |
| 6 Fine mapping | complete |
| 0 Study design | not started, narrative only |
| 3 Imputation | not started, partly demonstrable |
| 4B Survival GWAS | not started, demonstrable, the dataset has `survival.txt` |
| 4C Software comparison | not started, partly demonstrable |
| 6.5 Functional annotation | not started, narrative only |
| 7 Reporting and FAIR | not started, partly demonstrable |
| Statistical power | proposed, all numbers available |

### Manuscript

At version 6. No numeric placeholders remain. Outstanding: Figure 1, the workflow
diagram; Figure 5, fine mapping, now producible from Section 6; the Zenodo DOI,
authors' contributions and acknowledgements; conversion of references to Vancouver;
verification of the DOIs supplied for the power section.

Figures 2, 3 and 4 are produced by `docs/make_manuscript_figures.R` from committed
results.

### Open items

- Whether WSL2 remains a hard requirement for Windows readers. Raised with Murat
  Güler, unresolved. It affects who can follow the guide at all.
- `env/software_versions.md` is still a template and Section 7 of the manuscript
  commits the project to reporting exact versions.
- Repository weight: one 5.9 MB illustration video and two regenerable relatedness
  tables of 4.5 and 3.6 MB.
- The GWAS Catalog association file for GCST005434, needed for the replication lookup
  in Section 7 of the tutorial.
- The manuscript states in Results that the trans-ancestry credible set was not
  smaller than the single-ancestry set. Section 6 measures the opposite, 55 variants
  against 1. The Results text needs updating to match.

---

## 6. Practical notes

- The demonstration data are not in git. Fetch with `bash demo_dataset/download_data.sh`
  or from the release tag `v0.1-data`.
- Work outside OneDrive. Synchronisation corrupts `.git` directories.
- All random steps use seed 2026.
- A pull request is merged only when GitHub shows the label **Merged**. Opening one is
  not enough, and a fast-forward `git merge` locally does not update the remote branch.
- Quarto renders only `.qmd` and `.md` by configuration. Scripts placed in section
  folders are ignored, which is deliberate.
