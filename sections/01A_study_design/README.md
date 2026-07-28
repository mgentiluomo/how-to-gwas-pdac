# Section 1A — Genotyping technologies and file formats

**Lead:** Bogdan Ungureanu, Anca Riza Costache

From raw array output to an analysable PLINK dataset: what a Final Report contains,
how it is converted, which checks to run before anything else, and why a cohort
genotyped on several platforms needs a decision taken at the design stage.

## Contents

- `01A_study_design.qmd` — the website tutorial page.
- `scripts/` — four steps, PLINK 1.9, PLINK 2 and base R.
- `example_data/` — a real Illumina Final Report excerpt (see its README).
- `results/` — outputs of one locked-in run.

## Manuscript code tags

| Tag | Script |
|---|---|
| CODE 1A-00 | `scripts/00_make_example_raw_data.R` (builds the example raw data) |
| CODE 1A-01 | `scripts/01_inspect_final_report.sh` |
| CODE 1A-02 | `scripts/02_report_to_lgen.R` |
| CODE 1A-03 | `scripts/03_lgen_to_binary.sh` |
| CODE 1A-04 | `scripts/04_preflight_and_formats.sh` |

## Key numbers from the demonstration run

| | |
|---|---|
| Example data | 8,000 variants, 3 samples, synthetic, seed 2026 |
| Median GenCall | 0.73; 1.1% below the 0.15 failure threshold |
| Manifest match | 7,520 of 8,000 (94.0%); the remaining 6% are absent by construction |
| Positional names vs manifest coordinates | 1,860 positional names, none agreeing; 97 also differ in chromosome |
| Design strand | 3,194 of 7,520 (42.5%) on the minus strand, alleles complemented |
| Allele column contents, DEMO001 | 7,292 SNP calls, 102 no-calls, 606 indels dropped and counted |
| Merge | DEMO002 clean; DEMO003 conflicts on 142 variants, resolved by flipping |
| Unobserved alleles repaired from the manifest | 1,546 of 1,546 |
| Final merged dataset | 6,938 variants, 3 samples |

## Measurements quoted from a real array

These come from a genuine Illumina GSAMD-24v3 Final Report and are reported as
measurements; the file is not distributed.

| | |
|---|---|
| GenCall by marker type | rs-named 0.720 median with 1.2% failures; position-named custom and CNV content 0.441 median with 4.5% failures |
| Variants with positional names | 12,192 of 730,059 (1.7%) |
| Positional names agreeing with the manifest coordinate | 776 of 5,745; 4,969 differ by 384–543 kb, 75 also in chromosome |
| Minus-strand design | 43% |
| Indel probes | 7.7% of an 8,000-variant sample |
| Duplicated positions | 255 in 6,677 converted variants |

## Resources linked on this page

- Illumina DRAGEN Array, the current secondary-analysis path for Infinium arrays:
  <https://www.illumina.com/products/by-type/informatics-products/microarray-software/dragen-array.html>
- Broad Institute Birdsuite FAQ, for legacy Affymetrix 5.0/6.0 data: Birdseed versus
  Birdsuite calls, and the per-plate clustering requirement:
  <https://www.broadinstitute.org/birdsuite/birdsuite-faq>

Both are product or project pages rather than peer-reviewed work, so they are linked in
the tutorial and should not enter the manuscript reference list as citations.

## Verified references used on this page

- Verlouw JAM, Clemens E, de Vries JH, et al. A comparison of genotyping arrays.
  *Eur J Hum Genet* 2021;29:1611–1624. doi:10.1038/s41431-021-00917-7
- Johnson EO, Hancock DB, Levy JL, et al. Imputation across genotyping arrays for
  genome-wide association studies: assessment of bias and a correction strategy.
  *Hum Genet* 2013;132:509–522. doi:10.1007/s00439-013-1266-7, PMID:23334152

## Conventions

- Build GRCh38; PLINK 1.9 and 2, plus base R, no packages to install.
- Step 02 and 03 need PLINK 1.9 because the long format is a PLINK 1 input; every
  other section uses PLINK 2 only.
- Open a PR to `main` when ready; a maintainer reviews and merges.
