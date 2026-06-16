# Section 1B - Genotyping QC

**Lead:** Murat Güler
**Status:** ☐ not started ☐ in draft ☐ ready for review ☐ harmonised

## Website page

The learner-facing tutorial page for this section is:

```text
01B_genotyping_qc.qmd
```

Keep the narrative in that file so the website and section folder stay together.

## What goes in this folder
- `01B_genotyping_qc.qmd` - the website tutorial page.
- `scripts/` - annotated shell and R scripts for this QC section.
- `results/` - small committed outputs used in the documentation, such as PNG plots,
  small TSV summaries, or example text reports.
- Do not commit large genotype files or full runtime result folders.

## Conventions
- Build GRCh38; PLINK 2 + R; heavy comments for a non-bioinformatics audience.
- References as `[DOI/PMID]` in square brackets, highlighted yellow (Vancouver at the end).
- Open a PR to `main` when ready; a maintainer reviews and merges.
