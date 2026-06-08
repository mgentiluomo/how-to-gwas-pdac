# How to carry out a GWAS — companion website (Quarto)

Quarto source for the website accompanying the TRANSPAN WG1 methodology manuscript
*"How to carry out a Genome-Wide Association Study: a step-by-step annotated guide using
pancreatic cancer as a case study."*

## Requirements

- [Quarto](https://quarto.org) (bundled with recent RStudio)
- R + RStudio (for the executable R chunks)

## Preview in RStudio

Open the folder as a project and click **Render**, or in the *Terminal* tab:

```bash
quarto preview      # live preview, reloads on save
```

## Build

```bash
quarto render       # builds the static site into _site/
```

## Structure

```
_quarto.yml            # project + website config, sidebar navigation
index.qmd              # home page
getting_started.qmd    # Before you start: download + PLINK setup + first --freq test
styles.css             # code-tag pills and PMID highlight
```

Manuscript sections (QC, population stratification, ...) are added as further `.qmd` files
and listed in the `sidebar` section of `_quarto.yml`.

## Notes

- Update `repo-url` in `_quarto.yml` to the final University of Pisa GitHub repo.
- PLINK commands are shown for the reader to run locally; the site does not re-execute them.
  Light R chunks use Quarto `freeze` (run once, cached).
