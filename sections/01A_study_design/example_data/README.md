# Example raw data

Three Final Reports and an array manifest, all **synthetic**, generated from the
demonstration dataset by `../scripts/00_make_example_raw_data.R` with seed 2026.

## Why synthetic

A real Illumina Final Report contains the genome-wide genotypes of an identifiable
person; a genotype profile is identifying even without a name, and publishing one
requires consent for that specific purpose. A manufacturer's array manifest is not
ours to redistribute either.

Both are therefore generated from `pdac_demo`, whose genotypes derive from the public
HGDP and 1000 Genomes reference panels, in the same way and with the same seed as the
rest of the demonstration data.

## What is real and what is not

Genotypes come from `pdac_demo`. GenCall scores, strand assignments, indel calls,
no-calls, variant names and the manifest are simulated. The parameters were calibrated
on a real Illumina GSAMD-24v3 Final Report so that the artefacts behave as they do in
practice; no part of that file is reproduced here.

## Injected artefacts

| | |
|---|---|
| Variants, samples | 8,000 across DEMO001–DEMO003 |
| GenCall | bimodal, median 0.73, 1.1% below the 0.15 threshold |
| No-calls | about 1.4%, written as `-`, a few with a `NaN` score |
| Indel probes | 7.7%, reported as `I`/`D` |
| Minus-strand design | 42.6% |
| Positional names carrying an older build | 25%, of which 104 also on the wrong chromosome |
| Absent from the manifest | 6% |
| Strand-flipped block in DEMO003 | 200 variants, so that the merge in Step 03 fails and has to be diagnosed |

Regenerate with:

```bash
Rscript scripts/01A_study_design/00_make_example_raw_data.R 8000 2026
```
