# Linkage disequilibrium and what it decides

**Lead:** Manuel Gentiluomo

Why the resolution of a fine-mapping analysis is set by the population being
studied and not by the method chosen. Placed after fine mapping because it
explains a result the reader has already seen: a credible set of eleven variants
that no amount of statistical care could have made smaller in this sample.

## Contents

- `07_linkage_disequilibrium.qmd` — the website tutorial page.
- `scripts/01_ld_decay.sh` — one step, PLINK 2 only.
- `results/` — outputs of the canonical run.

## Key numbers from the demonstration run

Variants that cannot be told apart from the causal variant at *ABO*, by ancestry
group, from `pdac_demo_07_ld_causal.tsv`:

| Group | Pairs tested | r² above 0.8 | r² above 0.5 | Highest r² |
|---|---|---|---|---|
| European | 127 | 11 | 12 | 0.983 |
| East Asian | 123 | 11 | 11 | 1.000 |
| African | 132 | **1** | 10 | 0.985 |

Mean r² by distance, from `pdac_demo_07_ld_decay.tsv`:

| Distance | European | East Asian | African |
|---|---|---|---|
| 0 to 10 kb | 0.271 | 0.344 | 0.184 |
| 10 to 25 kb | 0.121 | 0.173 | 0.062 |
| 25 to 50 kb | 0.047 | 0.089 | 0.022 |
| 50 to 100 kb | 0.031 | 0.029 | 0.010 |
| 100 to 200 kb | 0.017 | 0.013 | 0.006 |
| 200 to 300 kb | 0.007 | 0.005 | 0.005 |

## Why this section carries weight

It supplies the reason behind the fine-mapping result rather than repeating it.
Eleven variants at *ABO* are statistically indistinguishable from the causal one
in Europeans; in Africans there is one. The credible set in Section 6 could not
have been narrowed by a better prior, a different method or more careful work,
because the information needed to separate those variants is not present in a
European sample. It is present in an African one.

That is also the honest argument for ancestral diversity in genetic studies, and
it is a methodological argument rather than an ethical one: shorter haplotypes
carry more information per variant, so the same number of people resolves a locus
further. The ethical argument stands on its own; this section shows that the two
point the same way.

## Conventions

- PLINK 2 only, no packages to install.
- r² computed within ancestry group, since LD is a property of a population and
  pooling groups measures neither.
- Distance bins are reported with the number of variant pairs behind each, so a
  reader can see which estimates rest on few pairs.
