# Project state

Last updated 28 July 2026. This file is the single place where the current state
of the project is written down. It is deliberately blunt about what is finished,
what is half-done, and what is known to be wrong, because the alternative is that
someone opens the repository and cannot tell.

---

## The short version

The repository is currently in a **fork**. Two mutually inconsistent versions of
the demonstration data exist, and only one of them is published.

| | Published site and manuscript | New dataset, not yet published |
|---|---|---|
| Samples after QC | 1,217 | 1,214 |
| Variants after QC | 414,695 | 414,415 |
| European analysis set | 224 cases / 387 controls | 219 / 390 |
| Causal variants | one, at MAF 0.086 | two, at MAF 0.377 and 0.489 |
| Generative odds ratio | documented 2.40, **actually ~3.4** | 2.60 and 1.50, **verified** |
| ABO credible set | 77 in EUR, 1 in meta | 11 in EUR, 9 in meta |

Everything currently on <https://mgentiluomo.github.io/how-to-gwas-pdac/> and in
manuscript v9.1 describes the left column and is internally consistent. The right
column exists as a validated simulation script and a completed pipeline run held
outside the repository, plus one section already committed here but kept out of
the sidebar.

**Do not publish anything from the right column until the whole chain is
converted, or the site will describe two datasets at once.**

---

## Why the dataset is being replaced

The demonstration dataset documented a causal odds ratio of 2.40 at *ABO*. It
does not contain one. Estimated from the two ancestry strata that were never used
for discovery, and therefore free of selection bias, the effect is **3.65 (95% CI
2.63 to 5.09)**; a model-free 2x2 allelic table gives 3.37, and the estimate is
unchanged with or without covariate adjustment, so non-collapsibility does not
explain it. The value 2.40 lies outside the confidence interval of every
unbiased estimate.

The most likely cause is that the effect was specified on a liability scale and
estimated on the log-odds scale. The ratio of the two on the log scale is about
1.36, which is the usual probit-to-logit conversion factor and matches what is
observed. The original simulation script is not in the repository, so this cannot
be confirmed from the code.

**What it invalidated.** The manuscript reported a winner's curse of 62%,
measured as 3.90 against a supposed truth of 2.40. Against the real effect the
selection bias is about 7%. The power section reported 9.2% power at the causal
variant; at the true effect it is roughly 74%. The claim that the locus was
"detected despite being underpowered" was therefore wrong: detection was
expected.

**How it was found.** It was not found by internal review. It was raised by an
external model asked to act as a hostile referee, and confirmed by direct
measurement. The failure was accepting a number from the dataset documentation
without ever comparing it against an estimate from the data it described.

**The rule this produces**, which belongs in Table 1 and in the reporting
section: *when a simulated dataset is used as ground truth, verify that the
analysis recovers the generative parameter before treating it as truth.*

---

## What was decided today

**Hardy-Weinberg and heterozygosity are computed within ancestry.** Both tests
assume a single randomly mating population. On the pooled cohort, HWE excluded
14,378 variants where the within-ancestry test finds five; heterozygosity found
two outliers where the within-ancestry calculation finds twenty. The two filters
fail in opposite directions and neither failure is visible in the final counts.
The pooled HWE test is retained before ancestry assignment as a diagnostic that
excludes nothing; the exclusion happens after the split.

This was verified quantitatively, not merely argued. Across 415,520 variants,
F_ST between the three groups predicts the observed HWE statistic with a
correlation of 0.85 and a regression slope of 0.94, and explains 69.6% of its
variance. Of 282,000 variants with F_ST below 0.05, eight are expected to fail.
The Wahlund effect accounts for the exclusions. **This verification is not yet in
the manuscript and should be added: it converts the most attackable claim in
Section 4 into the best documented one.**

**The new dataset has two causal variants**, so that the guide can show a true
positive recovered and a true positive missed with the truth known for both:

- `9:133273682:A:T`, *ABO* 9q34.2, MAF 0.377, generative OR **2.60**. Detected,
  *P* = 4.6e-9, estimated OR 2.33. About 99% probability of detection by design;
  2.30 was tried first and gave 91%, and the first draw fell in the missing 9%.
  The design parameter was raised, not the seed, because selecting a seed on
  significance is conditioning on significance.
- `5:1286401:C:A`, TERT/CLPTM1L 5p15.33, MAF 0.489, generative OR **1.50**. Not
  detected in Europeans (*P* = 0.15); reaches *P* = 0.03 when the three strata
  are combined, which is the meta-analysis argument made on a known true effect.

**The simulation validates itself.** It re-estimates every generative parameter
with the model the tutorial uses and refuses to write the dataset if any estimate
is inconsistent with its target. Survival is now defined for cases only, as time
from diagnosis, with staggered recruitment giving 22% censoring instead of the
previous 99.8% event rate. Sex and age effects use literature values (OR 1.30 for
male sex, 1.50 per decade).

**The fine-mapping result changed and is being reported honestly.** With a common
causal variant the European credible set is already narrow, 11 variants over
12 kb, and the meta-analysis reduces it only to 9. What the meta-analysis does is
move the causal variant from third place to first. The dramatic collapse from 77
to 1 was a property of the old, rarer causal variant and is not reproducible with
a common one. Choosing a causal variant to recover that result would be selecting
an ingredient of the experiment on the basis of how the experiment comes out.

---

## New content written today

**`sections/07_linkage_disequilibrium/`** — committed, deliberately **not** in
`_quarto.yml`, because it carries the new dataset's numbers. It measures LD decay
per ancestry in the ABO region and shows that 11 variants are indistinguishable
from the causal one in Europeans against 1 in Africans, which is why the credible
set is the width it is. It argues that ancestral diversity buys resolution and
not only equity. References verified: Rogers 2014 (doi:10.1534/genetics.114.166454),
Reich 2001 (PMID:11346797), Ardlie 2002 (PMID:11967554), Wall & Pritchard 2003
(doi:10.1038/nrg1123), Chapman & Thompson 2001 (PMID:11037333).

**`scripts/dev/make_checkpoints.sh`** and **`docs/helpers/checkpoint_entry_block.qmd`** —
the checkpoint generator and the two-part entry check. Checksums prove integrity;
expected counts prove provenance. Deliberately not claimed: that re-running the
pipeline reproduces the same bytes. It does not, and for the eigenvector files it
cannot, because eigenvector sign is arbitrary.

**`docs/helpers/code_tab_format_sample.qmd`** — the agreed format for R and
Python code tabs: PLINK steps call PLINK identically in both tabs, only the
analysis layer is genuinely twinned.

---

## Releases

| Tag | Contents | Status |
|---|---|---|
| `v0.1-data` | pdac_demo genotypes, phenotype, covariates, survival | **Will be superseded.** The phenotype, covariates and survival files change with the new simulation; the genotypes do not |
| `v0.2-checkpoints` | summary-statistics bundle, 8 files + manifest | Published and verified end to end. **Will need regenerating** for the same reason |

Pages must link to an explicit tag, never to `latest`, so that a page and the
data it was written against move together.

---

## What remains

**To convert the fork** (a full day, comparable to today):

1. Publish the new phenotype files under a new data tag.
2. Regenerate figures, result files and the seven existing pages.
3. Add the LD section to `_quarto.yml`.
4. Update the manuscript: the new numbers, the corrected winner's curse, the
   Wahlund verification, and the honest fine-mapping result.
5. Regenerate the checkpoint bundle under a new tag.

**Independent of the fork:**

- Section 6.5, functional annotation. Needs five minutes of VEP on the credible
  set; supplies the annotation that Table 3 requires, that table having been
  deleted from the manuscript for lack of it.
- Sections 0, 3, 4B, 4C and 7 of the tutorial, currently described in the
  manuscript but not demonstrated.
- The R and Python code tabs across the existing pages.
- A limitation the manuscript does not yet state: the demonstration has three
  cleanly separated continental groups, so within-ancestry filtering is easy. In
  a real, predominantly European PDAC cohort with continuous gradients and
  admixed individuals, the central recommendation of Section 4 is much harder to
  apply. RUTH is mentioned in passing but the easy case is not admitted.

**Deadline.** Acceptance is required by mid-October 2026 for the COST Action
budget. The fork conversion and the missing sections do not both fit comfortably
before submission.

---

## Standing practices

- References are verified before citation. Two fabricated entries were caught in
  an AI-generated bibliography early in the project; the 2.40 error is the same
  failure applied to a number instead of a citation.
- Numbers in the manuscript, the site and the released files must agree. Check by
  scanning for the old values, not by trusting that an edit was applied.
- Work outside OneDrive. Synchronisation corrupts `.git`.
- Windows PowerShell mangles non-ASCII characters pasted into the console. Put
  replacements involving them in a `.ps1` file instead.
- Repository files use CRLF. Multi-line string matching must normalise line
  endings or it silently fails.
