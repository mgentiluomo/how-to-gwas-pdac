# Project state

Last updated 28 July 2026. This file is the single place where the current
state of the project is written down. It is deliberately blunt about what is
finished, what is half-done, and what is known to be wrong, because the
alternative is that someone opens the repository and cannot tell.

---

## The short version

The repository is in a **fork**. Two mutually inconsistent versions of the
demonstration data exist, and only one of them is published.

| | Published site, manuscript and repository | New dataset, not yet published |
|---|---|---|
| Samples after QC | 1,217 | 1,214 |
| Variants after QC | 414,695 | 414,415 |
| European analysis set | 224 cases / 387 controls | 219 / 390 |
| Causal variants | one, at MAF 0.086 | two, at MAF 0.377 and 0.489 |
| Generative odds ratio | documented 2.40, **actually ~3.65** | 2.60 and 1.50, **verified** |
| ABO credible set | 77 in EUR, 1 in meta | 11 in EUR, 9 in meta |

Everything on the site, in manuscript v11 and in this repository describes the
left column and is internally consistent, verified across all versioned text
files. The right column exists as a validated simulation script and a completed
pipeline run held outside the repository, plus one section committed here but
excluded from rendering.

**Do not publish anything from the right column until the whole chain is
converted, or the resource will describe two datasets at once.**

---

## Why the dataset is being replaced

The demonstration dataset documented a causal odds ratio of 2.40 at *ABO*. It
does not contain one. Estimated from the two ancestry strata never used for
discovery, and therefore free of selection bias, the per-allele effect is
**3.65 (95% CI 2.63 to 5.09)**; a model-free 2x2 allelic table gives 3.37, and
the estimate is unchanged with or without covariate adjustment, so
non-collapsibility does not explain it. The value 2.40 lies outside the
confidence interval of every unbiased estimate.

The cause is documented in the dataset's own README, which states that the
phenotype was generated under a **liability-threshold model** while the analysis
estimates an odds ratio on the log-odds scale. The two scales differ by roughly
1.3 to 1.8, and the observed ratio of log-effects is 1.36. It was written down
from the beginning; nobody compared it against an estimate from the data.

**What it invalidated.** The manuscript reported a winner's curse of 62%,
measured as 3.90 against a supposed truth of 2.40. Against the real effect the
selection bias is about 7%. The power section reported 9.2% power at the causal
variant; at the true effect it is roughly 74%. The claim that the locus was
"detected despite being underpowered" was wrong: detection was expected.

**The rule this produces**, which belongs in Table 1 and in the reporting
section: *when a simulated dataset is used as ground truth, verify that the
analysis recovers the generative parameter before treating it as truth.*

---

## What was decided

**Hardy-Weinberg and heterozygosity are computed within ancestry.** On the
pooled cohort, HWE excluded 14,378 variants where the within-ancestry test finds
five; heterozygosity found two outliers where the within-ancestry calculation
finds twenty. The two filters fail in opposite directions and neither failure is
visible in the final counts. The pooled test is retained before ancestry
assignment as a diagnostic that excludes nothing.

This is verified quantitatively, and the verification is now in the manuscript.
Across 415,520 variants, F_ST between the three groups predicts the observed HWE
statistic with a correlation of 0.85, a regression slope of 0.94 against the
1.00 of exact agreement, and 69.6% of its variance. Among 282,000 variants with
F_ST below 0.05, eight failures are expected.

**The new dataset has two causal variants**, so the guide can show a true
positive recovered and a true positive missed with the truth known for both:

- `9:133273682:A:T`, *ABO* 9q34.2, MAF 0.377, generative OR **2.60**. Detected,
  *P* = 4.6e-9, estimated 2.33. Designed for about 99% detection probability;
  2.30 was tried first, gave 91%, and the first draw fell in the missing 9%. The
  design parameter was raised, not the seed, because selecting a seed on
  significance is conditioning on significance.
- `5:1286401:C:A`, TERT/CLPTM1L 5p15.33, MAF 0.489, generative OR **1.50**. Not
  detected in Europeans (*P* = 0.15); reaches *P* = 0.03 across the three strata,
  which is the meta-analysis argument made on a known true effect.

**The simulation validates itself.** It re-estimates every generative parameter
with the model the tutorial uses and refuses to write the dataset if any
estimate is inconsistent with its target. Survival is defined for cases only, as
time from diagnosis, with staggered recruitment giving 22% censoring instead of
the previous 99.8% event rate.

**The fine-mapping result changed and is reported honestly.** With a common
causal variant the European credible set is already narrow, 11 variants over
12 kb, and the meta-analysis reduces it only to 9. What the meta-analysis does
is move the causal variant from third place to first. The collapse from 77 to 1
was a property of the old, rarer causal variant. Choosing a causal variant to
recover that result would be selecting an ingredient of the experiment on the
basis of how the experiment comes out.

---

## The public resource was inconsistent, and is no longer

An external hostile review found that the repository's own documentation stated
figures from two superseded runs, and judged this cause for major revision
bordering on rejection. It was right. A repository-wide scan then found nine
further files the review had not cited. All are corrected:

- `README.md`: the status section said the sections were scaffolding when seven
  were published. Now lists eight published stages and five as **not built**
  rather than forthcoming.
- `demo_dataset/README.md`: stated 254/508, lambda 1.01 and *P* 2e-10 from a run
  two versions old. Now carries the canonical expected results and an explicit
  caveat that the liability-scale construction means the observed odds ratio is
  not the specified parameter.
- `HARMONIZATION_PLAN.md`: **removed** from the public branch. It asserted
  1,239 samples and 401,909 variants as canonical, and it exposed section
  assignments by name on a public repository. It remains in git history.
- Two committed result files, two scripts that regenerate them, and five
  section READMEs carried superseded figures. Corrected in both the results and
  the sources, because fixing a result without its script lets the wrong number
  reappear on the next run.
- Internal progress markers removed from all fourteen section READMEs.
  Attribution lines stay: the names are on the paper and the credit is due.
  Work status is coordination and does not belong on a public repository.

**The lesson, and it cost three repetitions to learn.** Verify **every versioned
file**, not the published artefact. Section READMEs and committed result files
never reach the rendered HTML, but they are the first thing a reviewer opens on
GitHub. The check is one line:

```bash
grep -rlE "<old figure>" --include='*.md' --include='*.qmd' --include='*.R' \
     --include='*.sh' --include='*.txt' .
```

---

## Content added

**`sections/07_linkage_disequilibrium/`** — committed, deliberately excluded
from `_quarto.yml`, because it carries the new dataset's numbers. Measures LD
decay per ancestry and shows 11 variants indistinguishable from the causal one
in Europeans against 1 in Africans, which is why the credible set is the width
it is. Argues that ancestral diversity buys resolution and not only equity.
References verified: Rogers 2014 (doi:10.1534/genetics.114.166454), Reich 2001
(PMID:11346797), Ardlie 2002 (PMID:11967554), Wall and Pritchard 2003
(doi:10.1038/nrg1123), Chapman and Thompson 2001 (PMID:11037333).

**Manuscript v11** adds: the abstract distinguishing worked from described
stages; the Wahlund verification in Section 4; and **Table 5**, a sensitivity
analysis over the relatedness cutoff, the HWE threshold, the MAF filter and the
number of principal components. The result holds in every configuration. The
paragraph gives more space to the two things the reassuring headline hides: the
largest movement comes from the number of principal components, the least
formalised decision in the pipeline; and the relatedness threshold barely
matters here only because the source panels contain trios and duplicates, which
does not generalise.

**`scripts/dev/make_checkpoints.sh`** and the entry-check template: checksums
prove integrity, expected counts prove provenance. Deliberately not claimed:
that re-running reproduces the same bytes. It does not, and for eigenvector
files it cannot, because eigenvector sign is arbitrary.

---

## Releases

| Tag | Contents | Status |
|---|---|---|
| `v0.1-data` | genotypes, phenotype, covariates, survival | **Will be superseded**: the phenotype files change with the new simulation, the genotypes do not |
| `v0.2-checkpoints` | summary-statistics bundle, 8 files + manifest | Published, verified end to end. **Will need regenerating** |

Pages link to an explicit tag, never to `latest`, so a page and the data it was
written against move together.

---

## What remains, in order

1. **Convert the fork.** A full day. New data tag; regenerate figures, result
   files and the seven pages; add the LD section to `_quarto.yml`; update the
   manuscript; regenerate the checkpoint bundle. Every file corrected today
   carries the current run's figures and will need doing again, but the list is
   now known and the `grep -r` finds them.
2. **Methodological refinements** raised by review and still open:
   stratum-specific rather than global HWE exclusion; heterozygosity as a flag
   rather than automatic removal; the five HWE failures described as flagged for
   investigation rather than demonstrated errors, since 1.24 are expected by
   chance at that threshold across three tests; the lambda claim narrowed to what
   it supports; I-squared with three strata presented as underpowered.
3. **Freeze and deposit before review, not after acceptance.** A mutable GitHub
   repository is not a reproducible artefact. Versioned release, Zenodo DOI in
   the manuscript, commit SHA, checksums, environment lockfile.
   `env/software_versions.md` is still a template and the manuscript promises
   exact versions.
4. **Sections not built**: functional annotation (needs one VEP run on the
   credible set, and supplies the table deleted from the manuscript for lack of
   it), study design, imputation, time-to-event, software comparison, reporting.
5. **A limitation not yet stated**: the demonstration has three cleanly
   separated continental groups, so within-ancestry filtering is easy. In a real,
   predominantly European cohort with continuous gradients and admixed
   individuals, the central recommendation of Section 4 is much harder to apply.

**Deadline.** Acceptance is required by mid-October 2026 for the COST Action
budget. Items 1 and 4 do not both fit comfortably before submission.

---

## Standing practices

- References are verified before citation. Two fabricated entries were caught in
  an AI-generated bibliography early in the project; the 2.40 error is the same
  failure applied to a number instead of a citation.
- Verify every versioned file, not the published artefact.
- Correct a result and the script that generates it in the same commit.
- Work outside OneDrive. Synchronisation corrupts `.git`.
- PowerShell reads `.ps1` as ANSI unless the file has a UTF-8 BOM. Any script
  containing en dashes, multiplication signs or accented characters must be
  saved with a BOM or it will fail to parse.
- `[IO.File]` uses .NET's working directory, not PowerShell's. Always pass
  `(Resolve-Path $p).Path`.
- Repository files use CRLF. Multi-line string matching must normalise line
  endings or it silently fails.
