# Project state

Last updated 28 July 2026. This file is the single place where the current
state of the project is written down. It is deliberately blunt about what is
finished, what is half-done, and what is known to be wrong, because the
alternative is that someone opens the repository and cannot tell.

---

## The short version

The repository is in a **fork**. Two mutually inconsistent versions of the
demonstration data exist, and only one of them is published.

> **31 July 2026.** The figures in the right-hand column below are no longer
> predictions. They are measured from the canonical run recorded in
> `RUN_MANIFEST.txt`, produced by `scripts/dev/run_all.sh` at commit `63d3287`
> against data release `v0.3-data`, on bio-plink3. Where they differ from what
> was written here before, the run supersedes the note: the earlier figures came
> from an unrecorded run made before three script defects were found. Lambda for
> the analysis set is **1.016**.

| | Published site, manuscript and repository | New dataset, not yet published |
|---|---|---|
| Samples after QC | 1,217 | **1,215** |
| Variants after QC | 414,695 | **414,662**, 414,655 after within-ancestry HWE |
| European analysis set | 224 cases / 387 controls | **219 / 392**, effective 562 |
| Causal variants | one, at MAF 0.086 | two; risk allele is the **reference** in both, at frequency 0.61 and 0.49 |
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
selection bias is small. The power section reported 9.2% power at the causal
variant; measured at the true effect it is **99.1%** in the analysis set and
100% across the strata. The claim that the locus was
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

- `9:133273682:A:T`, *ABO* 9q34.2, risk allele `A` (the reference, and the
  major allele) at frequency 0.61, generative OR **2.60**. Detected,
  *P* = 4.35e-9, estimated OR 0.428. Measured power 99.1% in the analysis set;
  2.30 was tried first, gave 91%, and the first draw fell in the missing 9%. The
  design parameter was raised, not the seed, because selecting a seed on
  significance is conditioning on significance.
- `5:1286401:C:A`, TERT/CLPTM1L 5p15.33, risk allele `C` (the reference) at
  frequency 0.49, generative OR **1.50**. Not
  detected in Europeans (*P* = 0.156); reaches only *P* = 0.021 across the three
  strata, with measured power 2.0% and 28.6% respectively, so it is **not**
  recovered; note the East Asian stratum alone reaches *P* = 0.007 while the
  combined result is weaker, because the African stratum points the other way,
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
| `v0.1-data` | genotypes, phenotype, covariates, survival | Superseded by `v0.3-data`. Left in place, not withdrawn |
| `v0.2-checkpoints` | checkpoint bundle, 8 files + manifest | Superseded by `v0.4-checkpoints`. Left in place, not withdrawn |
| `v0.3-data` | genotypes, new phenotype, covariates, survival, `truth.tsv`, `MANIFEST.tsv` | **Current.** Nine files, checksums verified after upload |
| `v0.4-checkpoints` | checkpoints and full summary statistics, 15 files | **Current.** Built from pipeline commit `84223a4`, verified against its manifest after upload |

Pages link to an explicit tag, never to `latest`, so a page and the data it was
written against move together.

---

## What remains, in order

> **31 July 2026, end of day.** The conversion is complete and published. The
> site at `mgentiluomo.github.io/how-to-gwas-pdac` serves the new dataset, the
> repository and the two current releases agree with it, and a sweep across
> every versioned file finds no figure from the previous run. The manuscript is
> the one part still describing the old one.

1. **The manuscript, version 12.** The only substantial work left. Six passages
   change in argument and not only in number, and the rewritten version of each
   already exists on the corresponding site page:

   - Abstract: one causal locus becomes two.
   - Results: every count.
   - The winner's curse passage: **delete and replace**. At 99.1% power there is
     no selection bias to show; the recovered odds ratio, 2.34, sits below the
     generative 2.60 rather than above it. The replacement argument, that
     selection acts on which variant is named rather than on the size of the
     effect, is written on the Section 4A page.
   - The power passage: the claim that the locus was detected despite being
     underpowered is **wrong**, not stale.
   - Fine mapping: the credible set does not collapse to one variant. It goes
     from 11 to 9, and what changes is the ranking of the causal variant, from
     third to first.
   - Heterogeneity: I² at this sample size points the wrong way. The causal
     variant has the highest I² in the region while three of its proxies have
     none.

   Also: add a Table 1 row on verifying that an analysis recovers the generative
   parameter before treating it as truth; regenerate Tables 2, 4 and 5; and fix
   the numbering, since there is a Table 1, 2, 4 and 5 and no Table 3.

2. **The four items that were already pending**, unchanged: the methodological
   refinements raised in review, the freeze and Zenodo deposit, the six sections
   not built, and the unstated limitation about cleanly separated continental
   groups.

## What was published, and when

| | |
|---|---|
| Site | `mgentiluomo.github.io/how-to-gwas-pdac`, published 31 July 2026 from `main` at `f191d1b` |
| Pages | Nine published stages, the linkage disequilibrium section among them |
| Data | `v0.3-data`, nine files, checksums verified after upload |
| Checkpoints | `v0.4-checkpoints`, fifteen files including the full summary statistics, verified against its manifest after upload |
| Canonical run | pipeline commit `84223a4` on bio-plink3, 31 steps, deterministic on re-run |

Re-running the whole pipeline on the same data and scripts changes only the
generation timestamps in ten summary files. Nothing else moves. That property
is worth preserving, and the timestamps are worth removing from those files so
that a rerun with no substantive change produces no diff at all.

## What was found during the conversion

Four defects that the previous, hand-run workflow had concealed:

- Three scripts assigned to `GROUPS`, a reserved bash array. Under `set -e` the
  assignment returns non-zero and the script dies silently. The lesson had been
  written into one script's comments months earlier and never applied to the
  others.
- The fine-mapping and figure scripts carried the previous causal variant as a
  hardcoded constant, so they reported that the causal variant was absent from a
  credible set that in fact contained it, and the manuscript figure marked a
  variant the dataset does not contain.
- The power script carried the previous odds ratio and frequency, which is where
  the 9.2% power figure came from.
- The risk allele at both causal loci is the **reference** allele, and the major
  allele at *ABO*. The simulation script's comment said otherwise. An
  association run therefore prints an odds ratio near 0.43 where the documented
  effect is 2.60.

None of these changed a result once corrected; all of them changed what the
resource *said* about its results. The runner found them because it executes
every step in order in a clean shell and stops at the first failure.

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
