#!/usr/bin/env bash

###############################################################################
# check_pending_rerun.sh
#
# Pre-publication check that the pages show the figures of the canonical run.
# Checks, each of which can fail the script:
#
#   1. Markers. No page (.qmd or .md) may contain "PENDING RERUN", a figure not yet
#      filled in from the canonical run.
#   2. Superseded values. No page may contain a figure listed in
#      scripts/dev/superseded_values.tsv, i.e. one that an earlier run put on
#      the page and a later run replaced. Removing the marker is not enough;
#      the old number must be gone too.
#   2b. Forbidden characters. No page (.qmd or .md) may contain an em dash
#      (U+2014), which the house style excludes.
#   2c. Reachable pages. Every sections/**/*.qmd page must be listed in the
#      sidebar of _quarto.yml, or excluded from rendering there ("!path").
#      A page that renders but is absent from the sidebar is published and
#      unreachable.
#   3. Inline figures. Pages that read their figures from the results files
#      (`r ...` with R/facts.R) are knitted into a temporary directory. A
#      missing results file, a renamed key or an expression left unevaluated
#      fails the check. Needs R with knitr; skipped with a warning otherwise.
#
# USAGE
#   bash scripts/dev/check_pending_rerun.sh            # all checks
#   bash scripts/dev/check_pending_rerun.sh --no-knit  # all but check 3
#
# EXIT STATUS
#   0  all checks run passed; the pages are safe to render and publish
#   1  at least one check failed; each problem is listed
#
# Run it before `quarto render` / `quarto publish`. run_all.sh calls it at the
# start as a warning only.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

KNIT=1
for a in "$@"; do
  case "$a" in
    --no-knit) KNIT=0 ;;
    -h|--help) sed -n '3,30p' "$0"; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

FAIL=0
LIST="$SCRIPT_DIR/superseded_values.tsv"

# Every .qmd and .md file in the repository, excluding build and cache
# directories: Quarto renders the .qmd pages and some .md files (CONTRIBUTING,
# docs/), and the READMEs are what GitHub shows.
mapfile -t PAGES < <(find . \( -name '*.qmd' -o -name '*.md' \) \
  -not -path './.git/*' -not -path './_freeze/*' -not -path './_site/*' \
  -not -path './.quarto/*' -not -path './tools/*' | sed 's|^\./||' | sort)

# ------------------------------------------------------------ 1. markers ----
MARKER="PENDING RERUN"
HITS="$(grep -nF "$MARKER" "${PAGES[@]}" || true)"
if [ -n "$HITS" ]; then
  N="$(printf '%s\n' "$HITS" | wc -l | tr -d ' ')"
  N_FILES="$(printf '%s\n' "$HITS" | cut -d: -f1 | sort -u | wc -l | tr -d ' ')"
  printf '%s\n' "$HITS"
  printf 'FAILED: %s line(s) in %s .qmd/.md file(s) still contain "%s".\n\n' \
    "$N" "$N_FILES" "$MARKER" >&2
  FAIL=1
else
  printf 'OK: no "%s" markers in .qmd or .md files.\n' "$MARKER"
fi

# --------------------------------------------------- 2. superseded values ----
if [ ! -f "$LIST" ]; then
  echo "FAILED: $LIST not found." >&2
  FAIL=1
else
  N_OLD=0; N_ENTRIES=0
  while IFS=$'\t' read -r glob mode value note; do
    case "$glob" in ''|'#'*|page) continue ;; esac
    N_ENTRIES=$((N_ENTRIES + 1))
    # shellcheck disable=SC2206  # the glob is meant to expand
    files=( $glob )
    [ -e "${files[0]}" ] || continue
    case "$mode" in
      F) hits="$(grep -nF -- "$value" "${files[@]}" /dev/null || true)" ;;
      E) hits="$(grep -nE -- "$value" "${files[@]}" /dev/null || true)" ;;
      *) echo "superseded_values.tsv: bad mode '$mode' for '$value'" >&2; FAIL=1; continue ;;
    esac
    if [ -n "$hits" ]; then
      while IFS= read -r h; do
        printf '%s\n    superseded: %s  (%s)\n' "$h" "$value" "$note"
        N_OLD=$((N_OLD + 1))
      done <<< "$hits"
    fi
  done < "$LIST"
  if [ "$N_OLD" -gt 0 ]; then
    printf 'FAILED: %s superseded figure(s) still on the pages.\n\n' "$N_OLD" >&2
    FAIL=1
  else
    printf 'OK: none of the %s superseded figures in %s is on a page.\n' \
      "$N_ENTRIES" "${LIST#"$PROJECT_ROOT"/}"
  fi
fi

# ------------------------------------------------ 2b. forbidden characters ----
EMDASH=$'\u2014'
HITS="$(grep -nF -- "$EMDASH" "${PAGES[@]}" || true)"
if [ -n "$HITS" ]; then
  printf '%s\n' "$HITS"
  printf 'FAILED: %s line(s) contain an em dash; the house style excludes it.\n\n' \
    "$(printf '%s\n' "$HITS" | wc -l | tr -d ' ')" >&2
  FAIL=1
else
  printf 'OK: no em dashes in .qmd or .md files.\n'
fi

# ---------------------------------------------------- 2c. reachable pages ----
mapfile -t SECTION_PAGES < <(printf '%s\n' "${PAGES[@]}" | grep '^sections/.*\.qmd$' || true)
N_MISSING=0
for p in "${SECTION_PAGES[@]}"; do
  if ! grep -qE "^[[:space:]]*-[[:space:]]*(href:[[:space:]]*)?\"?$p\"?[[:space:]]*$" _quarto.yml &&
     ! grep -qF -- "\"!$p\"" _quarto.yml; then
    printf '%s\n    not in the sidebar of _quarto.yml and not excluded from rendering\n' "$p"
    N_MISSING=$((N_MISSING + 1))
  fi
done
if [ "$N_MISSING" -gt 0 ]; then
  printf 'FAILED: %s section page(s) would render without a sidebar entry.\n\n' "$N_MISSING" >&2
  FAIL=1
else
  printf 'OK: all %s section pages are in the sidebar of _quarto.yml.\n' "${#SECTION_PAGES[@]}"
fi

# ------------------------------------------------------ 3. inline figures ----
# Only .qmd pages are knitted; a .md file never evaluates inline R.
mapfile -t INLINE < <(printf '%s\n' "${PAGES[@]}" | grep '\.qmd$' | tr '\n' '\0' | xargs -0 grep -lE '`r [^`]+`' || true)
if [ "$KNIT" -eq 0 ]; then
  printf 'SKIPPED: knitting %s page(s) with inline figures (--no-knit).\n' "${#INLINE[@]}"
elif [ "${#INLINE[@]}" -eq 0 ]; then
  printf 'OK: no page uses inline figures.\n'
elif ! command -v Rscript >/dev/null 2>&1 || ! Rscript -e 'quit(status = !requireNamespace("knitr", quietly = TRUE))' >/dev/null 2>&1; then
  printf 'WARNING: R with knitr not available; %s page(s) with inline figures not checked.\n' "${#INLINE[@]}" >&2
else
  TMP="$(mktemp -d)"
  trap 'rm -rf "$TMP"' EXIT
  if QUARTO_PROJECT_DIR="$PROJECT_ROOT" Rscript - "$TMP" "${INLINE[@]}" <<'RSCRIPT'
args <- commandArgs(TRUE); out <- args[1]; pages <- args[-1]
suppressMessages(library(knitr))
knitr::pat_md(); opts_knit$set(progress = FALSE)
bad <- 0L
# Inspect the values the inline expressions produce, not the page text, which
# may legitimately contain words such as NaN.
na_vals <- character()
default_inline <- knit_hooks$get("inline")
knit_hooks$set(inline = function(x) {
  if (length(x) == 0L || any(is.na(x)) || any(trimws(as.character(x)) %in% c("NA", "NaN", "Inf", "-Inf", "")))
    na_vals <<- c(na_vals, paste(format(x), collapse = " "))
  default_inline(x)
})
for (p in pages) {
  o <- file.path(out, paste0(gsub("[/.]", "_", p), ".md"))
  na_vals <- character()
  owd <- setwd(dirname(p))
  msg <- tryCatch({ suppressMessages(knit(basename(p), output = o, quiet = TRUE, envir = new.env())); NULL },
                  error = function(e) conditionMessage(e))
  setwd(owd)
  if (!is.null(msg)) {
    cat(sprintf("%s\n    knit failed: %s\n", p, msg)); bad <- bad + 1L; next
  }
  x <- readLines(o, warn = FALSE)
  left <- grep("`r [^`]+`", x)
  for (i in left) cat(sprintf("%s\n    unevaluated inline code in output line %d\n", p, i))
  for (v in na_vals) cat(sprintf("%s\n    an inline figure evaluated to NA/NaN/Inf or empty: '%s'\n", p, v))
  bad <- bad + length(left) + length(na_vals)
}
quit(status = if (bad > 0L) 1L else 0L)
RSCRIPT
  then
    printf 'OK: %s page(s) with inline figures knit cleanly from the results files.\n' "${#INLINE[@]}"
  else
    printf 'FAILED: inline figures could not be computed from the results files.\n\n' >&2
    FAIL=1
  fi
fi

if [ "$FAIL" -ne 0 ]; then
  printf 'Fix the pages from the canonical run before publishing.\n' >&2
  exit 1
fi
printf 'All checks passed.\n'
