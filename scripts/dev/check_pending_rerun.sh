#!/usr/bin/env bash

###############################################################################
# check_pending_rerun.sh
#
# Pre-publication check: fails if any .qmd page still contains the marker
# "PENDING RERUN", i.e. a figure that has not yet been filled in from the
# canonical run.
#
# USAGE
#   bash scripts/dev/check_pending_rerun.sh
#
# EXIT STATUS
#   0  no marker found; the pages are safe to render and publish
#   1  at least one marker found; each occurrence is listed as file:line
#
# Run it before `quarto render` / `quarto publish`. run_all.sh calls it at the
# start as a warning only, because the markers are expected until the run that
# replaces them has finished.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

MARKER="PENDING RERUN"

HITS="$(grep -rnF "$MARKER" --include='*.qmd' \
  --exclude-dir=.git --exclude-dir=_freeze --exclude-dir=_site \
  --exclude-dir=.quarto . || true)"

if [ -n "$HITS" ]; then
  N="$(printf '%s\n' "$HITS" | wc -l | tr -d ' ')"
  N_FILES="$(printf '%s\n' "$HITS" | cut -d: -f1 | sort -u | wc -l | tr -d ' ')"
  printf '%s\n' "$HITS"
  printf '\nFAILED: %s line(s) in %s .qmd file(s) still contain "%s".\n' \
    "$N" "$N_FILES" "$MARKER" >&2
  printf 'Fill them in from the canonical run before publishing.\n' >&2
  exit 1
fi

printf 'OK: no "%s" markers in .qmd files.\n' "$MARKER"
