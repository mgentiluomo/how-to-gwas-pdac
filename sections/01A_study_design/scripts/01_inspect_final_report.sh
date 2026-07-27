#!/usr/bin/env bash

################################################################################
# Section 1A: Genotyping technologies — Step 01: what comes off the machine
#
# PURPOSE:
#   Look at real raw genotyping output before converting it. Almost every
#   downstream problem that is expensive to fix, wrong strand, wrong build,
#   wrong allele coding, silently low-quality calls, is visible here and
#   invisible later.
#
# WHAT A FINAL REPORT IS:
#   An Illumina array is read by the scanner, the intensities are clustered by
#   the calling software (GenomeStudio), and the result is exported as a "Final
#   Report": a text file listing, for every variant and every sample, the two
#   called alleles and a confidence score.
#
#   Affymetrix/Thermo Fisher platforms produce the equivalent from .CEL files
#   via their own calling pipeline, with different column names and conventions.
#   The principles below are the same.
#
# THE GENCALL SCORE:
#   Every call carries a GenCall (GC) score between 0 and 1, measuring how
#   cleanly the sample fell into its genotype cluster. Illumina's guidance is
#   that calls with GC below about 0.15 should be treated as failures. Note what
#   this means: a low-GC call is not missing data in the file. It is a genotype
#   that looks like any other unless you check.
#
# INPUT:
#   - example_data/example_final_report.csv
#     A genuine Illumina Final Report from a GSA-24v3 array, the same array
#     family the demonstration dataset is built on. Excerpt of 8,000 variants
#     for one sample; the sample identifier has been replaced with DEMO001.
#
# OUTPUT:
#   - results/pdac_demo_01A_gc_summary.txt
#   - results/pdac_demo_01A_gc_distribution.png
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -d "$SCRIPT_DIR/../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
elif [ -d "$SCRIPT_DIR/../../../scripts/dev" ]; then
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
else
  PROJECT_ROOT="$(pwd)"
fi
cd "$PROJECT_ROOT"

REPORT="${1:-example_data/example_final_report.csv}"
OUT_DIR="${2:-results/raw}"
GC_MIN="${3:-0.15}"

mkdir -p "$OUT_DIR"

echo ""
echo "=== The header ==="
echo ""
sed -n '1,10p' "$REPORT"

echo ""
echo "The header tells you three things you must record now, because they are"
echo "hard to recover later:"
echo "  1. which array was used (the .bpm manifest name on the Content line);"
echo "  2. which version of the calling software produced the calls;"
echo "  3. how many variants and samples the file claims to contain."
echo ""

echo "=== The first rows of data ==="
echo ""
sed -n '11,16p' "$REPORT"

echo ""
echo "Note the column names. 'Allele1 - Plus' means the alleles are reported on"
echo "the plus strand of the reference genome. Exports also exist in TOP/BOT and"
echo "in A/B coding, which are NOT interchangeable. Converting a TOP/BOT export"
echo "as though it were plus-strand silently flips a large fraction of variants."
echo "If the column header does not say which convention is used, stop and ask"
echo "the genotyping facility before doing anything else."
echo ""

# --- GenCall score distribution ----------------------------------------------
Rscript --vanilla -e "
report  <- '${REPORT}'
out_dir <- '${OUT_DIR}'
gc_min  <- ${GC_MIN}

# Find where the data block starts: the line beginning with 'SNP Name'.
lines <- readLines(report, n = 50)
skip  <- grep('^SNP Name', lines)[1]

d <- read.table(report, sep = ';', header = TRUE, skip = skip - 1,
                check.names = FALSE, stringsAsFactors = FALSE)
names(d) <- make.names(names(d))
gc <- as.numeric(d[[grep('GC.Score', names(d))[1]]])
gc <- gc[is.finite(gc)]

png(file.path(out_dir, 'pdac_demo_01A_gc_distribution.png'),
    width = 1600, height = 1200, res = 200)
hist(gc, breaks = 50, col = '#4B6584', border = 'white',
     xlab = 'GenCall score', ylab = 'Number of calls',
     main = 'Confidence of the raw genotype calls')
abline(v = gc_min, col = '#D55E00', lty = 2, lwd = 2)
legend('topleft', bty = 'n', lty = 2, col = '#D55E00',
       legend = sprintf('GC = %.2f, Illumina failure threshold', gc_min))
dev.off()

con <- file(file.path(out_dir, 'pdac_demo_01A_gc_summary.txt'), open = 'wt')
w <- function(...) { cat(..., '\n', sep = ''); cat(..., '\n', sep = '', file = con) }
w('GenCall score summary')
w('File: ', report)
w('')
w('Calls in file:        ', length(gc))
w('Median GC:            ', sprintf('%.4f', median(gc)))
w('Mean GC:              ', sprintf('%.4f', mean(gc)))
w('Calls below ', gc_min, ':     ', sum(gc < gc_min),
  sprintf(' (%.2f%%)', 100 * mean(gc < gc_min)))
w('Calls below 0.50:     ', sum(gc < 0.50), sprintf(' (%.2f%%)', 100 * mean(gc < 0.50)))
w('Calls below 0.70:     ', sum(gc < 0.70), sprintf(' (%.2f%%)', 100 * mean(gc < 0.70)))
close(con)
"

echo ""
echo "Written:"
echo "  ${OUT_DIR}/pdac_demo_01A_gc_summary.txt"
echo "  ${OUT_DIR}/pdac_demo_01A_gc_distribution.png"
echo ""
echo "=== Why this file cannot be analysed as it stands ==="
echo ""
SIZE=$(du -h "$REPORT" | cut -f1)
echo "This excerpt holds 8,000 variants for ONE sample and already weighs ${SIZE}."
echo "The full report it came from carries 730,059 variants for 130 samples."
echo "At that scale a text export runs to gigabytes, and every analysis would"
echo "re-parse all of it. This is why the first real task of a GWAS is a format"
echo "conversion, not a statistical test."
echo ""
echo "=== NEXT STEP ==="
echo ""
echo "  bash scripts/01A_study_design/02_report_to_lgen.sh"
echo ""
