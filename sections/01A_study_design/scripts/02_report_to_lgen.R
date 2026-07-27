#!/usr/bin/env Rscript
################################################################################
# Section 1A: Genotyping technologies — Step 02: Final Report to PLINK, using
#                                                the array manifest
#
# PURPOSE:
#   Convert raw calling output into files PLINK can read, using the array
#   manifest as the source of coordinates, alleles and strand.
#
#   The route is:
#       Final Report (.csv) + manifest (.csv)
#           -> .lgen + .map + .fam   (this script)
#           -> .ped/.map -> .bed/.bim/.fam   (Step 03)
#
# WHY THE MANIFEST IS NOT OPTIONAL:
#   A Final Report names variants and reports their alleles. It does not say
#   where they are, which genome build they refer to, or which strand the
#   design used. All three come from the manifest.
#
#   It is tempting to shortcut this when variant names look like coordinates,
#   for instance "1:103380393". Do not. Measured on the example data of this
#   tutorial: of 5,745 variants whose names encode a position, only 776 agree
#   with the manifest coordinate. The other 4,969 differ by roughly half a
#   megabase, because the names carry GRCh37 positions while the A2 manifest is
#   GRCh38, and 75 variants disagree even on the chromosome. Coordinates taken
#   from names would be wrong for 86% of variants, and would look entirely
#   plausible. Nothing downstream would complain.
#
# WHY STRAND MATTERS HERE:
#   The manifest reports the design alleles in its SNP column, as "[A/G]", on
#   the design strand. The RefStrand column says whether that strand is the
#   plus strand of the reference. In the example data 3,229 of 7,508 variants,
#   43%, are on the minus strand, so their alleles must be complemented before
#   they can be compared with genotypes reported on the plus strand. Skipping
#   this gives wrong reference alleles for nearly half the array, silently.
#
# THE GENCALL FILTER:
#   Calls below the threshold are written as missing rather than deleted, so
#   the loss stays visible as missingness for the quality control of Section 1B
#   to find. The threshold is a parameter, its effect is counted, and the count
#   goes into the log. An unfiltered dataset carries unreliable calls that
#   nobody can quantify afterwards.
#
# INPUT:
#   - one or more Final Report .csv files
#   - the array manifest .csv for the product named in the report header
#
# OUTPUT, per sample:
#   - results/raw/lgen/<sample>.lgen / .map / .fam
#   - results/raw/pdac_demo_01A_manifest_check.txt
#   - results/raw/pdac_demo_01A_conversion_log.tsv
#
# Usage:
#   Rscript 02_report_to_lgen.R [report_dir_or_file] [manifest.csv] [gc_min]
################################################################################

args      <- commandArgs(trailingOnly = TRUE)
report_in <- if (length(args) >= 1) args[1] else "example_data/example_final_report.csv"
manifest  <- if (length(args) >= 2) args[2] else "example_data/manifest.csv"
gc_min    <- if (length(args) >= 3) as.numeric(args[3]) else 0.15

out_dir <- "results/raw"
lgen_dir <- file.path(out_dir, "lgen")
dir.create(lgen_dir, showWarnings = FALSE, recursive = TRUE)

check_file <- file.path(out_dir, "pdac_demo_01A_manifest_check.txt")
log_file   <- file.path(out_dir, "pdac_demo_01A_conversion_log.tsv")

con <- file(check_file, open = "wt")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = con) }

# --- helpers -----------------------------------------------------------------

read_report_header <- function(f) {
  h <- readLines(f, n = 40)
  data_line <- grep("^SNP Name", h)[1]
  if (is.na(data_line)) stop("Could not find the data header in ", f)
  get <- function(key) {
    i <- grep(paste0("^", key), h)[1]
    if (is.na(i)) return(NA_character_)
    trimws(sub("^[^;,]*[;,]+", "", h[i]))
  }
  list(skip = data_line - 1,
       sep  = if (grepl(";", h[data_line])) ";" else ",",
       content  = get("Content"),
       num_snps = suppressWarnings(as.numeric(get("Num SNPs"))),
       num_samples = suppressWarnings(as.numeric(get("Num Samples"))))
}

read_manifest_header <- function(f) {
  h <- readLines(f, n = 60)
  get <- function(key) {
    i <- grep(paste0("^", key), h)[1]
    if (is.na(i)) return(NA_character_)
    trimws(sub("^[^,]*,", "", h[i]))
  }
  # The column header is the first line that names both Name and MapInfo. This
  # works for a full manufacturer manifest, which carries a [Heading] block, and
  # for a plain extract that starts directly with the column names.
  hdr <- grep("(^|,)\"?Name\"?,", h)
  hdr <- hdr[sapply(hdr, function(i) grepl("MapInfo", h[i]))][1]
  if (is.na(hdr)) stop("Could not find the column header line in ", f)
  list(hdr_line = hdr,
       colnames = gsub('^"|"$', "", trimws(strsplit(h[hdr], ",")[[1]])),
       descriptor = get("Descriptor File Name"),
       loci_count = suppressWarnings(as.numeric(get("Loci Count"))),
       date = get("Date Manufactured"))
}

complement <- c(A = "T", T = "A", C = "G", G = "T"[0])  # placeholder, set below
complement <- c(A = "T", T = "A", C = "G", G = "C")

# --- 1. the manifest ---------------------------------------------------------

if (!file.exists(manifest)) {
  w("MANIFEST NOT FOUND: ", manifest)
  w("")
  w("This script cannot run without the array manifest for the product named in")
  w("the Final Report header. Manifests are distributed free by the array")
  w("manufacturer; download the CSV manifest for exactly the product and version")
  w("your samples were genotyped on, and pass its path as the second argument:")
  w("")
  w("  Rscript 02_report_to_lgen.R <report> <manifest.csv> [gc_min]")
  w("")
  w("Do not substitute a manifest from a different product or version. See the")
  w("notes at the top of this script for what happens if you do.")
  close(con)
  quit(status = 1)
}

mh <- read_manifest_header(manifest)
# quote = "" because manufacturer manifests contain unbalanced quotation marks
# in sequence fields, which silently truncate a default read.
man <- read.csv(manifest, header = FALSE, skip = mh$hdr_line, sep = ",",
                quote = "", comment.char = "", stringsAsFactors = FALSE,
                col.names = mh$colnames, check.names = FALSE)
for (j in seq_along(man)) if (is.character(man[[j]])) man[[j]] <- gsub('^"|"$', "", man[[j]])

needed <- c("Name", "Chr", "MapInfo", "SNP", "RefStrand")
missing_cols <- setdiff(needed, names(man))
if (length(missing_cols)) stop("Manifest is missing columns: ", paste(missing_cols, collapse = ", "))

# --- 2. the report -----------------------------------------------------------

report_files <- if (dir.exists(report_in)) {
  list.files(report_in, pattern = "\\.csv$", full.names = TRUE)
} else report_in
# A directory may hold the manifest as well. A Final Report is recognised by its
# data header; anything without one is not a report.
report_files <- setdiff(normalizePath(report_files), normalizePath(manifest))
is_report <- vapply(report_files, function(f) {
  any(grepl("^SNP Name", readLines(f, n = 40)))
}, logical(1))
report_files <- report_files[is_report]
if (!length(report_files)) stop("No Final Report files found at ", report_in)

rh <- read_report_header(report_files[1])
rep1 <- read.csv(report_files[1], header = TRUE, skip = rh$skip, sep = rh$sep,
                 check.names = FALSE, stringsAsFactors = FALSE)
names(rep1) <- make.names(names(rep1))
snp_col <- grep("^SNP.Name$", names(rep1), value = TRUE)[1]
report_snps <- unique(rep1[[snp_col]])

# --- 3. the four checks ------------------------------------------------------

w("Manifest check")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("")
w("Report:   ", report_files[1])
w("Manifest: ", manifest)
w("")
w("1. PRODUCT")
w("   report header 'Content':            ", rh$content)
w("   manifest 'Descriptor File Name':    ", mh$descriptor)
prod_ok <- !is.na(rh$content) && !is.na(mh$descriptor) &&
  identical(sub("\\.bpm$", "", basename(rh$content)), sub("\\.bpm$", "", mh$descriptor))
w("   identical: ", if (prod_ok) "yes" else "NO")
if (!prod_ok) {
  w("   The manifest describes a different product or version from the one that")
  w("   produced these calls. Variants present only in the other product will")
  w("   not be mapped. Continue only if you understand the difference.")
}
w("")

w("2. SIZE")
w("   variants in this report:  ", length(report_snps))
w("   loci in the manifest:     ", mh$loci_count)
w("")

w("3. NAME MATCHING")
matched <- sum(report_snps %in% man$Name)
pct <- 100 * matched / length(report_snps)
w("   matched as written:       ", matched, sprintf(" (%.1f%%)", pct))
alt <- sum(sub("^GSA-", "", man$Name) %in% report_snps)
w("   matched after stripping a 'GSA-' prefix from manifest names: ", alt)
if (alt > matched) {
  w("   The manifest prefixes its identifiers. Names are being stripped.")
  man$Name <- sub("^GSA-", "", man$Name)
  matched <- alt
  pct <- 100 * matched / length(report_snps)
}
if (pct < 90) {
  w("   WARNING: fewer than 90% of variants matched. Check that the manifest")
  w("   corresponds to the product in the report header before trusting the output.")
}
w("")

w("4. GENOME BUILD")
if ("GenomeBuild" %in% names(man)) {
  builds <- table(man$GenomeBuild[man$Name %in% report_snps])
  w("   manifest GenomeBuild: ", paste(names(builds), builds, sep = " (n=", collapse = "), "), ")")
} else {
  w("   the manifest has no GenomeBuild column; establish the build from its")
  w("   documentation and record it with the data")
}
w("")

# --- 4. a diagnostic that is worth running once ------------------------------
# Where variant names look like coordinates, do they agree with the manifest?
w("5. DIAGNOSTIC: names that look like coordinates")
sub_man <- man[man$Name %in% report_snps, ]
pos_like <- grepl("^([0-9]+|X|Y|XY|MT):[0-9]+", sub_man$Name)
if (any(pos_like)) {
  nm  <- sub_man[pos_like, ]
  parts <- do.call(rbind, strsplit(sub("^([^:]+):([0-9]+).*", "\\1|\\2", nm$Name), "|", fixed = TRUE))
  name_chr <- parts[, 1]; name_pos <- as.numeric(parts[, 2])
  agree <- sum(name_pos == as.numeric(nm$MapInfo), na.rm = TRUE)
  chrmis <- sum(name_chr != as.character(nm$Chr), na.rm = TRUE)
  w("   variants with a positional name: ", nrow(nm))
  w("   whose name position equals the manifest coordinate: ", agree)
  w("   whose name position differs:                        ", nrow(nm) - agree)
  w("   whose name chromosome differs:                      ", chrmis)
  w("   Positional names are not coordinates. They typically carry an older")
  w("   build. Always take positions from the manifest.")
}
w("")
close(con)

# --- 5. build the reference table --------------------------------------------
# Design alleles come from the SNP column as "[A/G]" on the design strand;
# RefStrand says whether that strand is the plus strand of the reference.

ref <- sub_man[, c("Name", "Chr", "MapInfo", "SNP", "RefStrand")]
ref$A1raw <- sub("^\\[([A-Z-]+)/.*", "\\1", ref$SNP)
ref$A2raw <- sub(".*/([A-Z-]+)\\]$", "\\1", ref$SNP)
flip <- ref$RefStrand == "-"
ref$A1 <- ifelse(flip, complement[ref$A1raw], ref$A1raw)
ref$A2 <- ifelse(flip, complement[ref$A2raw], ref$A2raw)

# Chromosome codes PLINK expects.
ref$Chr <- as.character(ref$Chr)
ref$Chr[ref$Chr == "X"]  <- "23"
ref$Chr[ref$Chr == "Y"]  <- "24"
ref$Chr[ref$Chr == "XY"] <- "25"
ref$Chr[ref$Chr == "MT"] <- "26"
ref$Chr[ref$Chr == "0" | ref$Chr == ""] <- "0"

write.csv(ref[, c("Name", "Chr", "MapInfo", "A1", "A2", "RefStrand")],
          file.path(out_dir, "pdac_demo_01A_reference_alleles.csv"), row.names = FALSE)

cat("\nReference table built for", nrow(ref), "variants;",
    sum(flip), "had their alleles complemented (RefStrand '-').\n\n")

# --- 6. convert each report --------------------------------------------------

log_rows <- list()

for (f in report_files) {
  h <- read_report_header(f)
  d <- read.csv(f, header = TRUE, skip = h$skip, sep = h$sep,
                check.names = FALSE, stringsAsFactors = FALSE)
  names(d) <- make.names(names(d))

  c_snp <- grep("^SNP.Name$", names(d), value = TRUE)[1]
  c_id  <- grep("^Sample.ID$", names(d), value = TRUE)[1]
  c_gc  <- grep("^GC.Score$",  names(d), value = TRUE)[1]
  c_a1  <- grep("^Allele1",    names(d), value = TRUE)[1]
  c_a2  <- grep("^Allele2",    names(d), value = TRUE)[1]

  if (!grepl("Plus", c_a1)) {
    cat("NOTE: allele columns are '", c_a1, "', not plus-strand. TOP/BOT and A/B\n",
        "exports use different conventions and are not handled here. Request a\n",
        "plus-strand export, or convert using the manifest IlmnStrand column.\n", sep = "")
  }

  d$..gc <- suppressWarnings(as.numeric(d[[c_gc]]))
  d$..a1 <- d[[c_a1]]
  d$..a2 <- d[[c_a2]]

  n_total  <- nrow(d)
  acgt     <- c("A", "C", "G", "T")

  # Three distinct things hide in the allele columns, and collapsing them is a
  # common silent error:
  #   "-"       a genuine no-call
  #   "I"/"D"   an insertion or deletion call; the array carries indel probes
  #   A/C/G/T   an ordinary SNP call
  no_call  <- d$..a1 == "-" | d$..a2 == "-"
  indel    <- (d$..a1 %in% c("I", "D")) | (d$..a2 %in% c("I", "D"))
  bad_gc   <- is.na(d$..gc) | d$..gc < gc_min

  # Decision, declared rather than implied: indels are dropped. PLINK's biallelic
  # A/C/G/T model and the manifest's SNP column do not represent them cleanly, and
  # most GWAS pipelines analyse SNPs only. If your study needs indels, handle them
  # separately from a VCF rather than forcing them through this route. Either way
  # the count is reported, so the loss is never invisible.
  bad_call <- !(d$..a1 %in% acgt) | !(d$..a2 %in% acgt)
  set_miss <- bad_gc | bad_call

  # Missing calls are simply omitted from the .lgen; PLINK treats absent
  # sample-variant pairs as missing, which keeps them visible in Section 1B.
  keep <- !set_miss & d[[c_snp]] %in% ref$Name
  n_unmapped <- sum(!(d[[c_snp]] %in% ref$Name))

  for (s in unique(d[[c_id]])) {
    sel <- keep & d[[c_id]] == s
    sid <- gsub("[/ ]", "_", s)

    lgen <- data.frame(FID = sid, IID = sid,
                       SNP = d[[c_snp]][sel],
                       A1 = d$..a1[sel], A2 = d$..a2[sel])
    write.table(lgen, file.path(lgen_dir, paste0(sid, ".lgen")),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

    m <- ref[match(lgen$SNP, ref$Name), c("Chr", "Name", "MapInfo")]
    m$morgan <- 0
    write.table(m[, c("Chr", "Name", "morgan", "MapInfo")],
                file.path(lgen_dir, paste0(sid, ".map")),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

    write.table(data.frame(sid, sid, 0, 0, 0, -9),
                file.path(lgen_dir, paste0(sid, ".fam")),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

    inS <- d[[c_id]] == s
    log_rows[[length(log_rows) + 1]] <- data.frame(
      sample     = sid,
      calls      = sum(inS),
      no_call    = sum(no_call & inS),
      indel      = sum(indel & inS),
      failed_gc  = sum(bad_gc & inS),
      unmapped   = sum(!(d[[c_snp]] %in% ref$Name) & inS),
      written    = sum(sel))
  }
  cat("converted ", f, ": ", n_total, " calls | no-call ", sum(no_call),
      " | indel ", sum(indel), " | GC<", gc_min, " ", sum(bad_gc),
      " | unmapped ", n_unmapped, "\n", sep = "")
}

log_df <- do.call(rbind, log_rows)
write.table(log_df, log_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nPer-sample conversion summary:\n")
print(log_df, row.names = FALSE)
cat("\nWritten:\n  ", check_file, "\n  ", log_file, "\n  ", lgen_dir, "/*.lgen|.map|.fam\n", sep = "")
cat("\n=== NEXT STEP ===\n\n")
cat("  bash scripts/01A_study_design/03_lgen_to_binary.sh\n\n")
