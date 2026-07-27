#!/usr/bin/env Rscript

################################################################################
# Statistical power — Step 02: survival analysis obeys a different rule
#
# PURPOSE:
#   Show why the power of a time-to-event scan is governed by the number of
#   observed events and not by the number of participants, and compute it for
#   the demonstration data.
#
# THE RULE:
#   In a Cox model the partial likelihood compares, at each event time, the
#   genotype of the individual having the event with those still at risk. Only
#   event times contribute a comparison. The information about the genetic
#   coefficient is therefore
#
#       information  =  D  x  Var(G),        Var(G) = 2p(1 - p)
#
#   with D the number of observed events. Censored individuals enter the risk
#   sets but add nothing to the test. A cohort of several thousand people with
#   few events is, for this purpose, a small study.
#
#   The practical consequence for a rare and rapidly fatal cancer is
#   uncomfortable: a survival scan is powered by raising the event count,
#   through enriching for individuals at higher risk or extending follow up,
#   and not by recruiting more people who will be censored.
#
# INPUT:
#   - demo_data/survival.txt          FID IID TIME EVENT
#   - results/pca/pdac_demo_02_eur_keep.txt
#
# OUTPUT:
#   - results/power/pdac_demo_07_survival_power.txt
#   - results/power/pdac_demo_07_survival_power.png
################################################################################

out_dir <- "results/power"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 5e-8
crit  <- qchisq(1 - alpha, df = 1)

info_cox  <- function(D, maf) D * 2 * maf * (1 - maf)
power_cox <- function(D, maf, hr)
  pchisq(crit, 1, ncp = info_cox(D, maf) * log(hr)^2, lower.tail = FALSE)
min_hr <- function(D, maf, target = 0.80) {
  lam <- uniroot(function(l) pchisq(crit, 1, ncp = l, lower.tail = FALSE) - target,
                 c(1e-6, 500))$root
  exp(sqrt(lam / info_cox(D, maf)))
}

surv_file <- "demo_data/survival.txt"
keep_file <- "results/pca/pdac_demo_02_eur_keep.txt"

con <- file(file.path(out_dir, "pdac_demo_07_survival_power.txt"), open = "wt")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = con) }

w("Power for a time-to-event scan")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("")

if (!file.exists(surv_file)) {
  w("No survival file found at ", surv_file, ".")
  close(con); quit(status = 0)
}

sv <- read.table(surv_file, header = TRUE, comment.char = "", check.names = FALSE)
names(sv)[1] <- sub("^#", "", names(sv)[1])
ev_col <- grep("EVENT", names(sv), value = TRUE)[1]

if (file.exists(keep_file)) {
  kp <- read.table(keep_file, header = FALSE)
  sv <- sv[sv$IID %in% kp$V2, ]
  w("Restricted to the European analysis set.")
}

N <- nrow(sv)
D <- sum(sv[[ev_col]] == 1, na.rm = TRUE)
w("  individuals with survival information: ", N)
w("  observed events: ", D, sprintf("  (%.1f%%)", 100 * D / N))
w("  censored: ", N - D)
w("")
w("The number that governs power is ", D, ", not ", N, ".")
w("")

mafs <- c(0.30, 0.20, 0.10, 0.05)
w("SMALLEST HAZARD RATIO DETECTABLE AT 80% POWER, WITH ", D, " EVENTS")
for (p in mafs) w(sprintf("  MAF %.2f: HR %.2f", p, min_hr(D, p)))
w("")
w("POWER AT HAZARD RATIOS OF A REALISTIC SIZE")
for (hr in c(1.2, 1.5, 2.0))
  w(sprintf("  HR %.1f at MAF 0.30: %5.1f%%", hr, 100 * power_cox(D, 0.30, hr)))
w("")
w("HOW MANY EVENTS WOULD BE NEEDED")
for (hr in c(1.3, 1.5, 2.0)) {
  need <- uniroot(function(d) power_cox(d, 0.30, hr) - 0.80, c(10, 5e6))$root
  w(sprintf("  to detect HR %.1f at MAF 0.30 with 80%% power: %s events",
            hr, format(round(need), big.mark = ",")))
}
w("")
w("Compare those requirements with the event count actually available. This is")
w("the calculation that decides whether a survival scan is an analysis or an")
w("exploratory exercise, and it belongs in the protocol, not in the discussion.")
close(con)

png(file.path(out_dir, "pdac_demo_07_survival_power.png"),
    width = 1600, height = 1200, res = 200)
par(mar = c(4.4, 4.6, 2.5, 1), mgp = c(2.6, 0.7, 0))
Ds <- seq(50, 5000, by = 25)
plot(Ds, sapply(Ds, function(d) min_hr(d, 0.30)), type = "l", lwd = 2,
     col = "#0072B2", log = "y", ylim = c(1.1, 8),
     xlab = "Number of observed events", ylab = "Smallest hazard ratio detectable at 80% power",
     main = "Survival power is bought with events, not participants")
lines(Ds, sapply(Ds, function(d) min_hr(d, 0.10)), lwd = 2, col = "#D55E00")
abline(v = D, lty = 2, col = "grey40")
text(D, 7, sprintf("this dataset: %d events", D), pos = 4, cex = 0.7)
rect(min(Ds), 1.1, max(Ds), 1.5, col = "#00000010", border = NA)
legend("topright", bty = "n", lwd = 2, cex = 0.8, col = c("#0072B2", "#D55E00"),
       legend = c("MAF 0.30", "MAF 0.10"))
dev.off()

cat("\nWritten to ", out_dir, "\n", sep = "")
