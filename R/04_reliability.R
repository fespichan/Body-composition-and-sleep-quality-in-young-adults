# 04 — Reliability on observed components (Supplementary Table S1).
# Primary 6-component estimate plus the 7-component with/without
# duration-adjustment sensitivity.
#
# Reads : data/raw_private/Raw_Data_Procedimiento.xlsx (canonical scored workbook)
# Writes: console only
#
# ⚠ DEPOSIT GAP: this step produces a PUBLISHED table but reads the PRIVATE
#   workbook, so it cannot be run from a public deposit as it stands.

suppressMessages({
  library(here)
  library(readxl)
  library(psych)
})
source(here("R", "lib_psqi_scoring.R"))
source(here("R", "lib_reliability.R"))

workbook <- here("data", "raw_private", "Raw_Data_Procedimiento.xlsx")
stopifnot(file.exists(workbook))
raw <- read_excel(workbook, sheet = "Raw Data", col_types = "text")
num <- function(x) suppressWarnings(as.numeric(x))

parsed <- data.frame(
  bedtime           = num(raw$P1_hora_parseada),
  wake_time         = num(raw$P3_hora_parseada),
  bedtime_ambiguous = raw$P1_ambiguo == "SI",
  sleep_hours       = num(raw$P4_horas_original),
  quality           = num(raw$P6_score),
  latency           = num(raw$P2_score),
  early_difficulty  = num(raw$P5a_score),
  disturbance_sum   = num(raw$P5b_score) + num(raw$P5c_score) + num(raw$P5d_score) +
    num(raw$P5e_score) + num(raw$P5f_score) + num(raw$P5g_score) +
    num(raw$P5h_score) + num(raw$P5i_score) + num(raw$P5j_score),
  medication        = num(raw$P7_score),
  somnolence        = num(raw$P8_score),
  enthusiasm        = num(raw$P9_score)
)

comp_adj <- score_psqi(parsed, adjust_duration = TRUE)
comp_raw <- score_psqi(parsed, adjust_duration = FALSE)

six <- c("C1", "C2", "C3", "C5", "C6", "C7")
seven <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7")

report <- function(label, comp, cols, target_alpha) {
  a <- cronbach_alpha(comp[cols])
  w <- mcdonald_omega(comp[cols])
  cat(sprintf("%-46s alpha=%.3f  omega=%.3f  n=%d   [target alpha~%.3f]\n",
              label, a$alpha, w$omega, a$n, target_alpha))
  invisible(a)
}

cat("Reliability on OBSERVED components\n")
cat("==================================\n")
report("6 components (excl. C4), full sample [PRIMARY]", comp_adj, six, 0.695)
report("7 components, complete cases, WITH adjustment",  comp_adj, seven, 0.711)
report("7 components, complete cases, WITHOUT adjustment", comp_raw, seven, 0.719)

# Confirm the 6-component estimate is invariant to the duration adjustment
# (the floor/ceiling change only C4, which is excluded here).
a6_adj <- cronbach_alpha(comp_adj[six])$alpha
a6_raw <- cronbach_alpha(comp_raw[six])$alpha
cat(sprintf("\n6-component alpha identical with/without adjustment: %s (%.5f vs %.5f)\n",
            isTRUE(all.equal(a6_adj, a6_raw)), a6_adj, a6_raw))
