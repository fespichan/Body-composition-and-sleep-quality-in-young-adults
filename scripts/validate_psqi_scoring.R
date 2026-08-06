# Provenance check: verify the R scoring layer (R/lib_psqi_scoring.R) reproduces
# the canonical scorer's components exactly. Reads the private workbook; not part
# of the public pipeline and not part of the numbered run order.

suppressMessages({
  library(here)
  library(readxl)
})
source(here("R", "lib_psqi_scoring.R"))

workbook <- here("data", "raw_private", "Raw_Data_Procedimiento.xlsx")
stopifnot(file.exists(workbook))
raw <- read_excel(workbook, sheet = "Raw Data", col_types = "text")

num <- function(x) suppressWarnings(as.numeric(x))

disturbance_sum <- num(raw$P5b_score) + num(raw$P5c_score) + num(raw$P5d_score) +
  num(raw$P5e_score) + num(raw$P5f_score) + num(raw$P5g_score) +
  num(raw$P5h_score) + num(raw$P5i_score) + num(raw$P5j_score)

parsed <- data.frame(
  bedtime           = num(raw$P1_hora_parseada),
  wake_time         = num(raw$P3_hora_parseada),
  bedtime_ambiguous = raw$P1_ambiguo == "SI",
  sleep_hours       = num(raw$P4_horas_original),
  quality           = num(raw$P6_score),
  latency           = num(raw$P2_score),
  early_difficulty  = num(raw$P5a_score),
  disturbance_sum   = disturbance_sum,
  medication        = num(raw$P7_score),
  somnolence        = num(raw$P8_score),
  enthusiasm        = num(raw$P9_score)
)

scored <- score_psqi(parsed)

# Reference: the canonical blue columns.
ref <- data.frame(
  C1 = num(raw$C1), C2 = num(raw$C2), C3 = num(raw$C3), C4 = num(raw$C4),
  C5 = num(raw$C5), C6 = num(raw$C6), C7 = num(raw$C7),
  PSQI_Total = num(raw$PSQI_Total)
)
ref_class <- ifelse(grepl("[Bb]uen", raw$Clasificacion), "good",
                    ifelse(grepl("[Mm]al", raw$Clasificacion), "poor", NA))

equal_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)

cat(sprintf("Rows: %d\n\n", nrow(raw)))
cat("Column agreement vs canonical blue columns:\n")
all_ok <- TRUE
for (col in c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "PSQI_Total")) {
  ok <- equal_na(scored[[col]], ref[[col]])
  cat(sprintf("  %-11s %3d / %3d match\n", col, sum(ok), nrow(raw)))
  if (any(!ok)) {
    all_ok <- FALSE
    bad <- which(!ok)
    print(data.frame(ide = raw$ide[bad], P1 = raw$P1_raw[bad],
                     ours = scored[[col]][bad], canonical = ref[[col]][bad]))
  }
}
ok_class <- equal_na(scored$Classification, ref_class)
cat(sprintf("  %-11s %3d / %3d match\n", "Classific.", sum(ok_class), nrow(raw)))

# Also confirm the duration adjustment reproduces the canonical adjusted hours.
adj_ok <- equal_na(adjust_sleep_hours(num(raw$P4_horas_original)),
                   num(raw$P4_horas_ajustada))
cat(sprintf("\nDuration adjustment vs P4_horas_ajustada: %d / %d match\n",
            sum(adj_ok), nrow(raw)))

cat(if (all_ok && all(ok_class)) "\nVALIDATION PASSED: components reproduce canonical scorer.\n"
    else "\nVALIDATION INCOMPLETE: see disagreements above.\n")
