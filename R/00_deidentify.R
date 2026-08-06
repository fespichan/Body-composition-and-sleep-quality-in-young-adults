############################################################
# 00_deidentify.R — produce the public, de-identified source file
#
# Reads : data/raw_private/sleep_raw_original.csv   (PII: ide, nombre)
# Writes: data/raw/antro_anon.csv                   (public pipeline entry point)
#
# This step was previously performed BY HAND and was the single largest gap in
# the reproducibility chain: the committed public file could not be regenerated
# from anything, and step 01 read the private file directly, so the pipeline
# could not run from a public deposit at all.
#
# WHAT IS REMOVED
#   nombre  participant name          -> dropped entirely
#   ide     internal record number    -> replaced by a positional anon_id (P001…)
#
# 🔴 THE ide -> anon_id MAPPING IS NEVER WRITTEN. anon_id is derived from row
#    position alone, so re-running reproduces the same labels without any
#    key existing on disk. Do not add a mapping file to the public tree.
#
# The 16 analysis columns are copied unchanged. The script asserts that every
# analysis cell survives the transformation, and — while the hand-made file is
# still present — that the regenerated file agrees with it value for value.
#
# ONLY THIS SCRIPT READS THE PRIVATE SOURCE. Everything downstream, starting
# with 01_clinical_audit.R, runs on data/raw/antro_anon.csv.
############################################################

ROOT <- here::here()

SRC <- file.path(ROOT, "data/raw_private/sleep_raw_original.csv")
OUT <- file.path(ROOT, "data/raw/antro_anon.csv")
ANALYSIS <- c("age", "sex", "height", "weight", "wc", "bf", "mm", "water", "hb",
              paste0("Comp.", 1:7))

stopifnot(file.exists(SRC))
priv <- read.csv(SRC, stringsAsFactors = FALSE, check.names = TRUE,
                 fileEncoding = "UTF-8-BOM")
message(sprintf("private source: %d rows x %d cols", nrow(priv), ncol(priv)))
stopifnot(all(c("ide", "nombre", ANALYSIS) %in% names(priv)))

## ---------- de-identify ----------
anon <- priv[, ANALYSIS, drop = FALSE]
anon <- cbind(anon_id = sprintf("P%03d", seq_len(nrow(priv))), anon,
              stringsAsFactors = FALSE)

## no identifier may survive
stopifnot(!any(c("ide", "nombre") %in% names(anon)))
leak <- vapply(anon, function(col) any(grepl("[A-Za-z]{4,}", as.character(col)) &
                                       !as.character(col) %in% c("female", "male")),
               logical(1))
leak["anon_id"] <- FALSE
stopifnot(!any(leak))
message("de-identified: ", ncol(anon), " columns, no free text outside sex/anon_id")

## ---------- compare with the hand-made file it replaces ----------
if (file.exists(OUT)) {
  old <- read.csv(OUT, stringsAsFactors = FALSE, check.names = TRUE)
  stopifnot(nrow(old) == nrow(anon))
  stopifnot(identical(old$anon_id, anon$anon_id))
  diffs <- 0L
  for (v in ANALYSIS) {
    a <- suppressWarnings(as.numeric(old[[v]]))
    b <- suppressWarnings(as.numeric(anon[[v]]))
    if (v == "sex") { d <- sum(old[[v]] != anon[[v]]) } else {
      d <- sum(xor(is.na(a), is.na(b)) |
               (!is.na(a) & !is.na(b) & abs(a - b) > 1e-9))
    }
    diffs <- diffs + d
  }
  message(sprintf("agreement with the existing hand-made file: %d of %d cells differ",
                  diffs, nrow(anon) * length(ANALYSIS)))
  stopifnot(diffs == 0L)
}

write.csv(anon, OUT, row.names = FALSE, na = "")
message(sprintf("WROTE: %s  (%d rows x %d cols)", OUT, nrow(anon), ncol(anon)))
message("Next step: R/01_clinical_audit.R, which reads this file.")
