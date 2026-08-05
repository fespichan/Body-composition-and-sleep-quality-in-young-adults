## 01_clinical_audit.R
## Clean port of reference/legacy_scripts/01_Clinical_Audit.R (do NOT edit legacy in place).
## Three-stage physiological data validation prior to multiple imputation.
##
## Input : data/raw/antro_anon.csv  (PUBLIC, de-identified; written by 00_deidentify.R)
## Output: data/processed/sleep_data_clean_for_imputation.csv
##
## Audit logic is IDENTICAL to legacy; the only deviations are:
##   (1) input is the PUBLIC de-identified file. Until 2026-08-04 this script read
##       data/raw_private/sleep_raw_original.csv directly, which meant the pipeline
##       could not be run from a public deposit at all -- the chain broke on its
##       first line. 00_deidentify.R now performs that step and is the ONLY script
##       that touches the private source. The keystone is byte-identical either way.
##   (2) BUG FIX: legacy L282 read the not-yet-written output file to compute
##       "audit impact"; here na_after is computed from the in-memory audited data.
##   (3) legacy dropped a `name` column that does not exist; the identifier in the
##       public file is `anon_id` -> drop that.
##   (4) xlsx exports omitted (writexl not required); writes the CSV keystone only.

suppressPackageStartupMessages(library(dplyr))
root <- here::here()

SRC <- file.path(root, "data/raw/antro_anon.csv")
datos <- read.csv(SRC, stringsAsFactors = FALSE, check.names = TRUE,
                  fileEncoding = "UTF-8-BOM")
n_original <- nrow(datos)
message(sprintf("Observations loaded: %d", n_original))

## pre-existing missingness (computed BEFORE any NA-setting) — bug-fix source of truth
na_original <- sapply(datos[, c("hb","bf","mm","water")], function(x) sum(is.na(x)))

datos <- datos %>% mutate(bmi = weight / ((height/100)^2))

## ---------- STAGE 1 — absolute range validation ----------
datos <- datos %>% mutate(
  hb    = ifelse(hb < 7  | hb > 22,  NA, hb),
  bf    = ifelse(bf < 3  | bf > 60,  NA, bf),
  mm    = ifelse(mm < 10 | mm > 70,  NA, mm),
  water = ifelse(water < 30 | water > 75, NA, water),
  bmi   = ifelse(bmi < 12 | bmi > 55, NA, bmi)
)
na_stage1 <- sapply(datos[, c("hb","bf","mm","water","bmi")], function(x) sum(is.na(x)))
message("\n--- STAGE 1 (absolute) cumulative NA ---"); print(na_stage1)

## ---------- STAGE 2 — sex-specific coherence ----------
datos <- datos %>% mutate(
  flag_bf = case_when(
    sex=="female" & !is.na(bf) & (bf<10 | bf>50) ~ 1,
    sex=="male"   & !is.na(bf) & (bf<5  | bf>45) ~ 1, TRUE ~ 0),
  flag_mm = case_when(
    sex=="female" & !is.na(mm) & (mm<15 | mm>50) ~ 1,
    sex=="male"   & !is.na(mm) & (mm<25 | mm>65) ~ 1, TRUE ~ 0),
  flag_water = case_when(
    sex=="female" & !is.na(water) & (water<35 | water>60) ~ 1,
    sex=="male"   & !is.na(water) & (water<43 | water>65) ~ 1, TRUE ~ 0),
  flag_hb = case_when(
    sex=="female" & !is.na(hb) & (hb<9.5  | hb>17.5) ~ 1,
    sex=="male"   & !is.na(hb) & (hb<10.5 | hb>19.0) ~ 1, TRUE ~ 0)
)
datos <- datos %>% mutate(
  bf    = ifelse(flag_bf==1, NA, bf),
  mm    = ifelse(flag_mm==1, NA, mm),
  water = ifelse(flag_water==1, NA, water),
  hb    = ifelse(flag_hb==1, NA, hb)
)
message("\n--- STAGE 2 (sex-specific) flags ---")
print(data.frame(flag_bf=sum(datos$flag_bf), flag_mm=sum(datos$flag_mm),
                 flag_water=sum(datos$flag_water), flag_hb=sum(datos$flag_hb)))

## ---------- STAGE 3 — relational coherence ----------
datos <- datos %>% mutate(
  expected_water = case_when(
    sex=="female" & !is.na(bf) ~ 70 - (bf*0.70),
    sex=="male"   & !is.na(bf) ~ 75 - (bf*0.60), TRUE ~ NA_real_),
  flag_rel_bf_mm = case_when(
    !is.na(bf) & !is.na(mm) & !is.na(weight) &
      ((bf/100*weight)+mm)/weight*100 > 90 ~ 1, TRUE ~ 0),
  flag_rel_water = case_when(
    !is.na(water) & !is.na(expected_water) & abs(water-expected_water) > 12 ~ 1, TRUE ~ 0),
  flag_rel_bmi_bf = case_when(
    !is.na(bmi) & !is.na(bf) & bmi>30 &
      ((sex=="female" & bf<18) | (sex=="male" & bf<12)) ~ 1,
    !is.na(bmi) & !is.na(bf) & bmi<20 & bf>35 ~ 1, TRUE ~ 0),
  flag_relational = ifelse(flag_rel_bf_mm==1 | flag_rel_water==1 | flag_rel_bmi_bf==1, 1, 0)
)
datos <- datos %>% mutate(
  water = ifelse(flag_rel_water==1, NA, water),
  bf    = ifelse(flag_rel_bmi_bf==1, NA, bf),
  bf    = ifelse(flag_rel_bf_mm==1, NA, bf),
  mm    = ifelse(flag_rel_bf_mm==1, NA, mm)
)
message("\n--- STAGE 3 (relational) flags ---")
print(data.frame(flag_rel_bf_mm=sum(datos$flag_rel_bf_mm), flag_rel_water=sum(datos$flag_rel_water),
                 flag_rel_bmi_bf=sum(datos$flag_rel_bmi_bf), flag_relational=sum(datos$flag_relational)))

## ---------- AUDIT REPORT ----------
na_final <- data.frame(
  Variable = c("hb","bf","mm","water","bmi"),
  NAs_total = sapply(datos[, c("hb","bf","mm","water","bmi")], function(x) sum(is.na(x))),
  Pct_missing = round(sapply(datos[, c("hb","bf","mm","water","bmi")],
                             function(x) mean(is.na(x))*100), 1))
message("\n--- FINAL missing summary ---"); print(na_final, row.names = FALSE)

affected <- datos %>% filter(flag_bf==1|flag_mm==1|flag_water==1|flag_hb==1|flag_relational==1)
message(sprintf("\nTotal affected participants: %d of %d (%.1f%%)",
                nrow(affected), n_original, nrow(affected)/n_original*100))

## BUG-FIXED audit impact: na_after from in-memory data (NOT a re-read of the output file)
na_after <- sapply(datos[, c("hb","bf","mm","water")], function(x) sum(is.na(x)))
message("\n--- AUDIT IMPACT ---")
message(sprintf("Pre-existing missing (hb+bf+mm+water): %d", sum(na_original)))
message(sprintf("Audit-generated missing:               %d", sum(na_after) - sum(na_original)))
message(sprintf("Total missing:                         %d", sum(na_after)))

## per-participant count of individual values set to missing by the audit (across hb,bf,mm,water)
was_na    <- is.na(read.csv(SRC, fileEncoding="UTF-8-BOM")[, c("hb","bf","mm","water")])
now_na    <- is.na(datos[, c("hb","bf","mm","water")])
audit_set <- now_na & !was_na
message(sprintf("Individual values set to NA by audit (hb+bf+mm+water): %d", sum(audit_set)))
message(sprintf("Participants with >=1 audit-set value:                 %d", sum(rowSums(audit_set)>0)))
message(sprintf("Participants flagged by relational stage only:         %d", sum(datos$flag_relational==1)))

## ---------- EXPORT (de-identified) ----------
out <- datos %>% dplyr::select(-bmi, -flag_bf, -flag_mm, -flag_water, -flag_hb,
                               -flag_rel_bf_mm, -flag_rel_water, -flag_rel_bmi_bf,
                               -flag_relational, -expected_water, -anon_id)
op <- file.path(root, "data/processed/sleep_data_clean_for_imputation.csv")
# data/processed/ is gitignored and therefore absent from a fresh checkout; an
# empty directory is not committed either. Without this the very first step of
# the pipeline fails on a clean clone.
dir.create(dirname(op), recursive = TRUE, showWarnings = FALSE)
write.csv(out, op, row.names = FALSE)
message(sprintf("\nWROTE: %s  (%d rows x %d cols)", op, nrow(out), ncol(out)))
