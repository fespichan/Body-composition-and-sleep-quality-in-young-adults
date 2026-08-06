############################################################
# Table 2 — global adjusted model (4 predictors)
#
#   poor sleep (PSQI > 5) ~ sex + age + bmic + hb
#
# PRIMARY      n = 206  (>= 1 observed PSQI component)
# SENSITIVITY  n = 289  (full cohort)
#
# Pooling: Rubin's rules via mice::pool on the m = 100 production
# imputations (data/processed/mice_imputation_m100.rds).
#
# Also emits: Nagelkerke R2 (mean +/- SD), good-sleeper events, EPV,
# and the Supplementary Table S2 VIF values (mean +/- SD).
#
# Every value published in the manuscript is HARD-ASSERTED at the end.
# The script regenerates them; it does not hardcode them.
#
# Provenance note: this script reconstructs the computation whose original
# output log was written on 2026-07-29 (then named table2_4pred_n206.txt,
# renamed to 05_table2_global_model.txt when the scripts were numbered). The
# generating script was never committed and was lost with the scratchpad;
# this file restores it. Assertions are against that log and the .docx.
############################################################

ROOT <- here::here()
suppressMessages({ library(mice); library(car) })

LOG <- file.path(ROOT, "outputs/reports/05_table2_global_model.txt")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)   # absent on a clean checkout
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }

t0 <- Sys.time()
say("=== TABLE 2 - GLOBAL ADJUSTED MODEL (4 predictors) ===  %s",
    format(t0, "%Y-%m-%d %H:%M:%S"))

mids  <- readRDS(file.path(ROOT, "data/processed/mice_imputation_m100.rds"))
orig  <- mids$data
n_imp <- mids$m
comps <- paste0("Comp.", 1:7)

## ---------- analytic samples, from the OBSERVED component pattern ----------
n_obs   <- rowSums(!is.na(orig[, comps]))
id_all  <- seq_len(nrow(orig))
id_206  <- id_all[n_obs >= 1]
n_all   <- length(id_all)
n_206   <- length(id_206)
n_112   <- sum(n_obs == 7)
n_92    <- sum(n_obs == 6)
n_2     <- sum(n_obs == 5)
n_83    <- sum(n_obs == 0)

say("\n--- analytic samples (from observed-component pattern) ---")
say("  full cohort                          n = %d", n_all)
say("  >=1 observed component  (PRIMARY)    n = %d", n_206)
say("    all seven observed                 n = %3d", n_112)
say("    six of seven observed              n = %3d   (missing Comp.4 only)", n_92)
say("    five of seven observed             n = %3d   (missing Comp.3 AND Comp.4)", n_2)
say("  no PSQI item completed (excluded)    n = %2d", n_83)
part <- orig[id_206, ]   # tally missing components WITHIN the analytic sample only
say("  partials n = %d; missing-component tally: Comp.4 x%d | Comp.3 x%d",
    n_92 + n_2, sum(is.na(part$Comp.4)), sum(is.na(part$Comp.3)))
stopifnot(n_206 == n_112 + n_92 + n_2, n_all - n_206 == n_83)
say("  [assert OK] %d = %d + %d + %d ; %d - %d = %d",
    n_206, n_112, n_92, n_2, n_all, n_206, n_83)

## ---------- long form + derived variables ----------
long <- complete(mids, "long")
long$bmi  <- long$weight / ((long$height / 100)^2)
long$bmic <- factor(ifelse(long$bmi < 25, "Normal",
                    ifelse(long$bmi < 30, "Overweight", "Obesity")),
                    levels = c("Normal", "Overweight", "Obesity"))
long$sex  <- factor(long$sex, levels = c("female", "male"))
long$.id  <- as.integer(long$.id)
long$psqi_total <- rowSums(long[, comps])
long$poor       <- as.integer(long$psqi_total > 5)

## ---------- one pooled fit ----------
fit_global <- function(rows_idx, label) {
  fits  <- vector("list", n_imp)
  nagel <- numeric(n_imp)
  good  <- numeric(n_imp)
  vifs  <- matrix(NA_real_, n_imp, 4,
                  dimnames = list(NULL, c("sex", "age", "bmic", "hb")))
  gvif_adj <- matrix(NA_real_, n_imp, 4,
                     dimnames = list(NULL, c("sex", "age", "bmic", "hb")))
  for (i in seq_len(n_imp)) {
    di <- long[long$.imp == i, ]
    di <- di[di$.id %in% rows_idx, ]
    m1 <- glm(poor ~ sex + age + bmic + hb, data = di, family = binomial())
    m0 <- glm(poor ~ 1, data = di, family = binomial())
    n  <- nobs(m1); L1 <- as.numeric(logLik(m1)); L0 <- as.numeric(logLik(m0))
    cs <- 1 - exp((2 / n) * (L0 - L1))
    nagel[i] <- cs / (1 - exp((2 / n) * L0))
    good[i]  <- sum(di$poor == 0)
    # car::vif() returns a MATRIX when any term has >1 df (here BMI category, 2 df):
    # columns are GVIF, Df, GVIF^(1/(2*Df)). Column 1 is the GVIF reported in
    # Table S2; column 3 is the scale-comparable quantity (Fox & Monette), which
    # is the one that can be read against the conventional VIF thresholds.
    v <- car::vif(m1)
    vifs[i, ] <- if (is.matrix(v)) v[, 1] else v
    if (is.matrix(v)) gvif_adj[i, ] <- v[, 3]
    fits[[i]] <- m1
  }
  pl <- summary(pool(as.mira(fits)), conf.int = TRUE)
  pl <- pl[pl$term != "(Intercept)", ]

  say("\n================================================================")
  say("%s   (n = %d)", label, length(rows_idx))
  say("================================================================")
  say("%-18s %6s  %-26s %8s", "Term", "OR", "95% CI", "p")
  for (i in seq_len(nrow(pl))) {
    say("%-18s %6.3f  (%.3f-%.3f)%s %8.4f",
        pl$term[i], exp(pl$estimate[i]), exp(pl$`2.5 %`[i]), exp(pl$`97.5 %`[i]),
        strrep(" ", 12), pl$p.value[i])
  }
  say("  Nagelkerke R2      %.3f +/- %.3f", mean(nagel), sd(nagel))
  say("  good-sleeper events %.1f (range %d-%d) | EPV = %.2f on %d coefficients",
      mean(good), min(good), max(good), mean(good) / nrow(pl), nrow(pl))
  say("  VIF (mean +/- SD): %s",
      paste(sprintf("%s %.2f+/-%.3f", colnames(vifs),
                    colMeans(vifs), apply(vifs, 2, sd)), collapse = " | "))
  say("  GVIF^(1/(2*Df)) scale-comparable: %s",
      paste(sprintf("%s %.3f", colnames(gvif_adj), colMeans(gvif_adj)), collapse = " | "))
  list(pl = pl, nagel = mean(nagel), vif = colMeans(vifs),
       vadj = colMeans(gvif_adj), good = mean(good))
}

primary <- fit_global(id_206, "PRIMARY - n = 206 (>=1 observed PSQI component)")
sens    <- fit_global(id_all, "SENSITIVITY - n = 289 (full cohort)")

## ---------- assertions against the manuscript as applied (BATCH-BF) ----------
say("\n--- assert against Table 2 as APPLIED to the manuscript (BATCH-BF) ---")
or_of <- function(res, term) exp(res$pl$estimate[res$pl$term == term])
targets <- list(
  list("sexmale",        or_of(primary, "sexmale"),        0.060,  "manuscript"),
  list("age",            or_of(primary, "age"),            1.007,  "manuscript"),
  list("bmicOverweight", or_of(primary, "bmicOverweight"), 2.568,  "manuscript"),
  list("bmicObesity",    or_of(primary, "bmicObesity"),    0.659,  "manuscript"),
  list("hb",             or_of(primary, "hb"),             1.433,  "manuscript"),
  list("VIF sex",        unname(primary$vif["sex"]),       2.68,   "Table S2  "),
  list("VIF age",        unname(primary$vif["age"]),       1.05,   "Table S2  "),
  list("VIF bmic",       unname(primary$vif["bmic"]),      1.11,   "Table S2  "),
  list("VIF hb",         unname(primary$vif["hb"]),        2.59,   "Table S2  "),
  list("NagelkerkeR2",   primary$nagel,                    0.226,  "manuscript"))
ok <- TRUE
for (t in targets) {
  d <- abs(t[[2]] - t[[3]])
  pass <- d < 0.006
  ok <- ok && pass
  say("  %-16s script %.4f  vs %s %.4f   delta %.4f  %s",
      t[[1]], t[[2]], t[[4]], t[[3]], d, if (pass) "ok" else "*** MISMATCH ***")
}
say("\n  ALL COEFFICIENTS REGENERATE FROM THIS SCRIPT: %s", if (ok) "YES" else "NO")
stopifnot(ok)

say("\nEND  %s  (%.1f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
