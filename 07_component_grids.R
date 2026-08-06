############################################################
# 07_component_grids.R — Component-level grids, OLD vs NEW, + BH-FDR
# (was dp3_component_grids.R; "DP3" was an internal batch code, not a method name)
#
# OLD spec (legacy 09_Component_Screening.R):
#     outcome ~ sex + age + bmic + bf + hb ,  ALL 289 rows
# NEW spec (4-predictor, per-component OBSERVED outcome sets):
#     outcome ~ sex + age + bmic + hb ,  rows where Comp.k observed
#
# Outcome per component: 0 = No_problem  vs  1-3 = Has_problem  (legacy rule)
# Pooling: Rubin's rules via mice::pool on m=100 production imputations.
# DIAGNOSTIC ONLY — writes a log, touches no manuscript.
############################################################

ROOT <- here::here()
suppressMessages({ library(mice) })

LOG <- file.path(ROOT, "outputs/reports/07_component_grids.txt")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)   # absent on a clean checkout
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }

t0 <- Sys.time()
say("=== COMPONENT-LEVEL GRIDS (OLD vs NEW) + BH-FDR ===  %s",
    format(t0, "%Y-%m-%d %H:%M:%S"))

mids <- readRDS(file.path(ROOT, "data/processed/mice_imputation_m100.rds"))
orig <- mids$data
n_imp <- mids$m
say("mids: n = %d, m = %d", nrow(orig), n_imp)

long <- complete(mids, "long")
long$bmi  <- long$weight / ((long$height / 100)^2)
long$bmic <- factor(ifelse(long$bmi < 25, "Normal",
                    ifelse(long$bmi < 30, "Overweight", "Obesity")),
                    levels = c("Normal", "Overweight", "Obesity"))
long$sex  <- factor(long$sex, levels = c("female", "male"))
long$.id  <- as.integer(long$.id)

comps <- paste0("Comp.", 1:7)
cnames <- c("Subjective Quality", "Sleep Latency", "Sleep Duration",
            "Sleep Efficiency", "Sleep Disturbances", "Medication Use",
            "Daytime Dysfunction")

# ---- one pooled component fit -------------------------------------------
fit_component <- function(comp, predictors, rows_idx, label) {
  models <- vector("list", n_imp)
  nfail <- 0L; nwarn <- 0L; maxabs <- 0
  for (i in seq_len(n_imp)) {
    di <- long[long$.imp == i, ]
    di <- di[di$.id %in% rows_idx, ]
    y  <- di[[comp]]
    di$outcome <- factor(ifelse(y == 0, "No_problem", "Has_problem"),
                         levels = c("No_problem", "Has_problem"))
    f <- as.formula(paste("outcome ~", predictors))
    m <- withCallingHandlers(
      tryCatch(glm(f, data = di, family = binomial()),
               error = function(e) { nfail <<- nfail + 1L; NULL }),
      warning = function(w) { nwarn <<- nwarn + 1L; invokeRestart("muffleWarning") })
    if (!is.null(m)) {
      maxabs <- max(maxabs, max(abs(coef(m)[-1]), na.rm = TRUE))
      models[[i]] <- m
    }
  }
  models <- Filter(Negate(is.null), models)
  if (length(models) < 50) { say("  !! %s / %s : only %d fits converged", comp, label, length(models)); return(NULL) }
  s <- summary(pool(as.mira(models)), exponentiate = TRUE, conf.int = TRUE)
  s <- s[s$term != "(Intercept)", ]
  data.frame(comp = comp, name = cnames[match(comp, comps)], spec = label,
             term = as.character(s$term), OR = s$estimate,
             lo = s$`2.5 %`, hi = s$`97.5 %`, p = s$p.value,
             nfit = length(models), nwarn = nwarn, maxabs = maxabs,
             stringsAsFactors = FALSE)
}

# ---- sample sets ---------------------------------------------------------
say("\n--- per-component OBSERVED-outcome sets (from keystone NA pattern) ---")
obs_rows <- list()
for (k in seq_along(comps)) {
  cc <- comps[k]
  idx <- which(!is.na(orig[[cc]]))
  obs_rows[[cc]] <- idx
  y <- orig[[cc]][idx]
  say("  %s (%-19s) observed N = %3d | No_problem = %3d | Has_problem = %3d",
      cc, cnames[k], length(idx), sum(y == 0), sum(y > 0))
}

# assert observed outcomes are invariant across imputations (they must be)
for (cc in comps) {
  idx <- obs_rows[[cc]]
  a <- long[[cc]][long$.imp == 1][match(idx, long$.id[long$.imp == 1])]
  b <- long[[cc]][long$.imp == n_imp][match(idx, long$.id[long$.imp == n_imp])]
  stopifnot(identical(a, b), identical(a, orig[[cc]][idx]))
}
say("  [assert OK] observed component values identical across imputations")

# ---- run both specifications --------------------------------------------
all_rows <- unique(long$.id)
say("\n--- fitting OLD spec (sex+age+bmic+bf+hb, all %d rows) ---", length(all_rows))
OLD <- do.call(rbind, lapply(comps, function(cc)
  fit_component(cc, "sex + age + bmic + bf + hb", all_rows, "OLD")))

say("--- fitting NEW spec (sex+age+bmic+hb, per-component observed rows) ---")
NEW <- do.call(rbind, lapply(comps, function(cc)
  fit_component(cc, "sex + age + bmic + hb", obs_rows[[cc]], "NEW")))

# ---- BH-FDR within each family ------------------------------------------
OLD$q <- p.adjust(OLD$p, method = "BH")
NEW$q <- p.adjust(NEW$p, method = "BH")

pr_grid <- function(D, title) {
  say("\n================================================================")
  say("%s   (family = %d cells)", title, nrow(D))
  say("================================================================")
  say("%-8s %-20s %-16s %9s %-18s %9s %9s %s",
      "Comp", "Name", "Term", "OR", "95% CI", "raw p", "BH q", "")
  for (cc in comps) {
    d <- D[D$comp == cc, ]
    if (!nrow(d)) next
    for (j in seq_len(nrow(d))) {
      flag <- ""
      if (d$p[j] < 0.05) flag <- "  <-- raw p<.05"
      if (d$q[j] < 0.05) flag <- paste0(flag, "  <== SURVIVES FDR")
      say("%-8s %-20s %-16s %9.3f  (%7.3f-%8.3f) %9.4f %9.4f%s",
          if (j == 1) cc else "", if (j == 1) d$name[1] else "", d$term[j],
          d$OR[j], d$lo[j], d$hi[j], d$p[j], d$q[j], flag)
    }
  }
  say("  min raw p = %.4f | min BH q = %.4f | cells with q<0.05: %d",
      min(D$p), min(D$q), sum(D$q < 0.05))
}

pr_grid(OLD, "OLD GRID — 5 predictors incl. body fat, n = 289 (published S3/S4 spec)")
pr_grid(NEW, "NEW GRID — 4 predictors, per-component observed sets")

# ---- old vs new, cell by cell -------------------------------------------
say("\n================================================================")
say("OLD vs NEW, matched cells (bf cells exist only in OLD)")
say("================================================================")
say("%-8s %-16s %10s %10s %10s %10s %10s %10s", "Comp", "Term",
    "OR_old", "p_old", "q_old", "OR_new", "p_new", "q_new")
for (cc in comps) {
  o <- OLD[OLD$comp == cc, ]; n <- NEW[NEW$comp == cc, ]
  for (tm in unique(c(o$term, n$term))) {
    oi <- o[o$term == tm, ]; ni <- n[n$term == tm, ]
    say("%-8s %-16s %10s %10s %10s %10s %10s %10s", cc, tm,
        if (nrow(oi)) sprintf("%.3f", oi$OR) else "-",
        if (nrow(oi)) sprintf("%.4f", oi$p)  else "-",
        if (nrow(oi)) sprintf("%.4f", oi$q)  else "-",
        if (nrow(ni)) sprintf("%.3f", ni$OR) else "-",
        if (nrow(ni)) sprintf("%.4f", ni$p)  else "-",
        if (nrow(ni)) sprintf("%.4f", ni$q)  else "-")
  }
}

# ---- the two published headlines ----------------------------------------
say("\n--- THE TWO PUBLISHED HEADLINE CELLS ---")
h <- rbind(
  cbind(where = "C6 sexmale",         OLD[OLD$comp == "Comp.6" & OLD$term == "sexmale", c("OR","lo","hi","p","q")]),
  cbind(where = "C6 sexmale",         NEW[NEW$comp == "Comp.6" & NEW$term == "sexmale", c("OR","lo","hi","p","q")]),
  cbind(where = "C7 bmicOverweight",  OLD[OLD$comp == "Comp.7" & OLD$term == "bmicOverweight", c("OR","lo","hi","p","q")]),
  cbind(where = "C7 bmicOverweight",  NEW[NEW$comp == "Comp.7" & NEW$term == "bmicOverweight", c("OR","lo","hi","p","q")]))
h$spec <- c("OLD (5-pred, n=289)", "NEW (4-pred, observed)", "OLD (5-pred, n=289)", "NEW (4-pred, observed)")
for (i in seq_len(nrow(h)))
  say("  %-18s %-24s OR %.3f (%.3f-%.3f)  raw p = %.4f   BH q = %.4f",
      h$where[i], h$spec[i], h$OR[i], h$lo[i], h$hi[i], h$p[i], h$q[i])
say("  published reference values: C6 sexmale OR 0.207 (0.055-0.773) p .0195 | C7 overweight OR 2.445 (1.018-5.874) p .0457")

# ---- FDR robustness to family definition --------------------------------
say("\n--- BH-FDR robustness to family definition (NEW grid) ---")
fams <- list(
  "all 35 cells (5 terms x 7 components)"        = rep(TRUE, nrow(NEW)),
  "28 cells, excluding age (4 terms x 7)"        = NEW$term != "age",
  "21 cells, exposures only (sex, 2 BMI)"        = NEW$term %in% c("sexmale","bmicOverweight","bmicObesity"),
  "7 cells, sex only"                            = NEW$term == "sexmale",
  "7 cells, hb only"                             = NEW$term == "hb")
for (nm in names(fams)) {
  sel <- fams[[nm]]
  qq <- p.adjust(NEW$p[sel], method = "BH")
  say("  %-42s k=%2d  min raw p=%.4f  min q=%.4f  survivors=%d",
      nm, sum(sel), min(NEW$p[sel]), min(qq), sum(qq < 0.05))
}

# ---- MANUSCRIPT TABLE 3 (component-level, as published) ------------------
# Multiplicity family = the 28 EXPOSURE cells (sex, 2 BMI contrasts, hb) x 7 components.
# Age is an ADJUSTMENT COVARIATE: its estimates are shown for transparency but it is
# not an exposure of interest, so it is excluded from the family and carries no q.
# q MUST be computed here from full-precision p-values -- never re-derived from the
# rounded values printed in this log (that gives 0.0644 instead of 0.0649).
say("\n=== MANUSCRIPT TABLE 3 (4-predictor, per-component observed sets) ===")
expo <- NEW$term != "age"
NEW$q28 <- NA_real_
NEW$q28[expo] <- p.adjust(NEW$p[expo], method = "BH")
say("  multiplicity family: %d exposure cells | min raw p = %.4f | min q = %.4f | survivors q<0.05 = %d",
    sum(expo), min(NEW$p[expo]), min(NEW$q28[expo]), sum(NEW$q28[expo] < 0.05))
stopifnot(sum(expo) == 28, sum(NEW$q28[expo] < 0.05) == 0)

lab3 <- c(sexmale = "Male sex", bmicOverweight = "Overweight",
          bmicObesity = "Obesity", hb = "Haemoglobin", age = "Age")
ord3 <- c("sexmale", "bmicOverweight", "bmicObesity", "hb", "age")
say("%-24s %-13s %-22s %7s %7s", "Component", "Term", "OR (95% CI)", "p", "q")
for (k in seq_along(comps)) {
  d <- NEW[NEW$comp == comps[k], ]
  d <- d[match(ord3, d$term), ]
  for (j in seq_len(nrow(d))) {
    say("%-24s %-13s %6.2f (%5.2f-%6.2f) %7.3f %7s",
        if (j == 1) sprintf("C%d %s", k, cnames[k]) else "",
        lab3[d$term[j]], d$OR[j], d$lo[j], d$hi[j], d$p[j],
        if (is.na(d$q28[j])) "-" else sprintf("%.3f", d$q28[j]))
  }
}
say("  Observed N per component: %s",
    paste(sprintf("C%d %d", 1:7, sapply(comps, function(cc) length(obs_rows[[cc]]))), collapse = " | "))

# ---- convergence / stability --------------------------------------------
say("\n--- fit stability (separation watch) ---")
for (D in list(OLD, NEW)) {
  lab <- D$spec[1]
  u <- D[!duplicated(D$comp), ]
  for (j in seq_len(nrow(u)))
    say("  %-4s %-7s fits=%3d  glm warnings=%4d  max|beta|=%8.2f%s",
        lab, u$comp[j], u$nfit[j], u$nwarn[j], u$maxabs[j],
        if (u$maxabs[j] > 10) "   <-- SEPARATION" else "")
}

say("\nSurvivors at q<0.05 — OLD: %d of %d | NEW: %d of %d",
    sum(OLD$q < 0.05), nrow(OLD), sum(NEW$q < 0.05), nrow(NEW))
say("\nEND  %s  (%.1f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
