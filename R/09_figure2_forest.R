############################################################
# 09_figure2_forest.R — Figure 2: forest plot of the adjusted model
#
#   poor sleep (PSQI > 5) ~ sex + age + BMI category + haemoglobin
#   PRIMARY analytic sample n = 206, pooled over m = 100 by Rubin's rules.
#
# Reads : data/processed/mice_imputation_m100.rds
# Writes: outputs/figures/Figure2_forest.tif
#         outputs/reports/09_figure2_forest.txt
#
# The plotted estimates are RECOMPUTED here from the imputations and asserted
# against the published Table 2 before anything is drawn. Figure and table
# therefore cannot drift apart.
#
# ⚠ TWO RULES THAT MUST BE PRESERVED — both are bugs that already occurred once:
#
#  (1) THREE-DECIMAL CONFIDENCE BOUNDS. The haemoglobin lower bound is 0.997.
#      Printed at two decimals it reads "1.00", which makes a null result
#      (p = .052, CI includes 1) look significant. fmt_ci() below forces three
#      decimals whenever a bound rounds to 1.00 at two.
#
#  (2) NEVER INDEX A LABEL VECTOR WITH A FACTOR. `labels[pl$term]` where
#      pl$term is a factor indexes by LEVEL CODE, not by name, and silently
#      shifts every label one row off its estimate. Always as.character().
#
# Provenance: the figure embedded in the manuscript was generated on 2026-07-26
# by a scratchpad script that was lost. This reconstructs it. ⚠ The embedded
# copy is 150 dpi; this writes 1200 dpi, so it is NOT byte-identical to the
# .docx copy. Re-embedding is a separate, approval-gated step.
############################################################

ROOT <- here::here()
suppressMessages(library(mice))

LOG <- file.path(ROOT, "outputs/reports/09_figure2_forest.txt")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
t0 <- Sys.time()
say("=== FIGURE 2 - FOREST PLOT (4 predictors, n = 206) ===  %s",
    format(t0, "%Y-%m-%d %H:%M:%S"))

## ---------- refit exactly as 05_table2_global_model.R does ----------
mids  <- readRDS(file.path(ROOT, "data/processed/mice_imputation_m100.rds"))
orig  <- mids$data
n_imp <- mids$m
comps <- paste0("Comp.", 1:7)
id_206 <- which(rowSums(!is.na(orig[, comps])) >= 1)
stopifnot(length(id_206) == 206)

long <- complete(mids, "long")
long$bmi  <- long$weight / ((long$height / 100)^2)
long$bmic <- factor(ifelse(long$bmi < 25, "Normal",
                    ifelse(long$bmi < 30, "Overweight", "Obesity")),
                    levels = c("Normal", "Overweight", "Obesity"))
long$sex  <- factor(long$sex, levels = c("female", "male"))
long$.id  <- as.integer(long$.id)
long$poor <- as.integer(rowSums(long[, comps]) > 5)

fits <- vector("list", n_imp)
for (i in seq_len(n_imp)) {
  di <- long[long$.imp == i, ]
  fits[[i]] <- glm(poor ~ sex + age + bmic + hb,
                   data = di[di$.id %in% id_206, ], family = binomial())
}
pl <- summary(pool(as.mira(fits)), conf.int = TRUE)
pl <- pl[pl$term != "(Intercept)", ]
pl$term <- as.character(pl$term)          # RULE (2): never leave this a factor
pl$or <- exp(pl$estimate); pl$lo <- exp(pl$`2.5 %`); pl$hi <- exp(pl$`97.5 %`)

## ---------- assert against the published Table 2 BEFORE drawing ----------
PUB <- c(sexmale = 0.060, age = 1.007, bmicOverweight = 2.568,
         bmicObesity = 0.659, hb = 1.433)
say("\n--- assert plotted estimates against published Table 2 ---")
ok <- TRUE
for (tm in names(PUB)) {
  got <- pl$or[pl$term == tm]; d <- abs(got - PUB[[tm]]); pass <- d < 0.006
  ok <- ok && pass
  say("  %-16s script %.4f  vs published %.4f   delta %.4f  %s",
      tm, got, PUB[[tm]], d, if (pass) "ok" else "*** MISMATCH ***")
}
stopifnot(ok)
say("  all five coefficients match; safe to draw")

## ---------- labels and formatting ----------
LAB <- c(sexmale        = "Male sex (vs female)",
         age            = "Age (per year)",
         bmicOverweight = "Overweight (vs normal)",
         bmicObesity    = "Obesity (vs normal)",
         hb             = "Hemoglobin (per g/dL)")   # American spelling, as in the text
ord <- c("sexmale", "age", "bmicOverweight", "bmicObesity", "hb")

# PALETTE. Colour encodes the DOMAIN a predictor belongs to -- demographic,
# adiposity, haematological -- and nothing else.
# 🔴 IT MUST NEVER ENCODE STATISTICAL SIGNIFICANCE. Colouring the "significant"
#    estimate differently is the same act as the asterisks removed from the
#    component tables: a dichotomous verdict smuggled in as formatting, in a
#    paper that reports confidence intervals and makes no such claims.
# Every row also carries a text label, so identity is never colour-alone.
# Hues validated for colour-vision deficiency (worst adjacent pair dE 11.0
# deutan, 25.8 normal vision) against a light surface.
DOM <- c(sexmale = "Demographic", age = "Demographic",
         bmicOverweight = "Adiposity", bmicObesity = "Adiposity",
         hb = "Haematological")
PAL <- c(Demographic = "#2979FF", Adiposity = "#FF6D00", Haematological = "#00BFA5")
pl  <- pl[match(ord, pl$term), ]
pl$label <- LAB[as.character(pl$term)]     # RULE (2) again: character, not factor

# RULE (1): a bound that would render as exactly 1.00 at two decimals is printed
# at three instead, so a CI that includes the null can never look like it excludes
# it. Applied PER BOUND, matching the manuscript, which prints "0.997-2.06".
fmt_row <- function(or, lo, hi) {
  b <- function(x) sprintf(if (abs(round(x, 2) - 1) < 1e-9) "%.3f" else "%.2f", x)
  sprintf("%.2f (%s–%s)", or, b(lo), b(hi))
}
fmt_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)
pl$txt <- mapply(fmt_row, pl$or, pl$lo, pl$hi)
pl$ptx <- vapply(pl$p.value, fmt_p, "")

say("\n--- as plotted ---")
for (i in seq_len(nrow(pl)))
  say("  %-24s %-22s p = %s", pl$label[i], pl$txt[i], pl$ptx[i])
stopifnot(any(grepl("0.997", pl$txt)))     # the three-decimal rule actually fired
say("  [assert OK] haemoglobin CI rendered at three decimals (0.997), not 1.00")

## ---------- draw ----------
OUT <- file.path(ROOT, "outputs/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUT, "Figure2_forest.tif")

k    <- nrow(pl)
ypos <- rev(seq_len(k))
xlim <- range(c(pl$lo, pl$hi, 1)); xlim <- c(xlim[1] * 0.75, xlim[2] * 1.35)
ticks <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8)
ticks <- ticks[ticks >= xlim[1] & ticks <= xlim[2]]

pl$dom <- DOM[as.character(pl$term)]        # character, never a factor
pl$col <- PAL[pl$dom]

FIG_W <- 7; FIG_H <- 4.2
tiff(f, width = FIG_W, height = FIG_H, units = "in", res = 1200, compression = "lzw")
op <- par(mar = c(5.8, 12.2, 1.6, 10.2), xaxs = "i")
plot(NA, xlim = log(xlim), ylim = c(0.4, k + 0.7), axes = FALSE, xlab = "", ylab = "")
abline(v = 0, lty = 2, col = "grey45", lwd = 1.1)
segments(log(pl$lo), ypos, log(pl$hi), ypos, lwd = 1.9, col = pl$col)
segments(log(pl$lo), ypos - 0.12, log(pl$lo), ypos + 0.12, lwd = 1.6, col = pl$col)
segments(log(pl$hi), ypos - 0.12, log(pl$hi), ypos + 0.12, lwd = 1.6, col = pl$col)
points(log(pl$or), ypos, pch = 21, bg = pl$col, col = "white", cex = 1.15, lwd = 0.9)
axis(1, at = log(ticks), labels = ticks, cex.axis = 0.62, tcl = -0.25, mgp = c(2, 0.4, 0))
mtext("Odds ratio (95% CI), log scale", side = 1, line = 1.9, cex = 0.68)
mtext(pl$label, side = 2, at = ypos, las = 1, adj = 1, line = 0.4, cex = 0.66)
# Right-hand columns. Drawn as TWO separate mtext columns rather than one
# space-padded string: padding only lines up in a monospace face, and mixing a
# monospace column into an otherwise proportional figure looks like a different
# typeface. Two columns keep one face throughout and stay aligned.
mtext(pl$txt, side = 4, at = ypos, las = 1, adj = 0, line = 0.4, cex = 0.62, font = 2)
mtext(pl$ptx, side = 4, at = ypos, las = 1, adj = 0, line = 7.4, cex = 0.62, font = 2)
mtext(c("OR (95% CI)", "p"), side = 4, at = k + 0.55, las = 1, adj = 0,
      line = c(0.4, 7.4), cex = 0.58, font = 3, col = "grey35")
mtext(sprintf("Adjusted model, n = %d, pooled over %d imputations", length(id_206), n_imp),
      side = 3, line = 0.2, adj = 0, cex = 0.62, font = 3)
# legend: colour = domain, stated explicitly so it cannot be read as a verdict.
# Placed BELOW the axis title. The y offset is expressed in margin lines converted
# to user units, so it cannot drift into the axis label if the y-range changes.
dm  <- unique(pl$dom)
# height of one margin line, in user units
lin <- diff(par("usr")[3:4]) / par("pin")[2] * (par("mai")[1] / par("mar")[1])
legend(x = mean(par("usr")[1:2]), y = par("usr")[3] - 3.4 * lin, xjust = 0.5,
       legend = dm, col = PAL[dm], pt.bg = PAL[dm], pch = 21,
       pt.cex = 1.1, pt.lwd = 0.9, horiz = TRUE, bty = "n",
       cex = 0.58, xpd = NA, x.intersp = 0.7, text.col = "grey20")
par(op); dev.off()

say("\n  written: %s (%.2f MB, %g x %g in, 1200 dpi, lzw)",
    f, file.size(f) / 1e6, FIG_W, FIG_H)
say("\nEND  %s  (%.2f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
