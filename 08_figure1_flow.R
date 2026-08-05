############################################################
# 08_figure1_flow.R — Figure 1: participant flow
#
#   289 assessed  ->  83 excluded (no PSQI item completed)  ->  206 analytic,
#   of which 112 had all seven components observed and 94 had a partial PSQI
#   (92 missing Component 4 only; 2 missing Components 3 and 4).
#
# Reads : data/processed/sleep_data_clean_for_imputation.csv  (the keystone)
# Writes: outputs/figures/Figure1_participant_flow.tif
#         outputs/reports/08_figure1_flow.txt
#
# Every count is DERIVED from the observed-component pattern and asserted; no
# number in the diagram is typed in by hand.
#
# Provenance: the figure embedded in the manuscript was generated on 2026-07-26
# by a scratchpad script that was lost. This file reconstructs it from the same
# counts. ⚠ The embedded copy is 150 dpi; this script writes 1200 dpi at 7 in,
# so the output is NOT byte-identical to what is currently in the .docx.
# Re-embedding is a separate, approval-gated step.
#
# Base graphics deliberately: the diagram is a handful of boxes and arrows, and
# a reader should be able to follow the geometry line by line.
############################################################

ROOT <- here::here()

LOG <- file.path(ROOT, "outputs/reports/08_figure1_flow.txt")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
t0 <- Sys.time()
say("=== FIGURE 1 - PARTICIPANT FLOW ===  %s", format(t0, "%Y-%m-%d %H:%M:%S"))

## ---------- counts, derived from the observed-component pattern ----------
d <- read.csv(file.path(ROOT, "data/processed/sleep_data_clean_for_imputation.csv"),
              stringsAsFactors = FALSE, check.names = TRUE)
comps <- paste0("Comp.", 1:7)
n_obs <- rowSums(!is.na(d[, comps]))

n_total   <- nrow(d)
n_excl    <- sum(n_obs == 0)
n_analytic<- sum(n_obs >= 1)
n_complete<- sum(n_obs == 7)
n_partial <- sum(n_obs %in% c(5, 6))
n_miss4   <- sum(n_obs == 6)
n_miss34  <- sum(n_obs == 5)
n_female  <- sum(d$sex == "female")
n_male    <- sum(d$sex == "male")

say("\n  assessed                          %3d  (female %d, male %d)",
    n_total, n_female, n_male)
say("  excluded: no PSQI item completed  %3d", n_excl)
say("  analytic sample                   %3d", n_analytic)
say("    all seven components observed   %3d", n_complete)
say("    partial PSQI                    %3d", n_partial)
say("      missing Component 4 only      %3d", n_miss4)
say("      missing Components 3 and 4    %3d", n_miss34)

stopifnot(
  n_total    == 289,
  n_excl     ==  83,
  n_analytic == 206,
  n_complete == 112,
  n_miss4    ==  92,
  n_miss34   ==   2,
  n_partial  ==  94,
  n_analytic == n_complete + n_partial,
  n_total    == n_analytic + n_excl,
  n_female + n_male == n_total,
  # the missing component is Component 4 for every partial
  sum(is.na(d$Comp.4[n_obs >= 1])) == n_partial,
  sum(is.na(d$Comp.3[n_obs >= 1])) == n_miss34)
say("\n  [assert OK] 289 = 206 + 83 ; 206 = 112 + 94 ; 94 = 92 + 2")
say("  [assert OK] Comp.4 missing in all %d partials; Comp.3 additionally in %d",
    n_partial, n_miss34)

## ---------- the diagram ----------
OUT <- file.path(ROOT, "outputs/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUT, "Figure1_participant_flow.tif")

# PALETTE. Two colours only, and they encode the two ROLES a box can have -- on
# the analytic path, or excluded. Colour never encodes a result here. Text stays
# black in every box: identity is carried by the wording, never by colour alone.
# Both hues validated for colour-vision deficiency against a light surface.
COL_PATH   <- "#0072B2"; FILL_PATH   <- "#E8F1F8"   # blue      - on the analytic path
COL_EXCL   <- "#D55E00"; FILL_EXCL   <- "#FBECE3"   # vermilion - excluded
COL_ARROW  <- "#33526B"

CEX_BOX  <- 0.85   # text inside the boxes
CEX_NOTE <- 0.70   # the footnote
PAD_X    <- 5.0    # horizontal padding inside a box, in user units
LINE_H   <- 4.6    # vertical space per line of text

# BOXES SIZE THEMSELVES FROM THEIR TEXT. Hand-tuned widths silently overflowed
# the moment the font size went up (that is exactly what happened at CEX_BOX
# 0.62 -> 0.85). Measuring the strings instead means the geometry can never
# disagree with the type size again, whatever it is set to.
box_w <- function(lines, cex = CEX_BOX)
  max(strwidth(lines, cex = cex)) + 2 * PAD_X
box_h <- function(lines) length(lines) * LINE_H + 3.0

box <- function(x, y, lines, cex = CEX_BOX, border = COL_PATH, fill = FILL_PATH) {
  w <- box_w(lines, cex); h <- box_h(lines)
  rect(x - w/2, y - h/2, x + w/2, y + h/2, border = border, lwd = 1.6, col = fill)
  n  <- length(lines)
  ys <- y + (h/2) - 1.5 - LINE_H * (seq_len(n) - 0.5)
  for (i in seq_len(n))
    text(x, ys[i], lines[i], cex = cex, font = if (i == 1) 2 else 1)
  invisible(c(l = x - w/2, r = x + w/2, t = y + h/2, b = y - h/2))
}
arrow_v <- function(x, y0, y1) arrows(x, y0, x, y1, length = 0.06, lwd = 1.4, col = COL_ARROW)

FIG_W <- 8.0; FIG_H <- 5.6
tiff(f, width = FIG_W, height = FIG_H, units = "in", res = 1200, compression = "lzw")
op <- par(mar = c(0.2, 0.2, 0.2, 0.2))
plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "", asp = NA)

L_ASSESS   <- c("Participants assessed", sprintf("n = %d", n_total),
                sprintf("%d female, %d male", n_female, n_male))
L_EXCL     <- c("Excluded", "no PSQI item completed",
                sprintf("n = %d", n_excl))
L_ANALYTIC <- c("Analytic sample", sprintf("n = %d", n_analytic),
                "at least one PSQI component observed")
L_COMPLETE <- c("Complete PSQI", sprintf("n = %d", n_complete), "all seven components")
L_PARTIAL  <- c("Partial PSQI", sprintf("n = %d", n_partial),
                sprintf("Component 4 missing (n = %d)", n_miss4),
                sprintf("Components 3 and 4 missing (n = %d)", n_miss34))

X_MAIN <- 33; X_EXCL <- 79; X_LEFT <- 16; X_RIGHT <- 56
b1 <- box(X_MAIN, 89, L_ASSESS)
b2 <- box(X_EXCL, 72, L_EXCL, border = COL_EXCL, fill = FILL_EXCL)
b3 <- box(X_MAIN, 57, L_ANALYTIC)
b4 <- box(X_LEFT, 22, L_COMPLETE)
b5 <- box(X_RIGHT, 22, L_PARTIAL)

# geometry guards: nothing may overlap and nothing may leave the canvas
stopifnot(b4["r"] < b5["l"],            # the two bottom boxes
          b1["r"] < b2["l"],            # assessed vs excluded
          b3["r"] < b2["l"],            # analytic vs excluded
          b4["l"] > 0, b5["r"] < 100, b2["r"] < 100, b1["l"] > 0)

arrow_v(X_MAIN, b1["b"], b3["t"])
arrows(X_MAIN, 72, b2["l"], 72, length = 0.06, lwd = 1.4, col = COL_EXCL)
arrow_v(X_MAIN, b3["b"], 38)
segments(X_LEFT, 38, X_RIGHT, 38, lwd = 1.4, col = COL_ARROW)
arrow_v(X_LEFT, 38, b4["t"]); arrow_v(X_RIGHT, 38, b5["t"])

text(50, 4, paste("Component 4 (habitual sleep efficiency) could not be computed when",
                  "reported bed and rise times were ambiguous."),
     cex = CEX_NOTE, font = 2)
par(op); dev.off()
say("  [assert OK] ninguna caja se solapa ni se sale del lienzo")

say("\n  written: %s (%.2f MB, %g x %g in, 1200 dpi, lzw)",
    f, file.size(f) / 1e6, FIG_W, FIG_H)
say("\nEND  %s  (%.2f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
