############################################################
# 10_supp_figures.R — Supplementary Figures S1-S4 and S8
#
#   S1  kernel density, observed vs imputed, physiological variables
#   S2  kernel density, observed vs imputed, PSQI components
#   S3  trace plots, MICE convergence, physiological variables
#   S4  trace plots, MICE convergence, PSQI components
#   S8  MFA individual factor map, Component 7 x BMI category
#
# Reads : data/processed/mice_imputation_m100.rds
# Writes: outputs/figures/FigureS1..S4, S8 (.tif)
#         outputs/reports/10_supp_figures.txt
#
# ⚠ FIGURE NUMBERING. The legacy names are on an older scheme and do NOT map
#   one-to-one onto the manuscript. Verified against the current captions:
#       manuscript S1  =  legacy Fig_S2c_DensityPlot_Physiological
#       manuscript S2  =  legacy Fig_S2d_DensityPlot_PSQI_Components
#       manuscript S3  =  legacy Fig_S2a_TracePlot_Physiological
#       manuscript S4  =  legacy Fig_S2b_TracePlot_PSQI_Components
#   Legacy also produced Fig_S1_Missing_Pattern, which is NOT in the manuscript
#   and is therefore not reproduced here. Legacy S5/S6/S7 were deleted with the
#   predictive section.
#
# 🔴 FIGURE S8 IS NOT A FAITHFUL REPRODUCTION OF THE EMBEDDED COPY, DELIBERATELY.
#   Legacy 14_MFA_Comp7.R builds its MFA with group = c(2,9,8,2) and
#   num.group.sup = c(4) — i.e. the SLEEP GROUP IS ACTIVE. That is exactly the
#   circularity removed from Figures 3 and 4 by 06_mfa.R: the outcome helps
#   build the axes on which the outcome is then shown to separate. Reproducing
#   it unchanged would reinstate, in the supplement, the defect corrected in the
#   main text. This script therefore applies the SAME outcome-blind
#   specification as 06_mfa.R. The figure it writes consequently DIFFERS from
#   the copy currently embedded in the supplementary .docx.
#   ⇒ Do not swap it in without a decision. Figure S8's survival is itself an
#     open question, since it belongs to the Component-7 predictive framing.
############################################################

ROOT <- here::here()
suppressPackageStartupMessages({
  library(mice); library(dplyr); library(lattice)
  library(FactoMineR); library(factoextra); library(ggplot2)
})

LOG <- file.path(ROOT, "outputs/reports/10_supp_figures.txt")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
t0 <- Sys.time()
say("=== SUPPLEMENTARY FIGURES S1-S4, S8 ===  %s", format(t0, "%Y-%m-%d %H:%M:%S"))

OUT <- file.path(ROOT, "outputs/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
mids <- readRDS(file.path(ROOT, "data/processed/mice_imputation_m100.rds"))
say("mids loaded: m = %d, maxit = %d, n = %d", mids$m, mids$iteration, nrow(mids$data))

vars_physio <- c("height", "weight", "wc", "bf", "mm", "water", "hb")
vars_psqi   <- paste0("Comp.", 1:7)
stopifnot(all(c(vars_physio, vars_psqi) %in% names(mids$data)))

# Only variables that were actually imputed can appear in a density/trace plot.
imputed <- names(mids$nmis)[mids$nmis > 0]
say("variables with imputed values: %s", paste(imputed, collapse = ", "))
dens_physio  <- intersect(vars_physio, imputed)
dens_psqi    <- intersect(vars_psqi,   imputed)
say("  S1/S3 physiological: %s", paste(dens_physio, collapse = ", "))
say("  S2/S4 PSQI:          %s", paste(dens_psqi,   collapse = ", "))
stopifnot(length(dens_physio) > 0, length(dens_psqi) > 0)

## ---------- S1 / S2 : density, observed vs imputed ----------
# 600 dpi, not 1200: these are dense antialiased panels that compress poorly.
dens <- function(vars, file, w, h) {
  fm <- as.formula(paste("~", paste(vars, collapse = " + ")))
  tiff(file, width = w, height = h, units = "in", res = 600, compression = "lzw")
  print(densityplot(mids, fm, layout = c(4, 2), thicker = 3))
  dev.off()
  say("  written: %-46s %.2f MB", basename(file), file.size(file) / 1e6)
}
say("\n--- S1 / S2 : kernel density, observed (blue) vs imputed (red) ---")
dens(dens_physio, file.path(OUT, "FigureS1_density_physiological.tif"), 8.0, 4.8)
dens(dens_psqi,   file.path(OUT, "FigureS2_density_PSQI_components.tif"), 9.5, 4.8)

## ---------- S3 / S4 : trace plots, convergence ----------
trace <- function(vars, file, w, h) {
  tiff(file, width = w, height = h, units = "in", res = 600, compression = "lzw")
  print(plot(mids, vars, layout = c(2, length(vars))))
  dev.off()
  say("  written: %-46s %.2f MB", basename(file), file.size(file) / 1e6)
}
say("\n--- S3 / S4 : MICE trace plots over %d iterations ---", mids$iteration)
trace(dens_physio, file.path(OUT, "FigureS3_trace_physiological.tif"), 6.7, 7.9)
trace(dens_psqi,   file.path(OUT, "FigureS4_trace_PSQI_components.tif"), 6.7, 7.9)

## ---------- S8 : MFA individual map, Component 7 x BMI category ----------
say("\n--- S8 : MFA individual factor map, Comp.7 x BMI category ---")
say("  SPEC: outcome-blind, matching 06_mfa.R (group = c(1,9,1,8),")
say("        num.group.sup = c(3,4)). Legacy 14 used c(2,9,8,2)/sup=4, i.e.")
say("        the sleep group ACTIVE. See the header note.")

d <- complete(mids, 1)
d <- d %>% mutate(
  bmi  = weight / (height / 100)^2,
  bmic = factor(ifelse(bmi < 25, "Normal", ifelse(bmi < 30, "Overweight", "Obesity")),
                levels = c("Normal", "Overweight", "Obesity")),
  psqi.total = Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 + Comp.7,
  dysf = factor(ifelse(Comp.7 == 0, "No_dysfunction", "Dysfunction"),
                levels = c("No_dysfunction", "Dysfunction")),
  sex = as.factor(sex))
say("  Comp.7 dysfunction: %d of %d (%.1f%%)",
    sum(d$dysf == "Dysfunction"), nrow(d), 100 * mean(d$dysf == "Dysfunction"))

vp  <- c("age", "height", "weight", "bmi", "wc", "bf", "mm", "water", "hb")
vs  <- c(paste0("Comp.", 1:7), "psqi.total")
df  <- d[, c("bmic", vp, "dysf", vs)]
res <- MFA(df, group = c(1, 9, 1, 8), type = c("n", "s", "n", "s"),
           name.group = c("BodyComp", "Physiological", "Comp7 class", "Sleep PSQI"),
           num.group.sup = c(3, 4), graph = FALSE)
eg <- get_eigenvalue(res)[1:2, 2]
say("  variance explained: Dim1 %.1f%% | Dim2 %.1f%%", eg[1], eg[2])

# Panel order is bmic FIRST, dysf second, to match the published caption
# ("Left panel: BMI category; right panel: daytime dysfunction").
# Five colours are required, not four: 3 BMI categories + 2 dysfunction levels.
# fviz_ellipses assigns the palette across the pooled levels ALPHABETICALLY
# (Dysfunction, No_dysfunction, Normal, Obesity, Overweight), not per variable.
p <- fviz_ellipses(res, c("bmic", "dysf"), geom = "point",
                   pointsize = 3.6, labelsize = 3, repel = FALSE,
                   palette = c("#E74C3C", "#2196F3",
                               "#7FBF3F", "#B0651F", "#FFD700")) +
  theme_minimal(base_size = 11) +
  theme(plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title  = element_text(face = "bold", size = 10),
        axis.text   = element_text(size = 9, color = "black"),
        strip.text  = element_text(face = "bold", size = 13),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position = "right") +
  labs(title = "Individual factor map - Component 7 x BMI category (body-composition space)")

# Same two post-hoc layer fixes as Figure 4 in 06_mfa.R, for consistency:
#  (1) alpha on the markers so the 95% confidence ellipses stay visible;
#  (2) drop the text layers from the legend -- a text layer contributes its
#      glyph, the letter "a", as a spurious legend key.
for (i in seq_along(p$layers)) {
  g <- p$layers[[i]]$geom
  if (inherits(g, "GeomPoint")) p$layers[[i]]$aes_params$alpha <- 0.55
  if (inherits(g, "GeomText") || inherits(g, "GeomLabel")) p$layers[[i]]$show.legend <- FALSE
}

f8 <- file.path(OUT, "FigureS8_MFA_Comp7_bmic.tif")
ggsave(f8, p, width = 7, height = 5, units = "in",
       dpi = 1200, compression = "lzw", bg = "white")
say("  written: %-46s %.2f MB", basename(f8), file.size(f8) / 1e6)

say("\nEND  %s  (%.2f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
