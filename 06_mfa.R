############################################################
# MFA — variable-selection structure, OUTCOME-BLIND respecification
#
# Clean port of reference/legacy_scripts/04_MFA_Analysis.R, which is
# READ-ONLY and cannot run standalone (it needs `data_stacked` in memory).
# This script reads the production mids from disk instead.
#
# It runs TWO specifications on IDENTICAL data so the contrast isolates
# one single change -- whether the outcome is ACTIVE or SUPPLEMENTARY:
#
#   SPEC A (legacy, as published):
#     groups = [psqi + bodycomp cat] [9 physio] [Comp.1-7 + psqi.total] [mean,std]
#     active = groups 1,2,3          sup = group 4
#     -> the MFA dimensions are built WITH the outcome (psqi.total loads ~0.89 on Dim1)
#
#   SPEC B (outcome-blind, proposed):
#     groups = [bodycomp cat] [9 physio] [psqi cat] [Comp.1-7 + psqi.total] [mean,std]
#     active = groups 1,2            sup = groups 3,4,5
#     -> dimensions built ONLY from predictors; the outcome is projected, not constitutive
#
# NOTE ON SCOPE: the whole sleep group is made supplementary, not only
# psqi.total. Comp.1-7 ARE the outcome's components; leaving them active would
# still build the dimensions from outcome information and the "did not use the
# outcome" claim would remain untrue.
#
# ** NO variable screening of any kind is performed here. ** The 4-predictor
# explanatory model is fixed. This script only re-specifies what constructs the
# MFA dimensions and redraws the correlation circles.
############################################################

ROOT <- here::here()
suppressPackageStartupMessages({
  library(mice); library(dplyr); library(FactoMineR); library(factoextra)
  library(ggplot2); library(patchwork)
})

LOG <- file.path(ROOT, "outputs/reports/06_mfa.txt")
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }
t0 <- Sys.time()
say("=== MFA OUTCOME-BLIND RESPECIFICATION ===  %s", format(t0, "%Y-%m-%d %H:%M:%S"))

## ---------- data: first imputation, as in legacy 04 ----------
mids <- readRDS(file.path(ROOT, "data/processed/mice_imputation_m100.rds"))
d <- complete(mids, 1)
say("first imputed dataset: n = %d", nrow(d))

d <- d %>%
  mutate(
    bmi        = weight / (height / 100)^2,
    psqi.total = Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6 + Comp.7,
    psqi       = factor(ifelse(psqi.total <= 5, "g.sleep", "p.sleep"),
                        levels = c("g.sleep", "p.sleep")),
    bmic = factor(ifelse(bmi < 25, "Normal", ifelse(bmi < 30, "Overweight", "Obesity")),
                  levels = c("Normal", "Overweight", "Obesity")),
    bfc = factor(case_when(
      sex == "female" & bf >= 5  & bf < 14 ~ "Essential.fat",
      sex == "female" & bf >= 14 & bf < 21 ~ "Athletes",
      sex == "female" & bf >= 21 & bf < 25 ~ "Fitness",
      sex == "female" & bf >= 25 & bf < 32 ~ "Average",
      sex == "female" & bf >= 32           ~ "Obese",
      sex == "male"   & bf >= 2  & bf < 6  ~ "Essential.fat",
      sex == "male"   & bf >= 6  & bf < 14 ~ "Athletes",
      sex == "male"   & bf >= 14 & bf < 18 ~ "Fitness",
      sex == "male"   & bf >= 18 & bf < 25 ~ "Average",
      sex == "male"   & bf >= 25           ~ "Obese")),
    hbc = factor(case_when(
      sex == "female" & hb <  12 ~ "anemia",
      sex == "female" & hb >= 12 ~ "normal",
      sex == "male"   & hb <  13 ~ "anemia",
      sex == "male"   & hb >= 13 ~ "normal")),
    sex = as.factor(sex))

vars_physio <- c("age", "height", "weight", "bmi", "wc", "bf", "mm", "water", "hb")
vars_psqi   <- c(paste0("Comp.", 1:7), "psqi.total")
configs <- c(bfc = "PSQI x Body Fat Category", bmic = "PSQI x BMI Category",
             hbc = "PSQI x Hemoglobin Category", sex = "PSQI x Sex")

run_spec <- function(cat_var, blind) {
  if (blind) {
    # mean/std DROPPED: they were the row-wise mean and SD of all 17 numeric
    # variables pooled across incompatible scales (height ~159 cm, Comp scores 0-3),
    # so `mean` is essentially a proxy for body weight (r = 0.86) with no
    # interpretable meaning. They were also computed partly FROM the outcome
    # (psqi.total is one of the 17). Removing them changes nothing: verified
    # identical eigenvalues and max |change| in physiological loadings = 0.0000.
    df <- d[, c(cat_var, vars_physio, "psqi", vars_psqi)]
    MFA(df, group = c(1, 9, 1, 8), type = c("n", "s", "n", "s"),
        name.group = c("BodyComp", "Physiological", "PSQI class", "Sleep PSQI"),
        num.group.sup = c(3, 4), graph = FALSE)
  } else {
    df <- d[, c("psqi", cat_var, vars_physio, vars_psqi)]
    num <- df[, c(vars_physio, vars_psqi)]
    df$mean <- rowMeans(num, na.rm = TRUE); df$std <- apply(num, 1, sd, na.rm = TRUE)
    MFA(df, group = c(2, 9, 8, 2), type = c("n", "s", "s", "s"),
        name.group = c("Biological & Behavioral", "Physiological", "Sleep PSQI", "Illustrative"),
        num.group.sup = c(4), graph = FALSE)
  }
}

## ---------- compare the two specifications ----------
say("\n--- physiological variable coordinates, SPEC A (active outcome) vs SPEC B (blind) ---")
cmp <- list()
for (cv in names(configs)) {
  A <- run_spec(cv, blind = FALSE); B <- run_spec(cv, blind = TRUE)
  ca <- get_mfa_var(A, "quanti.var")$coord[vars_physio, 1:2]
  cb <- get_mfa_var(B, "quanti.var")$coord[vars_physio, 1:2]
  eA <- get_eigenvalue(A)[1:2, 2]; eB <- get_eigenvalue(B)[1:2, 2]
  say("\n  config %-5s | SPEC A var%%: Dim1 %.1f Dim2 %.1f | SPEC B var%%: Dim1 %.1f Dim2 %.1f",
      cv, eA[1], eA[2], eB[1], eB[2])
  say("  %-8s %18s %18s", "", "SPEC A (Dim1,Dim2)", "SPEC B (Dim1,Dim2)")
  for (v in vars_physio)
    say("  %-8s %9.3f %8.3f %9.3f %8.3f", v, ca[v,1], ca[v,2], cb[v,1], cb[v,2])
  cmp[[cv]] <- list(A = ca, B = cb)
}

## ---------- does the cluster structure (and hence the selection) hold? ----------
say("\n--- CLUSTER STRUCTURE CHECK (adiposity vs lean mass) ---")
adip <- c("weight", "bmi", "wc", "bf"); lean <- c("mm", "water", "height", "hb")
for (cv in names(configs)) {
  B <- cmp[[cv]]$B
  # within-cluster mean pairwise cosine on the Dim1-Dim2 plane
  cosmat <- function(M) { u <- M / sqrt(rowSums(M^2)); u %*% t(u) }
  Cb <- cosmat(B)
  wa <- mean(Cb[adip, adip][upper.tri(diag(length(adip)))])
  wl <- mean(Cb[lean, lean][upper.tri(diag(length(lean)))])
  bt <- mean(Cb[adip, lean])
  say("  %-5s SPEC B: mean cos within adiposity %.3f | within lean %.3f | between %.3f  -> %s",
      cv, wa, wl, bt,
      if (wa > bt && wl > bt) "CLUSTERS HOLD" else "*** CLUSTERS DO NOT SEPARATE ***")
}
say("\n  Representatives implied by the one-per-cluster rule: adiposity -> BMI (category),")
say("  lean mass -> haemoglobin.  NO screening against the outcome was performed.")

## ---------- regenerate the 4 correlation circles, outcome-blind ----------
say("\n--- regenerating correlation circles (SPEC B, outcome supplementary) ---")
options(ggrepel.max.overlaps = Inf)   # every arrow must be labelled; see note below
circles <- list()
for (cv in names(configs)) {
  B <- run_spec(cv, blind = TRUE)
  # NOTE: fviz_mfa_var() does NOT forward max.overlaps to ggrepel, so it must be
  # set globally (see options() call above). Without it ggrepel silently drops 8
  # labels, and the 8 dropped ones are exactly the supplementary PSQI arrows --
  # including psqi.total, the variable this respecification exists to keep visible.
  circles[[cv]] <- fviz_mfa_var(B, "quanti.var", palette = "d3",
                                col.var.sup = "purple", repel = TRUE,
                                arrowsize = 1.0, labelsize = 4,
                                axes.linetype = "solid") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          axis.title = element_text(face = "bold", size = 10),
          axis.text  = element_text(size = 9, color = "black"),
          panel.background = element_rect(fill = "white", color = "black"),
          legend.position = "bottom") +
    labs(title = configs[[cv]])
}
combined <- (circles$bfc | circles$bmic) / (circles$hbc | circles$sex) +
  plot_annotation(tag_levels = "A")

OUT <- file.path(ROOT, "outputs/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUT, "Figure3_MFA_correlation_circles.tif")
# REPRODUCIBILITY: ggrepel places labels by a RANDOMISED optimisation, and it
# draws from the RNG at RENDER time, not build time. Without a seed set here
# this figure differs on every run (verified: three consecutive runs produced
# three different files, while Figure 4 -- which uses repel = FALSE -- was
# byte-identical every time). The seed makes the label layout deterministic;
# it changes nothing about the underlying analysis.
set.seed(20260803)
ggsave(f, combined, width = 7, height = 7, units = "in",
       dpi = 1200, compression = "lzw", bg = "white")
say("  written: %s (%.1f MB, 7x7 in, 1200 dpi, lzw)", f, file.size(f) / 1e6)

## ---------- Figure 4: individual factor map, PSQI x Sex, outcome-blind ----------
# Under SPEC B the axes are built from body composition ONLY, so any separation of
# good vs poor sleepers is a genuine property of the physiological space rather
# than an artefact of the outcome having helped define the axes.
say("\n--- regenerating individual factor map (SPEC B, PSQI x Sex) ---")
Bsex <- run_spec("sex", blind = TRUE)
eS <- get_eigenvalue(Bsex)[1:2, 2]
say("  variance explained: Dim1 %.1f%% | Dim2 %.1f%%", eS[1], eS[2])

sepr <- function(res, grp) {
  co <- res$ind$coord[, 1:2]
  cen <- aggregate(co, by = list(g = grp), FUN = mean)
  d0 <- sqrt(sum((cen[1, 2:3] - cen[2, 2:3])^2))
  w <- mean(sapply(split(as.data.frame(co), grp),
                   function(m) mean(sqrt(rowSums(scale(m, scale = FALSE)^2)))))
  d0 / w
}
say("  separation ratio (centroid distance / within-group spread):")
say("    PSQI  good vs poor : %.3f  %s", sepr(Bsex, d$psqi),
    if (sepr(Bsex, d$psqi) < 1) "-> groups OVERLAP" else "-> groups separate")
say("    Sex   male vs female: %.3f", sepr(Bsex, d$sex))

# pointsize raised from 1.6 to 3.6 (+2) at the author's request, 2026-08-04.
# (+4, i.e. 5.6, was tried first and rejected: at that size the markers cover the
#  95% confidence ellipses that the figure legend explicitly describes.)
p_ind <- fviz_ellipses(Bsex, c("psqi", "sex"), geom = "point",
                       pointsize = 3.6, labelsize = 3, repel = FALSE,
                       palette = c("#2196F3", "#C0FF3E", "#FFD700", "#E74C3C")) +
  theme_minimal(base_size = 11) +
  theme(plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title  = element_text(face = "bold", size = 10),
        axis.text   = element_text(size = 9, color = "black"),
        # panel headings ("psqi", "sex") enlarged at the author's request
        strip.text  = element_text(face = "bold", size = 13),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position = "right") +
  labs(title = "Individual factor map - PSQI x Sex (body-composition space)")

# Two post-hoc layer fixes that fviz_ellipses gives no argument for:
#  (1) ALPHA on the markers, so the 95% confidence ellipses stay visible through
#      the point cloud instead of being buried by it.
#  (2) The stray "a" in the legend. fviz_ellipses draws the group-mean labels
#      with a text layer, and a text layer contributes its glyph -- the letter
#      "a" -- as a legend key. Dropping that layer from the legend removes it;
#      the labels themselves stay on the plot.
for (i in seq_along(p_ind$layers)) {
  g <- p_ind$layers[[i]]$geom
  if (inherits(g, "GeomPoint")) p_ind$layers[[i]]$aes_params$alpha <- 0.55
  if (inherits(g, "GeomText") || inherits(g, "GeomLabel"))
    p_ind$layers[[i]]$show.legend <- FALSE
}
say("  point alpha 0.55; text layers dropped from the legend (removes the stray \"a\")")

f4 <- file.path(OUT, "Figure4_MFA_individual_map.tif")
ggsave(f4, p_ind, width = 7, height = 5, units = "in",
       dpi = 1200, compression = "lzw", bg = "white")
say("  written: %s (%.1f MB, 7x5 in, 1200 dpi, lzw)", f4, file.size(f4) / 1e6)

say("\nEND  %s  (%.1f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
