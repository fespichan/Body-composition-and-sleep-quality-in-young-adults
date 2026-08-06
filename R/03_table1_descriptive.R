############################################################
# 03_table1_descriptive.R — Table 1
#
#   Demographic, anthropometric and haematological
#   characteristics by sex, on the PRE-IMPUTATION cohort (N = 289).
#
# Reads : data/processed/sleep_data_clean_for_imputation.csv  (the keystone,
#         written by 01_clinical_audit.R)
# Writes: outputs/reports/03_table1_descriptive.txt
#
# Table 1 is stratified by SEX (0% missing) rather than by PSQI (61% missing),
# to use the whole cohort and avoid selection bias. Continuous variables are
# summarised as median (Q1, Q3) on OBSERVED values only; categorical as n (%)
# within the stratum. Tests: Wilcoxon rank-sum (continuous), Pearson's
# chi-squared (categorical) — the same tests named in the table footnote.
#
# Muscle mass is deliberately ABSENT: it is not a Table 1 row in the
# manuscript (see the published table), although it is measured and imputed.
#
# Every published cell is HARD-ASSERTED at the end. The script regenerates the
# values; it does not hardcode them. Replaces legacy 05_Table1_Descriptive.R,
# which depended on gtsummary and wrote a .docx; this version has no package
# dependency beyond base R so that a reader can follow every computation.
############################################################

ROOT <- here::here()

LOG <- file.path(ROOT, "outputs/reports/03_table1_descriptive.txt")
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
con <- file(LOG, open = "wt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = con) }

t0 <- Sys.time()
say("=== TABLE 1 - CHARACTERISTICS BY SEX (pre-imputation) ===  %s",
    format(t0, "%Y-%m-%d %H:%M:%S"))

## ---------- data ----------
d <- read.csv(file.path(ROOT, "data/processed/sleep_data_clean_for_imputation.csv"),
              stringsAsFactors = FALSE, check.names = TRUE)
comps <- paste0("Comp.", 1:7)
for (v in setdiff(names(d), "sex")) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
d$sex <- factor(d$sex, levels = c("female", "male"))

d$bmi <- d$weight / (d$height / 100)^2
d$bmic <- factor(ifelse(is.na(d$bmi), NA,
                 ifelse(d$bmi < 25, "Normal",
                 ifelse(d$bmi < 30, "Overweight", "Obesity"))),
                 levels = c("Normal", "Overweight", "Obesity"))
d$hbc <- factor(ifelse(is.na(d$hb), NA,
                ifelse((d$sex == "female" & d$hb < 12) |
                       (d$sex == "male"   & d$hb < 13), "Anemia", "Normal")),
                levels = c("Normal", "Anemia"))
# PSQI global score exists only when ALL seven components are observed
d$psqi.total <- ifelse(rowSums(is.na(d[, comps])) == 0, rowSums(d[, comps]), NA)

say("cohort loaded: n = %d  (female %d, male %d)",
    nrow(d), sum(d$sex == "female"), sum(d$sex == "male"))

## ---------- helpers ----------
# QUANTILE DEFINITION. The published table was produced with gtsummary, whose
# quantiles here correspond to stats::quantile TYPE 2 (the SAS definition:
# averaging at discontinuities), NOT R's default type 7. Verified empirically:
# type 2 reproduces 63 of 63 published continuous cells, type 7 only 58.
# The five cells type 7 gets wrong are wc female Q3, wc male Q1, bf male Q1,
# water female median and water female Q3.
QTYPE <- 2

# ROUNDING. Two published cells sit exactly on a half (water female median
# 45.250, Q3 47.050). R's sprintf uses round-half-to-even and prints 45.2 /
# 47.0; the published table shows 45.3 / 47.1, i.e. round-half-away-from-zero.
# This function reproduces the published convention.
round_half_up <- function(x, digits = 1) {
  z <- abs(x) * 10^digits
  sign(x) * floor(z + 0.5 + 1e-8) / 10^digits
}
f1 <- function(x) formatC(round_half_up(x, 1), format = "f", digits = 1)

mq <- function(x) {                       # median (Q1, Q3) on observed values
  x <- x[!is.na(x)]
  q <- stats::quantile(x, c(.25, .5, .75), type = QTYPE, names = FALSE)
  c(med = q[2], q1 = q[1], q3 = q[3])
}
fmt_mq  <- function(s) sprintf("%s (%s, %s)", f1(s["med"]), f1(s["q1"]), f1(s["q3"]))
fmt_p   <- function(p) if (is.na(p)) "" else if (p < 0.001) "<0.001" else
                       if (p >= 0.1) sprintf("%.1f", p) else sprintf("%.3f", p)

RES <- list()   # collected values, for the assertion block

cont_row <- function(var, label) {
  x <- d[[var]]; ok <- !is.na(x)
  s_all <- mq(x); s_f <- mq(x[d$sex == "female"]); s_m <- mq(x[d$sex == "male"])
  p <- suppressWarnings(stats::wilcox.test(x ~ d$sex)$p.value)
  RES[[var]] <<- list(N = sum(ok), all = s_all, f = s_f, m = s_m, p = p)
  say("%-26s %4d  %-20s %-20s %-20s %8s",
      label, sum(ok), fmt_mq(s_all), fmt_mq(s_f), fmt_mq(s_m), fmt_p(p))
}

cat_row <- function(var, label) {
  x <- d[[var]]; ok <- !is.na(x)
  tab <- table(x[ok], d$sex[ok])                    # levels x sex
  p <- suppressWarnings(stats::chisq.test(tab)$p.value)
  RES[[var]] <<- list(N = sum(ok), tab = tab, p = p)
  say("%-26s %4d  %-20s %-20s %-20s %8s", label, sum(ok), "", "", "", fmt_p(p))
  n_all <- sum(ok); n_f <- sum(ok & d$sex == "female"); n_m <- sum(ok & d$sex == "male")
  # percentage style follows the published table: one decimal below 10%, none above
  pc <- function(k, n) {
    p <- 100 * k / n
    sprintf("%d (%s%%)", k,
            if (p < 10) formatC(round_half_up(p, 1), format = "f", digits = 1)
            else        formatC(round_half_up(p, 0), format = "f", digits = 0))
  }
  for (lv in rownames(tab)) {
    a <- sum(tab[lv, ]); f <- tab[lv, "female"]; m <- tab[lv, "male"]
    say("  %-24s %4s  %-20s %-20s %-20s",
        lv, "", pc(a, n_all), pc(f, n_f), pc(m, n_m))
  }
}

## ---------- the table ----------
say("\n%-26s %4s  %-20s %-20s %-20s %8s",
    "Variable", "N", sprintf("Overall (N=%d)", nrow(d)),
    sprintf("female (N=%d)", sum(d$sex == "female")),
    sprintf("male (N=%d)", sum(d$sex == "male")), "p")
say("%s", strrep("-", 108))
cont_row("age",        "Age (years)")
cont_row("bmi",        "BMI (kg/m2)")
cat_row ("bmic",       "BMI Category")
cont_row("wc",         "Waist Circumference (cm)")
cont_row("bf",         "Body Fat (%)")
cont_row("water",      "Total Body Water (%)")
cont_row("hb",         "Hemoglobin (g/dL)")
cat_row ("hbc",        "Anemia Status")
cont_row("psqi.total", "PSQI Global Score")
say("%s", strrep("-", 108))
say("Median (Q1, Q3); n (%%). Wilcoxon rank-sum test; Pearson's Chi-squared test.")

## ---------- assertions against the published Table 1 ----------
# EVERY cell of the published table is checked, at the precision it is printed.
say("\n--- assert against Table 1 as published (every cell) ---")
ok_all <- TRUE; n_ok <- 0; n_tot <- 0
chk <- function(what, got, want, tol = 0) {
  pass <- !is.na(got) && abs(got - want) <= tol
  ok_all <<- ok_all && pass; n_tot <<- n_tot + 1; n_ok <<- n_ok + pass
  if (!pass) say("  %-32s script %9.3f  vs published %9.3f   *** MISMATCH ***",
                 what, got, want)
}
# --- observed N per variable ---
chk("N age", RES$age$N, 289); chk("N bmi", RES$bmi$N, 270); chk("N bmic", RES$bmic$N, 270)
chk("N wc", RES$wc$N, 269);   chk("N bf", RES$bf$N, 150);   chk("N water", RES$water$N, 261)
chk("N hb", RES$hb$N, 260);   chk("N hbc", RES$hbc$N, 260); chk("N psqi", RES$psqi.total$N, 112)
# --- all 63 continuous cells: median, Q1, Q3 x {overall, female, male} ---
PUB <- list(
  age        = list(all = c(21.0, 19.0, 24.0), f = c(21.0, 19.0, 24.0), m = c(21.0, 19.0, 24.0)),
  bmi        = list(all = c(25.2, 23.1, 28.2), f = c(25.3, 23.0, 28.2), m = c(25.2, 23.3, 28.1)),
  wc         = list(all = c(81.0, 75.5, 88.0), f = c(79.0, 74.0, 86.0), m = c(88.0, 80.0, 93.2)),
  bf         = list(all = c(29.6, 25.0, 31.8), f = c(31.0, 29.1, 32.3), m = c(21.5, 19.3, 23.2)),
  water      = list(all = c(46.5, 44.5, 50.9), f = c(45.3, 44.1, 47.1), m = c(53.0, 52.0, 54.4)),
  hb         = list(all = c(13.1, 12.3, 14.5), f = c(12.7, 12.0, 13.3), m = c(15.3, 14.5, 15.9)),
  psqi.total = list(all = c(9.0, 6.0, 12.0),   f = c(10.0, 7.0, 12.0),  m = c(7.0, 5.0, 9.0)))
for (v in names(PUB)) for (g in c("all", "f", "m")) {
  got <- RES[[v]][[g]]; want <- PUB[[v]][[g]]
  for (k in seq_len(3))
    chk(sprintf("%s/%s/%s", v, g, c("med", "Q1", "Q3")[k]),
        round_half_up(got[k], 1), want[k], 0.001)
}
# --- categorical counts ---
chk("bmic Normal f",     RES$bmic$tab["Normal", "female"],      92)
chk("bmic Normal m",     RES$bmic$tab["Normal", "male"],        35)
chk("bmic Overweight f", RES$bmic$tab["Overweight", "female"],  68)
chk("bmic Overweight m", RES$bmic$tab["Overweight", "male"],    32)
chk("bmic Obesity f",    RES$bmic$tab["Obesity", "female"],     35)
chk("bmic Obesity m",    RES$bmic$tab["Obesity", "male"],        8)
chk("hbc Normal f",      RES$hbc$tab["Normal", "female"],      143)
chk("hbc Normal m",      RES$hbc$tab["Normal", "male"],         67)
chk("hbc Anemia f",      RES$hbc$tab["Anemia", "female"],       44)
chk("hbc Anemia m",      RES$hbc$tab["Anemia", "male"],          6)
# --- p-values, at the precision the table prints them ---
chk("p age",   round_half_up(RES$age$p, 1),        0.8,   0.001)
chk("p bmi",   round_half_up(RES$bmi$p, 1),        0.8,   0.001)
chk("p bmic",  round_half_up(RES$bmic$p, 1),       0.3,   0.001)
chk("p hbc",   RES$hbc$p,                          0.008, 0.0005)
chk("p psqi",  RES$psqi.total$p,                   0.002, 0.0005)
for (v in c("wc", "bf", "water", "hb"))
  chk(sprintf("p %s < 0.001", v), as.numeric(RES[[v]]$p < 0.001), 1)

say("  %d of %d published cells regenerate exactly", n_ok, n_tot)
say("\n  ALL PUBLISHED TABLE 1 VALUES REGENERATE: %s", if (ok_all) "YES" else "NO")
stopifnot(ok_all)

say("\nEND  %s  (%.2f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))
close(con)
