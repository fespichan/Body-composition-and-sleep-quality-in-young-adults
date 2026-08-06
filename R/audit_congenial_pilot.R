## audit_congenial_pilot.R — RETIRED DIAGNOSTIC, not part of the numbered run order.
## DIAGNOSTIC PILOT (not a replacement analysis, not the production pipeline).
## Question: is the block-design predictor matrix attenuating the physio->PSQI
## associations (esp. hemoglobin) toward the null relative to a congenial matrix?
##
## Design: run TWO imputations on the SAME input (master_clean.csv), identical
## settings (m=20, maxit=10, seed=123456), differing ONLY in the predictor matrix:
##   BLOCK     = legacy two-block structure (physio & PSQI do NOT predict each other)
##   CONGENIAL = mice default full matrix (physio & PSQI predict each other)
## Then fit the Table 2 model on each and pool (Rubin). Compare vs observed-only.
##
## NOTE: master_clean.csv does NOT carry the manuscript's full relational-audit
## missingness (bf ~14% here vs ~48% published), so the published block-design
## Table 2 cannot be exactly reproduced here; it is an EXTERNAL reference only.
## The decisive contrast is BLOCK-pilot vs CONGENIAL-pilot on identical data.

suppressPackageStartupMessages({ library(mice) })
set.seed(123456)

d <- read.csv(here::here("data/processed/master_clean.csv"),
              stringsAsFactors = FALSE, check.names = TRUE, fileEncoding = "UTF-8-BOM")

comp <- paste0("Comp.", 1:7)
raw  <- c("sex","age","height","weight","wc","bf","mm","water","hb", comp)
dat  <- d[, raw]
dat$sex <- factor(dat$sex, levels = c("female","male"))
for (v in setdiff(raw, "sex")) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))

physio_vars <- c("hb","bf","mm","water","height","weight","wc","age")  # legacy: age in physio block
likert_vars <- comp

## ----- methods (identical for both runs) -----
ini  <- mice(dat, maxit = 0, printFlag = FALSE)
meth <- ini$method
meth[c("hb","bf","mm","water","height","weight","wc")] <- "rf"
meth[likert_vars] <- "pmm"
meth["sex"] <- ""; meth["age"] <- ""   # complete -> predictors, not imputed

## ----- CONGENIAL predictor matrix (mice default: all predict all, diag 0) -----
pred_cong <- ini$predictorMatrix   # default is 1 off-diagonal, 0 on diagonal

## ----- BLOCK predictor matrix (legacy two-block) -----
pred_block <- pred_cong
pred_block[,] <- 0
pred_block[, "sex"] <- 1; pred_block["sex", ] <- 0
pred_block[physio_vars, physio_vars] <- 1
pred_block[physio_vars, likert_vars] <- 0
pred_block[likert_vars, likert_vars] <- 1
pred_block[likert_vars, physio_vars] <- 0
diag(pred_block) <- 0

run_imp <- function(pred, label) {
  imp <- mice(dat, method = meth, predictorMatrix = pred,
              m = 20, maxit = 10, seed = 123456, printFlag = FALSE)
  long <- complete(imp, "long", include = FALSE)
  long$bmi  <- long$weight / ((long$height/100)^2)
  long$bmic <- factor(ifelse(long$bmi < 25, "Normal",
                      ifelse(long$bmi < 30, "Overweight", "Obesity")),
                      levels = c("Normal","Overweight","Obesity"))
  long$psqi_total <- rowSums(long[, comp])
  long$poor <- as.integer(long$psqi_total > 5)
  fits <- lapply(sort(unique(long$.imp)), function(k) {
    glm(poor ~ sex + age + bmic + bf + hb,
        data = long[long$.imp == k, ], family = binomial())
  })
  pl <- summary(pool(as.mira(fits)), conf.int = TRUE)
  out <- data.frame(term = as.character(pl$term),
                    OR = round(exp(pl$estimate), 3),
                    CI_low = round(exp(pl$`2.5 %`), 3),
                    CI_high = round(exp(pl$`97.5 %`), 3),
                    p = signif(pl$p.value, 3), row.names = NULL)
  cat(sprintf("\n===== %s MI pilot (m=20, maxit=10, seed=123456) =====\n", label))
  print(out); invisible(out)
}

b <- run_imp(pred_block, "BLOCK (legacy two-block)")
c_ <- run_imp(pred_cong,  "CONGENIAL (full matrix)")

## ----- observed-only full model on identical derived variables (n=89) -----
obs <- dat
obs$bmi  <- obs$weight / ((obs$height/100)^2)
obs$bmic <- factor(ifelse(obs$bmi < 25, "Normal",
                   ifelse(obs$bmi < 30, "Overweight","Obesity")),
                   levels = c("Normal","Overweight","Obesity"))
obs$psqi_total <- ifelse(rowSums(is.na(obs[, comp])) == 0, rowSums(obs[, comp]), NA)
obs$poor <- as.integer(obs$psqi_total > 5)
mobs <- glm(poor ~ sex + age + bmic + bf + hb, data = obs, family = binomial())
oo <- data.frame(term = names(coef(mobs)), OR = round(exp(coef(mobs)),3),
                 p = signif(summary(mobs)$coefficients[,4],3), row.names = NULL)
cat(sprintf("\n===== OBSERVED-ONLY full model (complete-case, n=%d, events=%d) =====\n",
            nrow(mobs$model), sum(mobs$model$poor==0)))
print(oo)

cat("\n===== PUBLISHED block MI Table 2 (EXTERNAL ref; different audit) =====\n")
cat("sexmale 0.08(0.02-0.37)p.002 | age 1.00(.95-1.05)p.98 | Overwt 1.89(.76-4.70)p.17\n")
cat("Obesity 0.81(.29-2.27)p.69 | bf 0.96(.87-1.06)p.42 | hb 1.23(.90-1.68)p.20\n")
