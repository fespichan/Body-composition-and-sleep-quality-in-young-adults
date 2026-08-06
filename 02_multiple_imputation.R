## 02_multiple_imputation.R  — PRODUCTION multiple imputation.
## Input : data/processed/sleep_data_clean_for_imputation.csv  (regenerated KEYSTONE)
## Output: data/processed/mice_imputation_m100.rds  (the mids object)
## Config: legacy block predictor matrix, m=100, maxit=50, rf physio / pmm PSQI,
##         sex not imputed, seed 123456.  saveRDS the mids object.
## NOTE  : mice `rf` now uses the ranger engine (published run likely randomForest),
##         so exact reproduction is NOT expected; criterion is CLOSE agreement.
## SCOPE : this script ONLY generates the imputations. The published Table 2 (and
##         Table S2) are produced by R/05_table2_global_model.R, which reads the mids
##         written here. A superseded 5-predictor (bf-containing, n=289) Table 2
##         block was removed from this script on 2026-08-02: the model it fitted
##         is no longer the published model and it printed target values that
##         contradict the manuscript.

suppressPackageStartupMessages(library(mice))
root <- here::here()
logf <- file.path(root, "outputs/reports/02_multiple_imputation.txt")
dir.create(dirname(logf), recursive = TRUE, showWarnings = FALSE)
say <- function(...) { line <- sprintf(...); cat(line, "\n"); cat(line, "\n", file = logf, append = TRUE) }
cat("", file = logf)  # truncate

say("=== MI PRODUCTION RUN ===")
say("START: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("mice version: %s | rf engine backend present -> ranger: %s, randomForest: %s",
    as.character(packageVersion("mice")),
    requireNamespace("ranger", quietly = TRUE),
    requireNamespace("randomForest", quietly = TRUE))

d <- read.csv(file.path(root, "data/processed/sleep_data_clean_for_imputation.csv"),
              stringsAsFactors = FALSE, check.names = TRUE)
comp <- paste0("Comp.", 1:7)
d$sex <- factor(d$sex, levels = c("female","male"))
for (v in setdiff(names(d), "sex")) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))

physio_vars <- c("hb","bf","mm","water","height","weight","wc","age")  # legacy: age in physio block
likert_vars <- comp

ini  <- mice(d, maxit = 0, printFlag = FALSE)
meth <- ini$method
meth[c("hb","bf","mm","water","height","weight","wc")] <- "rf"
meth[likert_vars] <- "pmm"
meth["sex"] <- ""; meth["age"] <- ""    # complete -> predictors, not imputed

pred <- ini$predictorMatrix
pred[,] <- 0
pred[, "sex"] <- 1; pred["sex", ] <- 0
pred[physio_vars, physio_vars] <- 1
pred[physio_vars, likert_vars] <- 0
pred[likert_vars, likert_vars] <- 1
pred[likert_vars, physio_vars] <- 0
diag(pred) <- 0

say("Imputing: m=100, maxit=50, seed=123456 (block predictor matrix) ...")
t0 <- Sys.time()
imp <- mice(d, method = meth, predictorMatrix = pred,
            m = 100, maxit = 50, seed = 123456, printFlag = FALSE)
say("Imputation wall-clock: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

rds <- file.path(root, "data/processed/mice_imputation_m100.rds")
dir.create(dirname(rds), recursive = TRUE, showWarnings = FALSE)   # absent on a clean checkout
saveRDS(imp, rds)
say("Saved mids object: %s", rds)

say("\nNext step: R/05_table2_global_model.R (reads this mids and produces Table 2 + Table S2).")

say("\n=== sessionInfo ===")
si <- capture.output(sessionInfo()); cat(si, sep = "\n", file = logf, append = TRUE); cat("\n", file = logf, append = TRUE)
say("\nEND: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
