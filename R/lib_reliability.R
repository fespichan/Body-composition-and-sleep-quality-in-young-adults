# Internal-consistency reliability on observed PSQI components.
#
# Operates on observed component scores only (no imputation). Each estimate is
# computed on the complete-case subset for the requested component set.

complete_cases <- function(components) {
  components[stats::complete.cases(components), , drop = FALSE]
}

cronbach_alpha <- function(components) {
  cc <- complete_cases(components)
  est <- suppressWarnings(psych::alpha(cc, check.keys = FALSE, warnings = FALSE))
  list(alpha = est$total$raw_alpha, n = nrow(cc))
}

# McDonald's omega (total). Maximum-likelihood factor extraction with a
# principal-axis fallback when ML fails to converge.
mcdonald_omega <- function(components) {
  cc <- complete_cases(components)
  fit <- tryCatch(
    psych::omega(cc, nfactors = 1, plot = FALSE, fm = "ml", rotate = "none"),
    error = function(e) NULL)
  if (is.null(fit)) {
    fit <- tryCatch(
      psych::omega(cc, nfactors = 1, plot = FALSE, fm = "pa", rotate = "none"),
      error = function(e) NULL)
  }
  list(omega = if (is.null(fit)) NA_real_ else fit$omega.tot, n = nrow(cc))
}
