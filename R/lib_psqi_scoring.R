# PSQI component scoring (Buysse et al., 1989).
#
# These functions implement the deterministic scoring layer only: recoding
# item-level responses into the seven components, the sleep-efficiency
# computation, and the duration plausibility bounds. Free-text parsing of
# reported clock times (bedtime, wake time, sleep hours) is handled upstream
# by the canonical scorer and supplied here as numeric inputs.

# Duration plausibility bounds applied before component 3 scoring.
# Floor and ceiling are physiological-plausibility limits (Hirshkowitz et al.,
# 2015), not statistical choices; they affect only the duration/efficiency
# components.
DURATION_FLOOR_HOURS   <- 4
DURATION_CEILING_HOURS <- 9

# A night's time in bed below/above these bounds is treated as non-computable
# for efficiency (e.g. unparseable or self-contradictory reported times).
HOURS_IN_BED_MIN <- 0
HOURS_IN_BED_MAX <- 16

adjust_sleep_hours <- function(hours) {
  ifelse(is.na(hours), NA_real_,
         pmin(pmax(hours, DURATION_FLOOR_HOURS), DURATION_CEILING_HOURS))
}

hours_in_bed <- function(bedtime, wake_time) {
  ifelse(is.na(bedtime) | is.na(wake_time), NA_real_,
         ifelse(bedtime > wake_time, (24 - bedtime) + wake_time,
                wake_time - bedtime))
}

sleep_efficiency <- function(adjusted_hours, time_in_bed) {
  ok <- !is.na(adjusted_hours) & !is.na(time_in_bed) & time_in_bed > 0
  ifelse(ok, pmin(adjusted_hours / time_in_bed * 100, 100), NA_real_)
}

recode_band <- function(x, breaks_upper, scores) {
  out <- rep(NA_real_, length(x))
  for (k in seq_along(breaks_upper)) {
    hit <- is.na(out) & !is.na(x) & x <= breaks_upper[k]
    out[hit] <- scores[k]
  }
  out[!is.na(x) & is.na(out)] <- scores[length(scores)]
  out
}

score_component_2 <- function(latency_score, early_difficulty_score) {
  recode_band(latency_score + early_difficulty_score,
              breaks_upper = c(0, 2, 4), scores = c(0, 1, 2, 3))
}

score_component_3 <- function(adjusted_hours) {
  ifelse(is.na(adjusted_hours), NA_real_,
         ifelse(adjusted_hours > 7, 0,
                ifelse(adjusted_hours >= 6, 1,
                       ifelse(adjusted_hours >= 5, 2, 3))))
}

score_component_4 <- function(efficiency, computable) {
  c4 <- ifelse(is.na(efficiency), NA_real_,
               ifelse(efficiency >= 85, 0,
                      ifelse(efficiency >= 75, 1,
                             ifelse(efficiency >= 65, 2, 3))))
  ifelse(computable, c4, NA_real_)
}

score_component_5 <- function(disturbance_sum) {
  recode_band(disturbance_sum, breaks_upper = c(0, 9, 18), scores = c(0, 1, 2, 3))
}

score_component_7 <- function(somnolence_score, enthusiasm_score) {
  recode_band(somnolence_score + enthusiasm_score,
              breaks_upper = c(0, 2, 4), scores = c(0, 1, 2, 3))
}

# Assemble the seven components and the global score from parsed inputs.
#
# `parsed` is a data frame with: bedtime, wake_time (numeric clock 0-24),
# bedtime_ambiguous (logical), sleep_hours (numeric, pre-adjustment),
# quality, latency, early_difficulty, disturbance_sum, medication,
# somnolence, enthusiasm (integer item scores 0-3).
score_psqi <- function(parsed, adjust_duration = TRUE) {
  adj <- if (adjust_duration) adjust_sleep_hours(parsed$sleep_hours) else parsed$sleep_hours
  tib <- hours_in_bed(parsed$bedtime, parsed$wake_time)
  eff <- sleep_efficiency(adj, tib)

  c4_computable <- !parsed$bedtime_ambiguous &
    !is.na(tib) & tib > HOURS_IN_BED_MIN & tib <= HOURS_IN_BED_MAX

  components <- data.frame(
    C1 = parsed$quality,
    C2 = score_component_2(parsed$latency, parsed$early_difficulty),
    C3 = score_component_3(adj),
    C4 = score_component_4(eff, c4_computable),
    C5 = score_component_5(parsed$disturbance_sum),
    C6 = parsed$medication,
    C7 = score_component_7(parsed$somnolence, parsed$enthusiasm)
  )

  total <- ifelse(rowSums(is.na(components)) == 0,
                  rowSums(components), NA_real_)
  components$PSQI_Total <- total
  components$Classification <- ifelse(is.na(total), NA_character_,
                                      ifelse(total <= 5, "good", "poor"))
  components
}
