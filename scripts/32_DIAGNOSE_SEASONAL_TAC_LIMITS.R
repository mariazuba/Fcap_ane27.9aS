# ============================================================
# 32_DIAGNOSE_SEASONAL_TAC_LIMITS.R
# Diagnose Q2 capture limits in seasonal robustness annual Catch/TAC flags.
# Reads saved validation and simulation objects; runs no simulations.
# Checks all flagged cases, including TAC overshoots, and saves results.
# Equality with a calculated bound is diagnostic evidence, not a trace
# proving which internal function was called. Q2 checks do not allocate
# the full annual shortfall to this season. No settings are changed.
# Usage from project root:
#   source("scripts/32_DIAGNOSE_SEASONAL_TAC_LIMITS.R")
# Select S-C2 with Sys.setenv(MSE_LIMITS_ROBUSTNESS_ID="S-C2").
# This script checks Q2 only; other seasons may also constrain catches.
# Requirements: here, dplyr, FLCore, FLFleet, FLBEIA.
# ============================================================

library(FLBEIA)
library(dplyr)
library(here)

robustness_id <- Sys.getenv("MSE_LIMITS_ROBUSTNESS_ID", "S-C1")
stopifnot(robustness_id %in% c("S-C1", "S-C2"))

validation_file <- here("outputs", "mse", "robustness", "validation",
                        "post_run", robustness_id, "21c_refined_new_validation_summary.rds")
stopifnot(file.exists(validation_file))
validation <- readRDS(validation_file)
TAC_flags <- validation$TAC_flags
stopifnot(is.data.frame(TAC_flags))
if (nrow(TAC_flags) == 0L) stop("No Catch/TAC flags to diagnose for ", robustness_id)
stopifnot(all(startsWith(TAC_flags$scenario_id, paste0(robustness_id, "__"))))
stopifnot(!anyDuplicated(TAC_flags[c("scenario_number", "year", "iter1000")]))

# Controls may retain all 1000 iterations while fleet outputs contain 100.
control_iteration <- function(x, local, global, block_size) {
  n <- dim(x)[6]
  if (n == 1L) return(1L)
  if (n == block_size) return(local)
  if (n == 1000L) return(global)
  stop("Unexpected control iteration dimension: ", n)
}

case_groups <- split(TAC_flags, paste(TAC_flags$scenario_id, TAC_flags$block, sep="__"))

Q2_limits <- bind_rows(lapply(case_groups, function(cases) {
  run_file <- here("data", "mse", "robustness_runs", cases$scenario_id[1],
                   sprintf("block_%02d.rds", cases$block[1]))
  stopifnot(file.exists(run_file))
  message("Checking ", cases$scenario_id[1], " | block ", cases$block[1])
  run <- readRDS(run_file)
  fleet <- run$fleets$SEINE
  biol <- run$biols$ANE
  catches <- fleet@metiers$ALL@catches$ANE
  threshold <- run$fleets.ctrl$catch.threshold

  bind_rows(lapply(seq_len(nrow(cases)), function(i) {
    yr <- as.character(cases$year[i])
    it <- cases$iter_local[i]
    control_it <- control_iteration(threshold, it, cases$iter1000[i], dim(biol@n)[6])
    rho <- as.numeric(threshold[,yr,,"2",,control_it])
    stopifnot(length(rho) == 1L, is.finite(rho), rho > 0, rho <= 1)
    N <- as.numeric(biol@n[,yr,,"2",,it])
    M <- as.numeric(biol@m[,yr,,"2",,it])
    W <- as.numeric(catches@landings.wt[,yr,,"2",,it])
    # The biomass comparison below assumes the current no-discard OM.
    discards <- as.numeric(catches@discards.n[,yr,,"2",,it])
    stopifnot(all(is.finite(discards)), all(abs(discards) < 1e-8))
    catch_Q2 <- sum(as.numeric(catches@landings[,yr,,"2",,it]))
    raw_limit <- sum(rho * N * W)
    limit_Q2 <- sum(rho * N * exp(-M/2) * W)
    effort_Q2 <- as.numeric(fleet@effort[,yr,,"2",,it])
    capacity_Q2 <- as.numeric(fleet@capacity[,yr,,"2",,it])
    stopifnot(all(is.finite(c(N, M, W, catch_Q2, limit_Q2, effort_Q2, capacity_Q2))))

    data.frame(scenario_number=cases$scenario_number[i],
               scenario_id=cases$scenario_id[i], block=cases$block[i],
               year=cases$year[i], iter_local=it, iter1000=cases$iter1000[i],
               rho=rho, catch_Q2=catch_Q2, upper_bound_rho=raw_limit,
               limit_Q2=limit_Q2, effort_Q2=effort_Q2, capacity_Q2=capacity_Q2,
               at_catch_limit=abs(catch_Q2-limit_Q2) <= 1e-6*max(1, limit_Q2),
               at_capacity=abs(effort_Q2-capacity_Q2) <= 1e-6*max(1, capacity_Q2))
  }))
}))

Q2_check <- Q2_limits %>%
  left_join(TAC_flags %>% select(scenario_number, year, iter1000, TAC, Catch, Catch_TAC_ratio),
            by=c("scenario_number", "year", "iter1000")) %>%
  mutate(direction=if_else(Catch_TAC_ratio < 1, "Below TAC", "Above TAC"))

Q2_summary <- Q2_check %>% count(direction, at_catch_limit, at_capacity, name="n_cases")
Q2_by_MP <- Q2_check %>% count(scenario_number, direction, at_catch_limit, at_capacity, name="n_cases")

above_TAC <- TAC_flags %>%
  filter(Catch_TAC_ratio > 1) %>%
  mutate(excess_t=Catch-TAC, excess_pct=100*(Catch_TAC_ratio-1)) %>%
  select(scenario_number, block, year, iter_local, iter1000,
         TAC, Catch, excess_t, excess_pct, HCR_branch)

stopifnot(nrow(Q2_check) == nrow(TAC_flags))
diagnostic_dir <- here("outputs", "mse", "robustness", "diagnostics", robustness_id)
dir.create(diagnostic_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(Q2_check, file.path(diagnostic_dir, paste0(robustness_id, "_Q2_capture_limits.csv")), row.names=FALSE)
write.csv(Q2_summary, file.path(diagnostic_dir, paste0(robustness_id, "_Q2_limits_summary.csv")), row.names=FALSE)
write.csv(Q2_by_MP, file.path(diagnostic_dir, paste0(robustness_id, "_Q2_limits_by_MP.csv")), row.names=FALSE)
write.csv(above_TAC, file.path(diagnostic_dir, paste0(robustness_id, "_TAC_overshoots.csv")), row.names=FALSE)
saveRDS(list(robustness_id=robustness_id, source_validation_file=validation_file, Q2_check=Q2_check,
             Q2_summary=Q2_summary, Q2_by_MP=Q2_by_MP, above_TAC=above_TAC,
             relative_tolerance=1e-6),
        file.path(diagnostic_dir, paste0("32_", robustness_id, "_TAC_limits_diagnostics.rds")))

print(Q2_summary)
print(tibble::as_tibble(above_TAC), n=Inf, width=Inf)
cat("\nDiagnostic outputs saved in:", diagnostic_dir, "\n")
