# ============================================================
# 30_COMPARE_ROBUSTNESS_WITH_REF.R
# Paired performance comparison of robustness MPs against REF
# Anchovy ane.27.9aS, Gulf of Cadiz
#
# Reads existing evaluations from scripts 22b and 29; no simulations
# or performance recalculation. Matches MPs by their original REF ID.
# Differences are robustness minus REF, in each metric's native units.
# Relative differences (%) are omitted when REF is zero.
# Probability differences remain fractions (0.01 = 1 percentage point).
# Annual maximum risk and ever-below risk are separate metrics.
# Trajectory comparisons retain iter1000, OM and replicate pairing.
# Conditional open-year metrics may be NA when there are no open years.
# No ranking, confidence intervals or figures are produced.
#
# R console from project root:
#   source("scripts/30_COMPARE_ROBUSTNESS_WITH_REF.R") # default R-OM
#   Sys.setenv(MSE_COMPARE_ROBUSTNESS_ID="S-C1")
#   source("scripts/30_COMPARE_ROBUSTNESS_WITH_REF.R")
# Inputs:
#   data/mse/performance_refined/refined_candidate_MP_performance.rds
#   data/mse/robustness_performance/<ID>/robustness_MP_performance.rds
# Outputs:
#   outputs/mse/robustness/comparison/<ID>/
#   30_performance_comparison.csv, 30_full_horizon_comparison.csv,
#   30_annual_comparison.csv, 30_trajectory_comparison.csv,
#   30_comparison_with_REF.rds
# Re-running replaces only these comparison outputs.
# ============================================================

rm(list=ls())
library(dplyr)
library(here)

robustness_id <- Sys.getenv("MSE_COMPARE_ROBUSTNESS_ID", "R-OM")
stopifnot(robustness_id %in% c("R-OM", "S-C1", "S-C2", "I0", "I1"))
ref_file <- here("data", "mse", "performance_refined", "refined_candidate_MP_performance.rds")
rob_file <- here("data", "mse", "robustness_performance", robustness_id, "robustness_MP_performance.rds")
output_dir <- here("outputs", "mse", "robustness", "comparison", robustness_id)
stopifnot(file.exists(ref_file), file.exists(rob_file))
ref <- readRDS(ref_file)
rob <- readRDS(rob_file)
stopifnot(identical(rob$robustness_id, robustness_id))
required_tables <- c("scenario_grid", "performance", "trajectory_metrics", "annual_metrics")
stopifnot(all(required_tables %in% names(ref)), all(required_tables %in% names(rob)))

keys <- rob$scenario_grid %>% dplyr::select(reference_scenario_number, reference_scenario_id, scenario_id, Fcap, Besc)
stopifnot(nrow(keys)==12, !anyDuplicated(keys$reference_scenario_id), !anyDuplicated(keys$scenario_id))
stopifnot(all(keys$reference_scenario_id %in% ref$scenario_grid$scenario_id))
ref_settings <- ref$scenario_grid[match(keys$reference_scenario_id, ref$scenario_grid$scenario_id), ]
stopifnot(all(keys$Fcap==ref_settings$Fcap), all(keys$Besc==ref_settings$Besc))

# Check equal MP settings and horizon boundaries before comparing values.
ref_perf <- ref$performance %>% filter(scenario_id %in% keys$reference_scenario_id)
rob_perf <- rob$performance
stopifnot(nrow(ref_perf)==48, nrow(rob_perf)==48)
setting_names <- c("Fcap", "Besc", "Blim", "Bpa", "first_year", "last_year")
stopifnot(all(setting_names %in% names(ref_perf)), all(setting_names %in% names(rob_perf)))
setting_check <- rob_perf %>%
  dplyr::select(reference_scenario_id, horizon, all_of(setting_names)) %>%
  left_join(ref_perf %>% dplyr::select(scenario_id, horizon, all_of(setting_names)),
            by=c("reference_scenario_id"="scenario_id", "horizon"), suffix=c("_rob", "_REF"))
stopifnot(nrow(setting_check)==48)
for (name in setting_names) stopifnot(all(setting_check[[paste0(name,"_rob")]]==setting_check[[paste0(name,"_REF")]]))

# One row per original REF MP and comparison key is mandatory.
compare_table <- function(ref_table, rob_table, extra_keys, excluded) {
  ref_table <- ref_table[ref_table$scenario_id %in% keys$reference_scenario_id, , drop=FALSE]
  join_keys <- c("reference_scenario_id", extra_keys)
  stopifnot(all(c("scenario_id", extra_keys) %in% names(ref_table)))
  stopifnot(all(join_keys %in% names(rob_table)))
  ref_table$reference_scenario_id <- ref_table$scenario_id
  stopifnot(!anyDuplicated(ref_table[join_keys]), !anyDuplicated(rob_table[join_keys]))
  stopifnot(nrow(ref_table)==nrow(rob_table))
  stopifnot(nrow(anti_join(rob_table, ref_table, by=join_keys))==0)
  stopifnot(nrow(anti_join(ref_table, rob_table, by=join_keys))==0)
  common <- intersect(names(ref_table), names(rob_table))
  metrics <- setdiff(common, c(join_keys, excluded))
  metrics <- metrics[vapply(metrics, function(x) is.numeric(ref_table[[x]]) && is.numeric(rob_table[[x]]), logical(1))]
  stopifnot(length(metrics)>0)
  left <- rob_table %>% dplyr::select(all_of(join_keys), all_of(metrics))
  right <- ref_table %>% dplyr::select(all_of(join_keys), all_of(metrics))
  out <- left_join(left, right, by=join_keys, suffix=c("_rob", "_REF"))
  for (metric in metrics) {
    x <- out[[paste0(metric,"_rob")]]
    baseline <- out[[paste0(metric,"_REF")]]
    out[[paste0(metric,"_delta")]] <- x-baseline
    relative <- rep(NA_real_, length(baseline))
    valid <- is.finite(x) & is.finite(baseline) & baseline!=0
    relative[valid] <- 100*(x[valid]-baseline[valid])/baseline[valid]
    out[[paste0(metric,"_relative_pct")]] <- relative
  }
  out <- out %>% left_join(keys, by="reference_scenario_id") %>% mutate(robustness_id=robustness_id, .before=1)
  stopifnot(nrow(out)==nrow(rob_table))
  out
}

excluded <- c("scenario_number", "reference_scenario_number", "Fcap", "Besc", "Blim", "Bpa", "first_year", "last_year")
performance_comparison <- compare_table(ref$performance, rob$performance, "horizon", excluded)
annual_comparison <- compare_table(ref$annual_metrics, rob$annual_metrics, c("horizon","year"), c(excluded,"n"))
trajectory_comparison <- compare_table(ref$trajectory_metrics, rob$trajectory_metrics,
                                      c("horizon","iter1000","om","replicate"), excluded)
stopifnot(nrow(performance_comparison)==48, nrow(annual_comparison)==720, nrow(trajectory_comparison)==48000)
full_comparison <- performance_comparison %>% filter(horizon=="Full") %>% arrange(Fcap, Besc)

# Compact table: raw probability differences, tonnes and original IAV units.
summary_full <- full_comparison %>% dplyr::select(
  reference_scenario_number, Fcap, Besc,
  max_P_Blim_REF, max_P_Blim_rob, max_P_Blim_delta,
  mean_catch_REF, mean_catch_rob, mean_catch_delta,
  mean_P_closed_REF, mean_P_closed_rob, mean_P_closed_delta,
  median_IAV_REF, median_IAV_rob, median_IAV_delta)
print(summary_full, n=Inf, width=Inf)
cat("\nAnnual risk criterion (maximum annual P_Blim <= 0.05):\n")
cat("REF MPs meeting criterion:", sum(full_comparison$max_P_Blim_REF<=0.05), "/ 12\n")
cat(robustness_id, "MPs meeting criterion:", sum(full_comparison$max_P_Blim_rob<=0.05), "/ 12\n")

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(performance_comparison, file.path(output_dir,"30_performance_comparison.csv"), row.names=FALSE)
write.csv(full_comparison, file.path(output_dir,"30_full_horizon_comparison.csv"), row.names=FALSE)
write.csv(annual_comparison, file.path(output_dir,"30_annual_comparison.csv"), row.names=FALSE)
write.csv(trajectory_comparison, file.path(output_dir,"30_trajectory_comparison.csv"), row.names=FALSE)
saveRDS(list(robustness_id=robustness_id, reference_file=ref_file, robustness_file=rob_file,
             difference_definition="robustness minus REF", performance=performance_comparison,
             performance_full=full_comparison, annual=annual_comparison,
             trajectory=trajectory_comparison, summary_full=summary_full),
        file.path(output_dir,"30_comparison_with_REF.rds"))
cat("\nCOMPARISON COMPLETED\nMPs: 12 | Horizons: 4 | Paired trajectory-horizon rows: 48000\n")
cat("Outputs saved in:", output_dir, "\n")

