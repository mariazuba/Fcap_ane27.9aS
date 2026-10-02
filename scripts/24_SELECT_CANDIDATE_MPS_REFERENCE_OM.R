# ============================================================
# 24_SELECT_CANDIDATE_MPS_REFERENCE_OM.R
# Initial screening of candidate MPs under the reference OM
#
# Stock / analysis:
#   European anchovy ane.27.9aS (Gulf of Cadiz)
#   Shortcut MSE in FLBEIA + FLasher
#
# Role in the workflow:
#   22_EVALUATE_CANDIDATE_MPS.R
#          |
#   23_PLOT_CANDIDATE_MPS.R
#          |
#   24_SELECT_CANDIDATE_MPS_REFERENCE_OM.R
#          |
#   Stage 3 robustness testing of a documented subset of MPs
#
# Purpose:
#   Apply the initial reference-OM screening used before robustness
#   testing. The script:
#     1. reads the definitive performance metrics from script 22;
#     2. evaluates the precautionary criterion using the maximum
#        annual probability of SSB < Blim over the Full horizon;
#     3. retains all MPs satisfying max P(SSB < Blim) < 0.05;
#     4. reports stock and fishery consequences for those MPs;
#     5. produces a candidate table for subsequent documented
#        shortlist selection.
#
#   The script DOES NOT:
#     - rerun FLBEIA;
#     - recalculate MSE performance metrics;
#     - rank or score MPs;
#     - automatically select a preferred MP;
#     - define the final subset for robustness testing.
#
# Methodological basis:
#   Adapted to ane.27.9aS from the reference-OM screening logic in
#   Trochta et al. (2024), IMR-PINRO 2024-17:
#     baseline OM -> precautionary risk screening ->
#     comparison of stock/fishery consequences ->
#     subset for robustness testing.
#
# Inputs:
#   data/mse/performance/candidate_MP_performance.rds
#
# Outputs:
#   outputs/mse/candidate_MP_selection/
#     Table24a_reference_OM_precautionary_screen.csv
#     Table24b_reference_OM_eligible_MPs.csv
#     Table24c_reference_OM_consequences.csv
#     Table24d_reference_OM_by_Besc_Fcap.csv
#     reference_OM_screening_summary.rds
#
# Working Document outputs:
#   Table24a documents the precautionary screening.
#   Table24c summarises stock and fishery consequences for MPs
#   passing the reference-OM risk screen.
#
# Important:
#   Blim is the biological reference point used for the risk screen.
#   Besc is an HCR parameter and is not used as a substitute for Blim.
#
#   The 5% threshold is applied to the maximum ANNUAL probability
#   of SSB < Blim over the Full projection horizon (2025-2054).
#
#   P_ever_Blim is retained as complementary information only and
#   is not used for the initial screening.
#
#   Fishery consequences are reported after the precautionary screen
#   but are not combined into an automatic score or ranking.
# ============================================================

rm(list=ls())

# ============================================================
# 0. PACKAGES
# ============================================================

library(dplyr)
library(here)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

performance_file <- here("data","mse","performance","candidate_MP_performance.rds")
output_dir <- here("outputs","mse","candidate_MP_selection")

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

risk_threshold <- 0.05
selection_horizon <- "Full"
expected_scenarios <- 20
expected_first_year <- 2025
expected_last_year <- 2054

stopifnot(file.exists(performance_file))

# ============================================================
# 2. LOAD DEFINITIVE PERFORMANCE OUTPUTS FROM SCRIPT 22
# ============================================================

performance_outputs <- readRDS(performance_file)

scenario_grid <- performance_outputs$scenario_grid
performance <- performance_outputs$performance

stopifnot(!is.null(scenario_grid))
stopifnot(!is.null(performance))
stopifnot(nrow(scenario_grid)==expected_scenarios)
stopifnot(dplyr::n_distinct(performance$scenario_number)==expected_scenarios)

cat("\n============================================================\n")
cat("REFERENCE-OM CANDIDATE MP SCREENING\n")
cat("============================================================\n")
cat("Candidate MPs loaded:",nrow(scenario_grid),"\n")
cat("Selection horizon:",selection_horizon,"\n")
cat("Risk threshold:",risk_threshold,"\n")

# ============================================================
# 3. EXTRACT FULL-HORIZON PERFORMANCE
# ============================================================

performance_full <- performance %>%
  filter(horizon==selection_horizon) %>%
  arrange(Besc,Fcap)

stopifnot(nrow(performance_full)==expected_scenarios)
stopifnot(dplyr::n_distinct(performance_full$scenario_number)==expected_scenarios)
stopifnot(all(performance_full$first_year==expected_first_year))
stopifnot(all(performance_full$last_year==expected_last_year))
stopifnot(!anyDuplicated(performance_full$scenario_number))
stopifnot(all(is.finite(performance_full$max_P_Blim)))
stopifnot(all(performance_full$max_P_Blim>=0))
stopifnot(all(performance_full$max_P_Blim<=1))

# ============================================================
# 4. PRECAUTIONARY SCREEN UNDER THE REFERENCE OM
# ============================================================
#
# Primary screening criterion:
#   max annual P(SSB < Blim) < 0.05
#
# n_years_P_Blim_gt_005 and P_ever_Blim are retained as supporting
# diagnostics. They do not replace the primary criterion.
# ============================================================

table_precautionary <- performance_full %>%
  transmute(
    scenario_number=scenario_number,
    scenario_id=scenario_id,
    Fcap=Fcap,
    Besc=Besc,
    Blim=Blim,
    Bpa=Bpa,
    horizon=horizon,
    first_year=first_year,
    last_year=last_year,
    max_P_Blim=max_P_Blim,
    n_years_P_Blim_gt_005=n_years_P_Blim_gt_5,
    P_ever_Blim=P_ever_Blim,
    pass_reference_OM_risk=max_P_Blim<risk_threshold
  ) %>%
  arrange(Besc,Fcap)

stopifnot(nrow(table_precautionary)==expected_scenarios)
stopifnot(!any(is.na(table_precautionary$pass_reference_OM_risk)))

n_pass <- sum(table_precautionary$pass_reference_OM_risk)
n_fail <- sum(!table_precautionary$pass_reference_OM_risk)

cat("\nPrecautionary screening\n")
cat("-----------------------\n")
cat("MPs passing max annual P(SSB < Blim) < 0.05:",n_pass,"\n")
cat("MPs failing max annual P(SSB < Blim) < 0.05:",n_fail,"\n")
cat("Maximum max_P_Blim across candidate MPs:",max(table_precautionary$max_P_Blim),"\n")

# ============================================================
# 5. MPs ELIGIBLE FOR CONSEQUENCE COMPARISON
# ============================================================
#
# Passing this screen means that the MP proceeds to comparison of
# stock and fishery consequences under the reference OM.
#
# It does NOT mean that the MP has been selected for robustness
# testing or recommended as a management procedure.
# ============================================================

eligible_ids <- table_precautionary %>%
  filter(pass_reference_OM_risk) %>%
  pull(scenario_number)

table_eligible <- table_precautionary %>%
  filter(pass_reference_OM_risk) %>%
  arrange(Besc,Fcap)

stopifnot(nrow(table_eligible)==n_pass)
stopifnot(all(table_eligible$scenario_number %in% eligible_ids))

# ============================================================
# 6. STOCK AND FISHERY CONSEQUENCES FOR ELIGIBLE MPs
# ============================================================
#
# These metrics are read directly from script 22.
#
# Risk:
#   max_P_Blim
#   n_years_P_Blim_gt_005
#   P_ever_Blim
#
# Stock:
#   min_SSB
#   median_mean_SSB_open
#
# Fishery:
#   mean_catch
#   mean_catch_open
#   mean_P_closed
#   mean_years_closed
#   median_years_closed
#   mean_IAV
#   median_IAV
#
# Conditional metrics retain their explicit denominators.
# ============================================================

table_consequences <- performance_full %>%
  filter(scenario_number %in% eligible_ids) %>%
  transmute(
    scenario_number=scenario_number,
    scenario_id=scenario_id,
    Fcap=Fcap,
    Besc=Besc,
    max_P_Blim=max_P_Blim,
    n_years_P_Blim_gt_005=n_years_P_Blim_gt_5,
    P_ever_Blim=P_ever_Blim,
    min_SSB=min_SSB,
    n_trajectories_SSB_open=n_trajectories_SSB_open,
    P_trajectories_SSB_open=P_trajectories_SSB_open,
    median_mean_SSB_open=median_mean_SSB_open,
    mean_catch=mean_catch,
    n_trajectories_open=n_trajectories_open,
    P_trajectories_open=P_trajectories_open,
    mean_catch_open=mean_catch_open,
    median_mean_catch_open=median_mean_catch_open,
    mean_P_closed=mean_P_closed,
    mean_years_closed=mean_years_closed,
    median_years_closed=median_years_closed,
    mean_IAV=mean_IAV,
    median_IAV=median_IAV
  ) %>%
  arrange(Besc,Fcap)

stopifnot(nrow(table_consequences)==n_pass)
stopifnot(all(table_consequences$mean_P_closed>=0))
stopifnot(all(table_consequences$mean_P_closed<=1))
stopifnot(all(table_consequences$P_trajectories_open>=0))
stopifnot(all(table_consequences$P_trajectories_open<=1))
stopifnot(all(table_consequences$P_trajectories_SSB_open>=0))
stopifnot(all(table_consequences$P_trajectories_SSB_open<=1))

# ============================================================
# 7. STRUCTURED COMPARISON BY Besc AND Fcap
# ============================================================
#
# This table preserves the complete eligible candidate set and
# facilitates comparison within and among Besc values.
#
# No score, ranking, Pareto filter or automatic shortlist is applied.
# ============================================================

table_by_Besc_Fcap <- table_consequences %>%
  select(
    scenario_number,
    Fcap,
    Besc,
    max_P_Blim,
    P_ever_Blim,
    mean_catch,
    mean_catch_open,
    mean_P_closed,
    mean_years_closed,
    median_mean_SSB_open,
    median_IAV
  ) %>%
  arrange(Besc,Fcap)

# ============================================================
# 8. SCREENING SUMMARY
# ============================================================

screening_summary <- list(
  risk_threshold=risk_threshold,
  selection_horizon=selection_horizon,
  first_year=expected_first_year,
  last_year=expected_last_year,
  n_candidate_MPs=expected_scenarios,
  n_pass_reference_OM_risk=n_pass,
  n_fail_reference_OM_risk=n_fail,
  eligible_scenario_numbers=eligible_ids,
  precautionary_screen=table_precautionary,
  eligible_MPs=table_eligible,
  consequences=table_consequences,
  comparison_by_Besc_Fcap=table_by_Besc_Fcap
)

# ============================================================
# 9. SAVE OUTPUTS
# ============================================================

write.csv(
  table_precautionary,
  file.path(output_dir,"Table24a_reference_OM_precautionary_screen.csv"),
  row.names=FALSE
)

write.csv(
  table_eligible,
  file.path(output_dir,"Table24b_reference_OM_eligible_MPs.csv"),
  row.names=FALSE
)

write.csv(
  table_consequences,
  file.path(output_dir,"Table24c_reference_OM_consequences.csv"),
  row.names=FALSE
)

write.csv(
  table_by_Besc_Fcap,
  file.path(output_dir,"Table24d_reference_OM_by_Besc_Fcap.csv"),
  row.names=FALSE
)

saveRDS(
  screening_summary,
  file.path(output_dir,"reference_OM_screening_summary.rds")
)

# ============================================================
# 10. FINAL VALIDATION CHECKS
# ============================================================

stopifnot(file.exists(file.path(output_dir,"Table24a_reference_OM_precautionary_screen.csv")))
stopifnot(file.exists(file.path(output_dir,"Table24b_reference_OM_eligible_MPs.csv")))
stopifnot(file.exists(file.path(output_dir,"Table24c_reference_OM_consequences.csv")))
stopifnot(file.exists(file.path(output_dir,"Table24d_reference_OM_by_Besc_Fcap.csv")))
stopifnot(file.exists(file.path(output_dir,"reference_OM_screening_summary.rds")))
stopifnot(nrow(table_precautionary)==n_pass+n_fail)
stopifnot(nrow(table_eligible)==nrow(table_consequences))
stopifnot(nrow(table_consequences)==nrow(table_by_Besc_Fcap))

cat("\n============================================================\n")
cat("REFERENCE-OM SCREENING COMPLETED\n")
cat("============================================================\n")
cat("Candidate MPs evaluated:",expected_scenarios,"\n")
cat("MPs passing precautionary screen:",n_pass,"\n")
cat("MPs failing precautionary screen:",n_fail,"\n")
cat("Output directory:",output_dir,"\n")
cat("\nNo automatic robustness shortlist has been selected.\n")
cat("The next step is to review the consequences of the eligible MPs\n")
cat("and document the subset to be carried forward to Stage 3 robustness.\n")

