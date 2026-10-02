# ============================================================
# 22b_DIAGNOSE_BLIM_METHOD_SENSITIVITY.R
# Diagnostic sensitivity of biological risk to Blim derivation
#
# Stock:
#   European anchovy ane.27.9aS (Gulf of Cadiz)
#   Shortcut MSE in FLBEIA + FLasher conditioned from SS3
#
# Role in workflow:
#   22_EVALUATE_CANDIDATE_MPS.R
#     -> 22b_DIAGNOSE_BLIM_METHOD_SENSITIVITY.R
#
# Purpose:
#   Diagnose whether the very low estimated probability of
#   SSB < Blim under the candidate MPs is associated with the
#   position of the adopted Blim relative to the lower tail of
#   the simulated SSB distribution.
#
#   The analysis recalculates diagnostic biomass thresholds from
#   the Blim derivation method considered during the benchmark,
#   rather than using the historical rounded threshold values.
#
#   This is a diagnostic sensitivity analysis only.
#   It does NOT redefine the adopted biological reference points
#   and does NOT modify the candidate MP performance evaluation.
#
# Final benchmark reference points:
#   Bpa = Bloss = 6561 t
#
#   Adopted Blim:
#
#     Blim = Bpa * exp(-1.645 * sigmaB)
#
#   with sigmaB = 0.20, giving approximately 4721 t.
#
# Diagnostic formulations:
#   sigmaB = 0.30
#   sigmaB = 0.20  [REFERENCE]
#   sigmaB = 0.10
#   Bpa = Bloss    [upper diagnostic threshold only]
#
# The 0.2 * B0 formulation is not included because the purpose
# here is specifically to diagnose the region around and above
# the adopted Blim.
#
# Main diagnostics:
#   1. Annual P(SSB < threshold)
#   2. Annual lower-tail SSB quantiles:
#        q01, q05, q10, median
#   3. Maximum annual probability by evaluation horizon
#   4. Number of years with annual P > 0.05
#   5. Distance of the lower SSB tail from each threshold
#
# IMPORTANT:
#   - No new MSE simulations are performed.
#   - Existing MP output blocks are read directly.
#   - Biological risk remains evaluated annually.
#   - Bpa/Bloss is used only as an upper diagnostic threshold.
# ============================================================


rm(list = ls())


library(FLCore)
library(FLBEIA)
library(dplyr)
library(tidyr)
library(purrr)
library(here)


# ============================================================
# DIRECTORIES AND SETTINGS
# ============================================================

scenario_file <- here("data","mse","scenarios","candidate_MP_grid.csv")
scenario_root <- here("data","mse","MP_runs")
diagnostic_dir <- here("data","mse","diagnostics","Blim_method_sensitivity")

dir.create(diagnostic_dir,recursive=TRUE,showWarnings=FALSE)

block_ids <- 1:10
risk_threshold <- 0.05

horizons <- list(
  Short=2025:2034,
  Medium=2035:2044,
  Long=2045:2054,
  Full=2025:2054)


# ============================================================
# FINAL BENCHMARK REFERENCE POINTS
# ============================================================

Bpa_final <- 6561
Bloss_final <- Bpa_final
Blim_final <- 4721
z95 <- 1.645

# ============================================================
# CALCULATE DIAGNOSTIC THRESHOLDS
# ============================================================

Blim_thresholds <- tibble(threshold_id=c(
    "Blim_sigma030",
    "Blim_reference",
    "Blim_sigma010",
    "Bpa_Bloss"),
  sigmaB=c(
    0.30,
    0.20,
    0.10,
    NA_real_),
  threshold=c(
    Bpa_final * exp(-z95 * 0.30),
    Bpa_final * exp(-z95 * 0.20),
    Bpa_final * exp(-z95 * 0.10),
    Bloss_final),
  role=c(
    "diagnostic",
    "reference",
    "diagnostic",
    "upper_diagnostic"))


# ============================================================
# CHECK CALCULATED THRESHOLDS
# ============================================================

cat("\n============================================================\n")
cat("DIAGNOSTIC BIOMASS THRESHOLDS\n")
cat("============================================================\n")

print(Blim_thresholds)

Blim_reference_calc <- Blim_thresholds %>%
  filter(threshold_id=="Blim_reference") %>%
  pull(threshold)

cat("\nCalculated reference Blim:",round(Blim_reference_calc,2),"t\n")
cat("Final benchmark Blim:",Blim_final,"t\n")

stopifnot(abs(Blim_reference_calc-Blim_final)<10)

# ============================================================
# LOAD CANDIDATE MP GRID
# ============================================================

stopifnot(file.exists(scenario_file))
scenario_grid <- read.csv(scenario_file,stringsAsFactors=FALSE)

required_columns <- c(
  "scenario_number",
  "scenario_id",
  "Fcap",
  "Besc",
  "Blim",
  "Bpa",
  "first_projection_year",
  "last_projection_year")

stopifnot(all(required_columns %in% names(scenario_grid)))

scenario_grid <- scenario_grid %>%
  arrange(scenario_number)

n_scenarios <- nrow(scenario_grid)

stopifnot(n_distinct(scenario_grid$scenario_number)==n_scenarios)
stopifnot(n_distinct(scenario_grid$scenario_id)==n_scenarios)

cat("\n============================================================\n")
cat("CANDIDATE MP GRID\n")
cat("============================================================\n")

cat("Candidate MPs:",n_scenarios,"\n")
print(scenario_grid %>%
    dplyr::select(scenario_number,scenario_id,Fcap,Besc,Blim,Bpa))

# ============================================================
# CHECK FINAL REFERENCE POINTS IN SCENARIO GRID
# ============================================================

stopifnot(all(abs(scenario_grid$Blim-Blim_final)<1))
stopifnot(all(abs(scenario_grid$Bpa-Bpa_final)<1))

cat("\nReference-point checks passed:","Blim =",Blim_final,"t | Bpa =",Bpa_final,"t\n")

# ============================================================
# HELPERS
# ============================================================

add_global_map <- function(dat,map_b) {
  
  out <- dat %>%
    mutate(iter_local=as.integer(as.character(iter))) %>%
    left_join(map_b %>% dplyr::select(iter_local,iter1000,om,replicate),by="iter_local")
  
  stopifnot(!anyNA(out$iter1000))
  stopifnot(!anyNA(out$om))
  stopifnot(!anyNA(out$replicate))
  
  out
}


standardise_dynamic <- function(x,map_b,proj_years) {
  
  as.data.frame(x) %>%
    rename(value=data) %>%
    mutate(year=as.integer(as.character(year))) %>%
    add_global_map(map_b) %>%
    filter(year %in% proj_years) %>%
    transmute(year,iter1000,om,replicate,value)
}


# ============================================================
# EXTRACT SSB FOR ONE CANDIDATE MP
# ============================================================

extract_candidate_SSB <- function(scenario_row) {
  
  scenario_number <- scenario_row$scenario_number[[1]]
  scenario_id <- scenario_row$scenario_id[[1]]
  Fcap <- scenario_row$Fcap[[1]]
  Besc <- scenario_row$Besc[[1]]
  first_proj <- scenario_row$first_projection_year[[1]]
  last_proj <- scenario_row$last_projection_year[[1]]
  proj_years <- first_proj:last_proj
  
  cat("\n============================================================\n")
  cat("SCENARIO:",scenario_number,"\n")
  cat("ID:",scenario_id,"\n")
  cat("Fcap:",Fcap,"| Besc:",Besc,"\n")
  cat("============================================================\n")
  
  
  # ==========================================================
  # FILES
  # ==========================================================
  
  mp_dir <- file.path(scenario_root,scenario_id)
  mp_files <- file.path(mp_dir,sprintf("block_%02d.rds",block_ids))
  map_files <- file.path(mp_dir,sprintf("block_%02d_map.csv",block_ids))
  
  if(!dir.exists(mp_dir)) {stop("Scenario directory not found: ",mp_dir)}
  if(!all(file.exists(mp_files))) {stop("Missing RDS blocks for scenario ",scenario_number)}
  if(!all(file.exists(map_files))) {stop("Missing map files for scenario ",scenario_number)}
  
  # ==========================================================
  # CHECK COMPLETE ITERATION MAP
  # ==========================================================
  
  iteration_maps <- map_dfr(block_ids, ~read.csv(map_files[.x]) %>%mutate(block=.x))
  
  stopifnot(nrow(iteration_maps)==1000)
  stopifnot(n_distinct(iteration_maps$iter1000)==1000)
  stopifnot(n_distinct(iteration_maps$om)==100)
  stopifnot(all(table(iteration_maps$om)==10))
  stopifnot(all(sort(unique(iteration_maps$replicate))==1:10))
  
  # ==========================================================
  # EXTRACT SSB FROM ONE BLOCK
  # ==========================================================
  
  extract_block_SSB <- function(block_id) {
    
    OM_b <- readRDS(mp_files[block_id])
    map_b <- read.csv(map_files[block_id])
    biol <- OM_b$biols$ANE
    Mat <- predict(biol@mat)
    SSB_q2 <- quantSums(biol@n[,,,"2",,] *biol@wt[,,,"2",,] *Mat[,,,"2",,])
    B_b <- standardise_dynamic(SSB_q2,map_b,proj_years)
    
    B_b
  }
  
  
  # ==========================================================
  # EXTRACT COMPLETE SSB ENSEMBLE
  # ==========================================================
  
  B_MP <- map_dfr(block_ids,extract_block_SSB)
  
  # ==========================================================
  # VALIDATE EXTRACTED SSB
  # ==========================================================
  
  expected_rows <- length(proj_years)*1000
  
  stopifnot(nrow(B_MP)==expected_rows)
  stopifnot(n_distinct(B_MP$iter1000)==1000)
  stopifnot(n_distinct(B_MP$om)==100)
  stopifnot(all(sort(unique(B_MP$replicate))==1:10))
  stopifnot(min(B_MP$year)==first_proj)
  stopifnot(max(B_MP$year)==last_proj)
  stopifnot(!anyNA(B_MP$value))
  stopifnot(all(is.finite(B_MP$value)))
  
  stopifnot(all(B_MP$value>=0))
  
  cat("SSB extraction passed:",nrow(B_MP),"rows |",n_distinct(B_MP$iter1000),"trajectories\n")
  
  # ==========================================================
  # ADD SCENARIO INFORMATION
  # ==========================================================
  
  B_MP <- B_MP %>%
    mutate(
      scenario_number=scenario_number,
      scenario_id=scenario_id,
      Fcap=Fcap,
      Besc=Besc,
      .before=1)
  
  B_MP
}


# ============================================================
# EXTRACT SSB FOR ALL CANDIDATE MPs
# ============================================================

SSB_all <- map_dfr(seq_len(nrow(scenario_grid)),~extract_candidate_SSB(scenario_grid[.x,,drop=FALSE]))

cat("\n============================================================\n")
cat("ALL SSB TRAJECTORIES EXTRACTED\n")
cat("============================================================\n")

cat("Scenarios:",n_distinct(SSB_all$scenario_number),"\n")
cat("Rows:",nrow(SSB_all),"\n")

# ============================================================
# STRUCTURAL CHECKS
# ============================================================

stopifnot(n_distinct(SSB_all$scenario_number)==n_scenarios)
stopifnot(all(table(SSB_all$scenario_number)==30000))
stopifnot(!anyNA(SSB_all$value))
stopifnot(all(is.finite(SSB_all$value)))
stopifnot(all(SSB_all$value>=0))

cat("Structural checks passed\n")


# ============================================================
# ANNUAL LOWER-TAIL SSB DIAGNOSTICS
# ============================================================

SSB_annual_tail <- SSB_all %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    year) %>%
  summarise(
    n=dplyr::n(),
    SSB_min=min(value),
    SSB_q01=quantile(value,0.01),
    SSB_q05=quantile(value,0.05),
    SSB_q10=quantile(value,0.10),
    SSB_median=median(value),
    SSB_mean=mean(value),
    .groups="drop")


stopifnot(all(SSB_annual_tail$n==1000))


# ============================================================
# ANNUAL PROBABILITY BELOW EACH DIAGNOSTIC THRESHOLD
# ============================================================

SSB_annual_threshold <- SSB_all %>%
  tidyr::crossing(
    Blim_thresholds) %>%
  mutate(
    below_threshold=value<threshold) %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    year,
    threshold_id,
    sigmaB,
    threshold,
    role) %>%
  summarise(
    n=dplyr::n(),
    n_below=sum(below_threshold),
    P_below=mean(below_threshold),
    .groups="drop")


stopifnot(all(SSB_annual_threshold$n==1000))
stopifnot(all(SSB_annual_threshold$P_below>=0 &SSB_annual_threshold$P_below<=1))


# ============================================================
# CHECK MONOTONICITY OF THRESHOLD RISK
# ============================================================

threshold_order <- Blim_thresholds %>%
  arrange(threshold) %>%
  pull(threshold_id)

risk_monotonicity <- SSB_annual_threshold %>%
  dplyr::select(
    scenario_number,
    year,
    threshold_id,
    P_below) %>%
  mutate(
    threshold_id=factor(
      threshold_id,
      levels=threshold_order)) %>%
  arrange(
    scenario_number,
    year,
    threshold_id) %>%
  group_by(
    scenario_number,
    year) %>%
  summarise(
    monotonic=all(diff(P_below)>=-1e-12),
    .groups="drop")


stopifnot(all(risk_monotonicity$monotonic))

cat("Threshold-risk monotonicity check passed\n")


# ============================================================
# DEFINE EVALUATION HORIZONS
# ============================================================

horizon_table <- imap_dfr(horizons,~tibble(horizon=.y,year=.x))

# ============================================================
# ADD HORIZONS TO ANNUAL DIAGNOSTICS
# ============================================================

SSB_annual_threshold_h <- SSB_annual_threshold %>%inner_join(horizon_table,by="year")
SSB_annual_tail_h <- SSB_annual_tail %>%inner_join(horizon_table,by="year")

# ============================================================
# THRESHOLD SENSITIVITY BY HORIZON
# ============================================================

Blim_sensitivity <- SSB_annual_threshold_h %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    horizon,
    threshold_id,
    sigmaB,
    threshold,
    role) %>%
  summarise(
    max_annual_P_below=max(P_below),
    mean_annual_P_below=mean(P_below),
    n_years_P_gt_5=sum(P_below>risk_threshold),
    year_max_P=year[which.max(P_below)],
    .groups="drop")


# ============================================================
# LOWER-TAIL SSB BY HORIZON
# ============================================================

SSB_tail_horizon <- SSB_annual_tail_h %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    horizon) %>%
  summarise(
    min_SSB=min(SSB_min),
    min_annual_q01=min(SSB_q01),
    min_annual_q05=min(SSB_q05),
    min_annual_q10=min(SSB_q10),
    min_annual_median=min(SSB_median),
    .groups="drop")


# ============================================================
# DISTANCE OF LOWER SSB TAIL FROM EACH THRESHOLD
# ============================================================

SSB_tail_distance <- SSB_annual_tail_h %>%
  tidyr::crossing(
    Blim_thresholds) %>%
  mutate(
    q01_minus_threshold=SSB_q01-threshold,
    q05_minus_threshold=SSB_q05-threshold,
    q10_minus_threshold=SSB_q10-threshold) %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    horizon,
    threshold_id,
    sigmaB,
    threshold,
    role) %>%
  summarise(
    min_q01_minus_threshold=min(q01_minus_threshold),
    min_q05_minus_threshold=min(q05_minus_threshold),
    min_q10_minus_threshold=min(q10_minus_threshold),
    .groups="drop")


# ============================================================
# COMBINE HORIZON DIAGNOSTICS
# ============================================================

Blim_diagnostic_by_horizon <- Blim_sensitivity %>%
  left_join(SSB_tail_distance,
    by=c(
      "scenario_number",
      "scenario_id",
      "Fcap",
      "Besc",
      "horizon",
      "threshold_id",
      "sigmaB",
      "threshold",
      "role")) %>%
  left_join(SSB_tail_horizon,
    by=c(
      "scenario_number",
      "scenario_id",
      "Fcap",
      "Besc",
      "horizon")) %>%
  arrange(scenario_number,horizon,threshold)


# ============================================================
# FULL-HORIZON DIAGNOSTIC
# ============================================================

Blim_diagnostic_full <- Blim_diagnostic_by_horizon %>% filter(horizon=="Full") %>% arrange(Fcap,Besc,threshold)

# ============================================================
# REFERENCE BLIM DIAGNOSTIC
# ============================================================

reference_diagnostic <- Blim_diagnostic_full %>%filter(threshold_id=="Blim_reference")

cat("\n============================================================\n")
cat("REFERENCE BLIM DIAGNOSTIC\n")
cat("============================================================\n")

print(reference_diagnostic %>%
    dplyr::select(
      scenario_number,
      Fcap,
      Besc,
      threshold,
      max_annual_P_below,
      n_years_P_gt_5,
      min_SSB,
      min_annual_q01,
      min_annual_q05,
      min_q01_minus_threshold,
      min_q05_minus_threshold))


# ============================================================
# FULL-HORIZON COMPARISON ACROSS THRESHOLDS
# ============================================================

comparison_full <- Blim_diagnostic_full %>%
  dplyr::select(
    scenario_number,
    Fcap,
    Besc,
    threshold_id,
    max_annual_P_below) %>%
  pivot_wider(
    names_from=threshold_id,
    values_from=max_annual_P_below,
    names_prefix="maxP_") %>%
  arrange(Fcap,Besc)

stopifnot(nrow(comparison_full)==n_scenarios)

cat("\n============================================================\n")
cat("MAXIMUM ANNUAL RISK ACROSS DIAGNOSTIC THRESHOLDS\n")
cat("============================================================\n")

print(comparison_full,n=Inf)

# ============================================================
# LOWER-TAIL SUMMARY - FULL HORIZON
# ============================================================

tail_summary_full <- SSB_tail_horizon %>%filter(horizon=="Full") %>%arrange(Fcap,Besc)

cat("\n============================================================\n")
cat("LOWER SSB TAIL - FULL HORIZON\n")
cat("============================================================\n")

print(tail_summary_full,n=Inf)

# ============================================================
# GLOBAL DIAGNOSTIC SUMMARY
# ============================================================

global_threshold_summary <- Blim_diagnostic_full %>%
  group_by(threshold_id,sigmaB,threshold,role) %>%
  summarise(
    max_risk_across_MPs=max(max_annual_P_below),
    n_MPs_with_any_risk=sum(max_annual_P_below>0),
    n_MPs_above_5pct=sum(max_annual_P_below>risk_threshold),
    min_q01_distance=min(min_q01_minus_threshold),
    min_q05_distance=min(min_q05_minus_threshold),
    .groups="drop") %>%
  arrange(threshold)


cat("\n============================================================\n")
cat("GLOBAL THRESHOLD DIAGNOSTIC\n")
cat("============================================================\n")

print(global_threshold_summary,n=Inf)

# ============================================================
# FINAL VALIDATION
# ============================================================

stopifnot(n_distinct(Blim_diagnostic_full$scenario_number)==n_scenarios)
stopifnot(all(Blim_diagnostic_full$max_annual_P_below>=0 &Blim_diagnostic_full$max_annual_P_below<=1))
stopifnot(all(Blim_diagnostic_full$n_years_P_gt_5>=0))
stopifnot(all(is.finite(Blim_diagnostic_full$min_SSB)))
stopifnot(all(is.finite(Blim_diagnostic_full$min_annual_q01)))
stopifnot(all(is.finite(Blim_diagnostic_full$min_annual_q05)))

cat("\nFINAL VALIDATION PASSED\n")


# ============================================================
# SAVE OUTPUTS
# ============================================================

diagnostic_outputs <- list(
  reference_points=Blim_thresholds,
  annual_lower_tail=SSB_annual_tail,
  annual_threshold_risk=SSB_annual_threshold,
  threshold_sensitivity_by_horizon=Blim_sensitivity,
  lower_tail_by_horizon=SSB_tail_horizon,
  tail_distance=SSB_tail_distance,
  diagnostic_by_horizon=Blim_diagnostic_by_horizon,
  diagnostic_full=Blim_diagnostic_full,
  comparison_full=comparison_full,
  tail_summary_full=tail_summary_full,
  global_threshold_summary=global_threshold_summary)


saveRDS(diagnostic_outputs,file.path(diagnostic_dir,"Blim_method_sensitivity.rds"))

write.csv(Blim_thresholds,file.path(diagnostic_dir,"Blim_diagnostic_thresholds.csv"),row.names=FALSE)
write.csv(SSB_annual_tail,file.path(diagnostic_dir,"SSB_annual_lower_tail.csv"),row.names=FALSE)
write.csv(SSB_annual_threshold,file.path(diagnostic_dir,"SSB_annual_threshold_risk.csv"), row.names=FALSE)
write.csv(Blim_diagnostic_by_horizon,file.path(diagnostic_dir,"Blim_diagnostic_by_horizon.csv"),row.names=FALSE)
write.csv(Blim_diagnostic_full,file.path(diagnostic_dir,"Blim_diagnostic_full.csv"),row.names=FALSE)
write.csv(comparison_full,file.path(diagnostic_dir,"Blim_sensitivity_comparison_full.csv"),row.names=FALSE)
write.csv(tail_summary_full,file.path(diagnostic_dir,"SSB_lower_tail_full.csv"),row.names=FALSE)
write.csv(global_threshold_summary,file.path(diagnostic_dir,"Blim_global_threshold_summary.csv"),row.names=FALSE)


# ============================================================
# END
# ============================================================

cat("\n============================================================\n")
cat("BLIM METHOD SENSITIVITY DIAGNOSTIC COMPLETED\n")
cat("============================================================\n")

cat("Candidate MPs evaluated:",n_scenarios,"\n")
cat("Reference Blim:",round(Blim_reference_calc,2),"t\n")
cat("Outputs saved in:",diagnostic_dir,"\n")
