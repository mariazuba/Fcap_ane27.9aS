# ============================================================
# 24c_DIAGNOSE_REFINED_MP_TRADEOFFS.R
# Marginal trade-off diagnostics for the definitive 28-MP grid
#
# Stock / analysis:
#   European anchovy ane.27.9aS (Gulf of Cadiz)
#   Shortcut MSE in FLBEIA + FLasher conditioned from SS3
#
# Role in the workflow:
#   22b_EVALUATE_REFINED_CANDIDATE_MPS.R
#     -> 24c_DIAGNOSE_REFINED_MP_TRADEOFFS.R
#     -> definition of candidate MPs for robustness testing
#
# Purpose:
#   Quantify marginal changes in validated Reference-OM performance
#   across the definitive 7 x 4 Fcap-Besc grid.
#
#   A. Within each Besc, compare consecutive Fcap values.
#   B. Within each Fcap, compare consecutive Besc values.
#   C. Express changes both as raw increments and per unit change
#      in Fcap or per 1000 t change in Besc.
#
#   This script DOES NOT:
#     - rerun FLBEIA;
#     - recalculate MSE performance metrics;
#     - introduce an acceptance threshold;
#     - rank MPs;
#     - select or reject MPs automatically.
#
# Input:
#   data/mse/performance_refined/
#     refined_candidate_MP_performance_full.csv
#
# Outputs:
#   outputs/mse/candidate_MP_selection/refined_diagnostics/
#     Table24c_01_reference_metrics.csv
#     Table24c_02_marginal_Fcap.csv
#     Table24c_03_marginal_Besc.csv
#     Table24c_04_Fcap_150_180.csv
#     Table24c_05_Besc_6561_8000.csv
#     refined_MP_tradeoffs_reference_OM.rds
#
# Interpretation:
#   Positive delta values mean that the metric increased from the
#   immediately preceding candidate level.
#
#   No delta is automatically interpreted as beneficial or adverse.
#   In particular:
#     - higher mean_catch = greater average yield;
#     - higher mean_P_closed / mean_years_closed = greater closure burden;
#     - higher max_P_Blim = greater annual biological risk;
#     - higher P_ever_Blim = greater probability of at least one Blim
#       crossing over the evaluated horizon;
#     - IAV must be interpreted jointly with closure metrics because
#       consecutive zero-catch years can reduce IAV.
#
# Important:
#   The 0.05 annual Blim-risk threshold is reported from script 22b
#   but is not redefined here.
# ============================================================

rm(list=ls())

# ============================================================
# 0. PACKAGES
# ============================================================

library(dplyr)
library(here)

# ============================================================
# 1. PATHS AND EXPECTED DESIGN
# ============================================================

input_file <- here("data","mse","performance_refined","refined_candidate_MP_performance_full.csv")
output_dir <- here("outputs","mse","candidate_MP_selection","refined_diagnostics")

expected_Fcap <- c(1.00,1.25,1.50,1.60,1.70,1.75,1.80)
expected_Besc <- c(5000,6561,7250,8000)
expected_scenarios <- length(expected_Fcap)*length(expected_Besc)
risk_threshold <- 0.05

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

stopifnot(file.exists(input_file))

# ============================================================
# 2. LOAD VALIDATED REFERENCE-OM PERFORMANCE
# ============================================================

performance_ref <- read.csv(input_file,stringsAsFactors=FALSE)

required_columns <- c(
  "scenario_number",
  "scenario_id",
  "Fcap",
  "Besc",
  "max_P_Blim",
  "n_years_P_Blim_gt_5",
  "P_ever_Blim",
  "max_P_Bpa",
  "P_ever_Bpa",
  "min_SSB",
  "mean_catch",
  "mean_catch_open",
  "mean_P_closed",
  "mean_years_closed",
  "median_years_closed",
  "median_mean_SSB_open",
  "mean_IAV",
  "median_IAV"
)

missing_columns <- setdiff(required_columns,names(performance_ref))

stopifnot(length(missing_columns)==0)
stopifnot(nrow(performance_ref)==expected_scenarios)
stopifnot(n_distinct(performance_ref$scenario_number)==expected_scenarios)
stopifnot(n_distinct(performance_ref$scenario_id)==expected_scenarios)
stopifnot(setequal(sort(unique(performance_ref$Fcap)),expected_Fcap))
stopifnot(setequal(sort(unique(performance_ref$Besc)),expected_Besc))
stopifnot(!anyDuplicated(performance_ref[c("Fcap","Besc")]))

performance_ref <- performance_ref %>%
  arrange(Besc,Fcap)

cat("\n============================================================\n")
cat("REFINED MP TRADE-OFF DIAGNOSTICS\n")
cat("============================================================\n")
cat("Candidate MPs:",nrow(performance_ref),"\n")
cat("Fcap values:",paste(expected_Fcap,collapse=", "),"\n")
cat("Besc values:",paste(expected_Besc,collapse=", "),"\n")
cat("Annual Blim-risk threshold:",risk_threshold,"\n")

# ============================================================
# 3. PRECAUTIONARITY SUMMARY
# ============================================================

precautionarity_summary <- performance_ref %>%
  summarise(
    n_MPs=dplyr::n(),
    max_observed_P_Blim=max(max_P_Blim),
    n_MPs_max_P_Blim_gt_005=sum(max_P_Blim>risk_threshold),
    n_MPs_with_years_gt_005=sum(n_years_P_Blim_gt_5>0),
    max_P_ever_Blim=max(P_ever_Blim),
    min_SSB_over_grid=min(min_SSB)
  )

stopifnot(precautionarity_summary$n_MPs==expected_scenarios)

# ============================================================
# 4. MARGINAL CHANGES WITH INCREASING Fcap
# ============================================================
#
# Within each Besc, each Fcap is compared with the immediately
# preceding Fcap value. Because the refined grid has unequal Fcap
# increments, both raw deltas and deltas per 0.10 Fcap are reported.
# ============================================================

table_marginal_Fcap <- performance_ref %>%
  group_by(Besc) %>%
  arrange(Fcap,.by_group=TRUE) %>%
  mutate(
    Fcap_previous=lag(Fcap),
    delta_Fcap=Fcap-Fcap_previous,
    delta_mean_catch=mean_catch-lag(mean_catch),
    delta_mean_P_closed=mean_P_closed-lag(mean_P_closed),
    delta_mean_years_closed=mean_years_closed-lag(mean_years_closed),
    delta_median_IAV=median_IAV-lag(median_IAV),
    delta_median_mean_SSB_open=median_mean_SSB_open-lag(median_mean_SSB_open),
    delta_max_P_Blim=max_P_Blim-lag(max_P_Blim),
    delta_P_ever_Blim=P_ever_Blim-lag(P_ever_Blim),
    delta_mean_catch_per_0.10_Fcap=delta_mean_catch/delta_Fcap*0.10,
    delta_mean_P_closed_per_0.10_Fcap=delta_mean_P_closed/delta_Fcap*0.10,
    delta_mean_years_closed_per_0.10_Fcap=delta_mean_years_closed/delta_Fcap*0.10,
    delta_max_P_Blim_per_0.10_Fcap=delta_max_P_Blim/delta_Fcap*0.10
  ) %>%
  filter(!is.na(Fcap_previous)) %>%
  ungroup() %>%
  select(
    scenario_number,
    scenario_id,
    Fcap_previous,
    Fcap,
    delta_Fcap,
    Besc,
    mean_catch,
    mean_P_closed,
    mean_years_closed,
    median_IAV,
    max_P_Blim,
    P_ever_Blim,
    delta_mean_catch,
    delta_mean_P_closed,
    delta_mean_years_closed,
    delta_median_IAV,
    delta_median_mean_SSB_open,
    delta_max_P_Blim,
    delta_P_ever_Blim,
    delta_mean_catch_per_0.10_Fcap,
    delta_mean_P_closed_per_0.10_Fcap,
    delta_mean_years_closed_per_0.10_Fcap,
    delta_max_P_Blim_per_0.10_Fcap
  ) %>%
  arrange(Besc,Fcap)

stopifnot(nrow(table_marginal_Fcap)==length(expected_Besc)*(length(expected_Fcap)-1))
stopifnot(all(table_marginal_Fcap$delta_Fcap>0))

# ============================================================
# 5. Fcap MARGINAL YIELD RELATIVE TO CLOSURE BURDEN
# ============================================================
#
# These are descriptive ratios, not scores or selection criteria.
# NA is retained when the denominator does not increase.
# ============================================================

table_marginal_Fcap <- table_marginal_Fcap %>%
  mutate(
    catch_gain_per_0.01_Pclosed=if_else(
      delta_mean_P_closed>0,
      delta_mean_catch/(delta_mean_P_closed/0.01),
      NA_real_
    ),
    catch_gain_per_extra_closed_year=if_else(
      delta_mean_years_closed>0,
      delta_mean_catch/delta_mean_years_closed,
      NA_real_
    )
  )

# ============================================================
# 6. MARGINAL CHANGES WITH INCREASING Besc
# ============================================================
#
# Within each Fcap, each Besc is compared with the immediately
# preceding Besc value. Because Besc increments are unequal,
# raw deltas and deltas per 1000 t Besc are reported.
# ============================================================

table_marginal_Besc <- performance_ref %>%
  group_by(Fcap) %>%
  arrange(Besc,.by_group=TRUE) %>%
  mutate(
    Besc_previous=lag(Besc),
    delta_Besc=Besc-Besc_previous,
    delta_mean_catch=mean_catch-lag(mean_catch),
    delta_mean_P_closed=mean_P_closed-lag(mean_P_closed),
    delta_mean_years_closed=mean_years_closed-lag(mean_years_closed),
    delta_median_IAV=median_IAV-lag(median_IAV),
    delta_median_mean_SSB_open=median_mean_SSB_open-lag(median_mean_SSB_open),
    delta_max_P_Blim=max_P_Blim-lag(max_P_Blim),
    delta_P_ever_Blim=P_ever_Blim-lag(P_ever_Blim),
    delta_mean_catch_per_1000t_Besc=delta_mean_catch/delta_Besc*1000,
    delta_mean_P_closed_per_1000t_Besc=delta_mean_P_closed/delta_Besc*1000,
    delta_mean_years_closed_per_1000t_Besc=delta_mean_years_closed/delta_Besc*1000,
    delta_max_P_Blim_per_1000t_Besc=delta_max_P_Blim/delta_Besc*1000
  ) %>%
  filter(!is.na(Besc_previous)) %>%
  ungroup() %>%
  select(
    scenario_number,
    scenario_id,
    Fcap,
    Besc_previous,
    Besc,
    delta_Besc,
    mean_catch,
    mean_P_closed,
    mean_years_closed,
    median_IAV,
    median_mean_SSB_open,
    max_P_Blim,
    P_ever_Blim,
    delta_mean_catch,
    delta_mean_P_closed,
    delta_mean_years_closed,
    delta_median_IAV,
    delta_median_mean_SSB_open,
    delta_max_P_Blim,
    delta_P_ever_Blim,
    delta_mean_catch_per_1000t_Besc,
    delta_mean_P_closed_per_1000t_Besc,
    delta_mean_years_closed_per_1000t_Besc,
    delta_max_P_Blim_per_1000t_Besc
  ) %>%
  arrange(Fcap,Besc)

stopifnot(nrow(table_marginal_Besc)==length(expected_Fcap)*(length(expected_Besc)-1))
stopifnot(all(table_marginal_Besc$delta_Besc>0))

# ============================================================
# 7. TARGETED REFINED-GRID VIEWS
# ============================================================
#
# These subsets are descriptive and do not imply selection.
# They isolate:
#   - the refined upper Fcap region 1.50-1.80; and
#   - the central Besc region 6561-8000.
# ============================================================

diagnostic_Fcap_150_180 <- performance_ref %>%
  filter(Fcap>=1.50) %>%
  arrange(Besc,Fcap)

diagnostic_Besc_6561_8000 <- performance_ref %>%
  filter(Besc %in% c(6561,7250,8000)) %>%
  arrange(Fcap,Besc)

stopifnot(nrow(diagnostic_Fcap_150_180)==5*length(expected_Besc))
stopifnot(nrow(diagnostic_Besc_6561_8000)==3*length(expected_Fcap))

# ============================================================
# 8. SAVE DIAGNOSTIC OUTPUTS
# ============================================================

write.csv(
  performance_ref,
  file.path(output_dir,"Table24c_01_reference_metrics.csv"),
  row.names=FALSE
)

write.csv(
  table_marginal_Fcap,
  file.path(output_dir,"Table24c_02_marginal_Fcap.csv"),
  row.names=FALSE
)

write.csv(
  table_marginal_Besc,
  file.path(output_dir,"Table24c_03_marginal_Besc.csv"),
  row.names=FALSE
)

write.csv(
  diagnostic_Fcap_150_180,
  file.path(output_dir,"Table24c_04_Fcap_150_180.csv"),
  row.names=FALSE
)

write.csv(
  diagnostic_Besc_6561_8000,
  file.path(output_dir,"Table24c_05_Besc_6561_8000.csv"),
  row.names=FALSE
)

diagnostic_outputs <- list(
  input_file=input_file,
  expected_Fcap=expected_Fcap,
  expected_Besc=expected_Besc,
  risk_threshold=risk_threshold,
  precautionarity_summary=precautionarity_summary,
  reference_performance=performance_ref,
  marginal_Fcap=table_marginal_Fcap,
  marginal_Besc=table_marginal_Besc,
  diagnostic_Fcap_150_180=diagnostic_Fcap_150_180,
  diagnostic_Besc_6561_8000=diagnostic_Besc_6561_8000
)

saveRDS(
  diagnostic_outputs,
  file.path(output_dir,"refined_MP_tradeoffs_reference_OM.rds")
)

# ============================================================
# 9. CONSOLE DIAGNOSTICS
# ============================================================

cat("\n============================================================\n")
cat("PRECAUTIONARITY SUMMARY\n")
cat("============================================================\n")

print(precautionarity_summary)

cat("\n============================================================\n")
cat("MARGINAL CHANGES WITH INCREASING Fcap\n")
cat("============================================================\n")

print(
  table_marginal_Fcap %>%
    select(
      Besc,
      Fcap_previous,
      Fcap,
      delta_Fcap,
      delta_mean_catch,
      delta_mean_catch_per_0.10_Fcap,
      delta_mean_P_closed,
      delta_mean_years_closed,
      delta_max_P_Blim
    ),
  n=Inf
)

cat("\n============================================================\n")
cat("MARGINAL CHANGES WITH INCREASING Besc\n")
cat("============================================================\n")

print(
  table_marginal_Besc %>%
    select(
      Fcap,
      Besc_previous,
      Besc,
      delta_Besc,
      delta_mean_catch,
      delta_mean_catch_per_1000t_Besc,
      delta_mean_P_closed,
      delta_mean_years_closed,
      delta_max_P_Blim
    ),
  n=Inf
)

cat("\n============================================================\n")
cat("REFINED Fcap REGION: 1.50-1.80\n")
cat("============================================================\n")

diagnostic_Fcap_150_180 %>%
  select(
    scenario_number,
    Fcap,
    Besc,
    max_P_Blim,
    P_ever_Blim,
    mean_catch,
    mean_P_closed,
    mean_years_closed,
    median_IAV
  ) %>%
  as_tibble() %>%
  print(n=Inf)

cat("\n============================================================\n")
cat("CENTRAL Besc REGION: 6561-8000\n")
cat("============================================================\n")

diagnostic_Besc_6561_8000 %>%
  select(
    scenario_number,
    Fcap,
    Besc,
    max_P_Blim,
    P_ever_Blim,
    mean_catch,
    mean_P_closed,
    mean_years_closed,
    median_IAV
  ) %>%
  as_tibble() %>%
  print(n=Inf)
# ============================================================
# 10. FINAL VALIDATION CHECKS
# ============================================================

stopifnot(nrow(performance_ref)==28)
stopifnot(nrow(table_marginal_Fcap)==24)
stopifnot(nrow(table_marginal_Besc)==21)
stopifnot(nrow(diagnostic_Fcap_150_180)==20)
stopifnot(nrow(diagnostic_Besc_6561_8000)==21)
stopifnot(all(performance_ref$max_P_Blim>=0 & performance_ref$max_P_Blim<=1))
stopifnot(all(performance_ref$P_ever_Blim>=0 & performance_ref$P_ever_Blim<=1))
stopifnot(all(performance_ref$mean_P_closed>=0 & performance_ref$mean_P_closed<=1))
stopifnot(all(performance_ref$mean_years_closed>=0))
stopifnot(all(performance_ref$mean_catch>=0))

cat("\n============================================================\n")
cat("24c REFINED TRADE-OFF DIAGNOSTICS COMPLETED\n")
cat("============================================================\n")
cat("Candidate MPs:",nrow(performance_ref),"\n")
cat("Fcap marginal comparisons:",nrow(table_marginal_Fcap),"\n")
cat("Besc marginal comparisons:",nrow(table_marginal_Besc),"\n")
cat("Maximum annual P(SSB < Blim):",max(performance_ref$max_P_Blim),"\n")
cat("MPs with annual P(SSB < Blim) > 0.05:",sum(performance_ref$max_P_Blim>risk_threshold),"\n")
cat("All structural checks: PASS\n")
cat("Outputs saved in:",output_dir,"\n")
