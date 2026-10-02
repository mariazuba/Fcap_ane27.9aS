# ============================================================
# 24b_DIAGNOSE_MARGINAL_MP_TRADEOFFS.R
# Marginal trade-off diagnostics for candidate MPs
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
#   24b_DIAGNOSE_MARGINAL_MP_TRADEOFFS.R
#          |
#   Document candidate range for Stage 3 robustness
#
# Purpose:
#   Quantify marginal changes in already validated Reference-OM
#   performance metrics when:
#
#     A. Fcap increases within a fixed Besc; and
#     B. Besc increases within a fixed Fcap.
#
#   The diagnostic is intended to identify whether the candidate grid
#   contains regions of diminishing fishery returns or strong increases
#   in closure burden before defining the set of MPs carried forward
#   to robustness testing.
#
#   The script DOES NOT:
#     - rerun FLBEIA;
#     - recalculate MSE performance metrics;
#     - define a new performance metric;
#     - apply an automatic acceptance threshold;
#     - rank MPs;
#     - automatically discard Fcap=2 or Besc=10000;
#     - select a preferred MP.
#
# Input:
#   outputs/mse/candidate_MP_selection/
#     Table24c_reference_OM_consequences.csv
#
# Outputs:
#   outputs/mse/candidate_MP_selection/diagnostics/
#     Table24b_01_marginal_Fcap.csv
#     Table24b_02_marginal_Besc.csv
#     Table24b_03_Fcap_reference_metrics.csv
#     Table24b_04_Besc_reference_metrics.csv
#     marginal_MP_tradeoffs_reference_OM.rds
#
# Interpretation:
#   Positive delta values always mean that the metric increased from
#   the previous candidate level to the current candidate level.
#
#   No delta is interpreted automatically as beneficial or adverse.
#   Interpretation depends on the metric:
#     - higher mean_catch = greater average yield;
#     - higher mean_P_closed / mean_years_closed = greater closure burden;
#     - higher max_P_Blim = greater biological risk;
#     - higher median_mean_SSB_open = higher SSB in open years;
#     - IAV must be interpreted together with closure metrics because
#       consecutive zero-catch years can reduce IAV.
#
# Important:
#   This is a diagnostic of marginal consequences, not a formal
#   selection rule. Any exclusion before Stage 3 must be documented
#   separately from these calculations.
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

input_file <- here("outputs","mse","candidate_MP_selection","Table24c_reference_OM_consequences.csv")
output_dir <- here("outputs","mse","candidate_MP_selection","diagnostics")

expected_Fcap <- c(1.00,1.25,1.50,1.75,2.00)
expected_Besc <- c(5000,6561,8000,10000)
expected_scenarios <- length(expected_Fcap)*length(expected_Besc)

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

stopifnot(file.exists(input_file))

# ============================================================
# 2. LOAD REFERENCE-OM CONSEQUENCES FROM SCRIPT 24
# ============================================================

performance_ref <- read.csv(input_file,stringsAsFactors=FALSE)

required_columns <- c(
  "scenario_number",
  "scenario_id",
  "Fcap",
  "Besc",
  "max_P_Blim",
  "n_years_P_Blim_gt_005",
  "P_ever_Blim",
  "min_SSB",
  "median_mean_SSB_open",
  "mean_catch",
  "mean_catch_open",
  "mean_P_closed",
  "mean_years_closed",
  "median_years_closed",
  "mean_IAV",
  "median_IAV"
)

missing_columns <- setdiff(required_columns,names(performance_ref))

stopifnot(length(missing_columns)==0)
stopifnot(nrow(performance_ref)==expected_scenarios)
stopifnot(dplyr::n_distinct(performance_ref$scenario_number)==expected_scenarios)
stopifnot(setequal(sort(unique(performance_ref$Fcap)),expected_Fcap))
stopifnot(setequal(sort(unique(performance_ref$Besc)),expected_Besc))
stopifnot(!anyDuplicated(performance_ref[c("Fcap","Besc")]))

performance_ref <- performance_ref %>%
  arrange(Besc,Fcap)

cat("\n============================================================\n")
cat("MARGINAL MP TRADE-OFF DIAGNOSTICS\n")
cat("============================================================\n")
cat("Candidate MPs:",nrow(performance_ref),"\n")
cat("Fcap values:",paste(sort(unique(performance_ref$Fcap)),collapse=", "),"\n")
cat("Besc values:",paste(sort(unique(performance_ref$Besc)),collapse=", "),"\n")

# ============================================================
# 3. MARGINAL CHANGES WITH INCREASING Fcap
# ============================================================
#
# Comparison:
#   within each Besc, each Fcap is compared with the immediately
#   preceding Fcap value.
#
# Example:
#   Fcap 1.75 versus Fcap 1.50 at the same Besc.
#
# The first Fcap within each Besc has no preceding candidate and is
# therefore omitted from the marginal table.
# ============================================================

table_marginal_Fcap <- performance_ref %>%
  group_by(Besc) %>%
  arrange(Fcap,.by_group=TRUE) %>%
  mutate(
    Fcap_previous=lag(Fcap),
    delta_Fcap=Fcap-lag(Fcap),
    delta_mean_catch=mean_catch-lag(mean_catch),
    delta_mean_catch_open=mean_catch_open-lag(mean_catch_open),
    delta_mean_P_closed=mean_P_closed-lag(mean_P_closed),
    delta_mean_years_closed=mean_years_closed-lag(mean_years_closed),
    delta_max_P_Blim=max_P_Blim-lag(max_P_Blim),
    delta_P_ever_Blim=P_ever_Blim-lag(P_ever_Blim),
    delta_median_mean_SSB_open=median_mean_SSB_open-lag(median_mean_SSB_open),
    delta_median_IAV=median_IAV-lag(median_IAV)
  ) %>%
  filter(!is.na(Fcap_previous)) %>%
  select(
    scenario_number,
    Besc,
    Fcap_previous,
    Fcap,
    delta_Fcap,
    mean_catch,
    delta_mean_catch,
    mean_catch_open,
    delta_mean_catch_open,
    mean_P_closed,
    delta_mean_P_closed,
    mean_years_closed,
    delta_mean_years_closed,
    max_P_Blim,
    delta_max_P_Blim,
    P_ever_Blim,
    delta_P_ever_Blim,
    median_mean_SSB_open,
    delta_median_mean_SSB_open,
    median_IAV,
    delta_median_IAV
  ) %>%
  ungroup() %>%
  arrange(Besc,Fcap)

stopifnot(nrow(table_marginal_Fcap)==length(expected_Besc)*(length(expected_Fcap)-1))
stopifnot(all(abs(table_marginal_Fcap$delta_Fcap-0.25)<1e-10))

# ============================================================
# 4. MARGINAL YIELD GAIN PER ADDITIONAL CLOSURE BURDEN
# ============================================================
#
# These are descriptive ratios based on existing metrics.
#
# They are NOT performance criteria and are NOT used to rank MPs.
#
# catch_gain_per_1pct_closure:
#   additional tonnes of mean catch associated with one percentage
#   point increase in mean annual closure probability.
#
# catch_gain_per_extra_closed_year:
#   additional tonnes of mean catch associated with one additional
#   mean closed year over the Full horizon.
#
# Ratios are NA when the corresponding denominator is zero.
# ============================================================

table_marginal_Fcap <- table_marginal_Fcap %>%
  mutate(
    delta_P_closed_percentage_points=100*delta_mean_P_closed,
    catch_gain_per_1pct_closure=ifelse(abs(delta_P_closed_percentage_points)>1e-12,delta_mean_catch/delta_P_closed_percentage_points,NA_real_),
    catch_gain_per_extra_closed_year=ifelse(abs(delta_mean_years_closed)>1e-12,delta_mean_catch/delta_mean_years_closed,NA_real_)
  )

# ============================================================
# 5. MARGINAL CHANGES WITH INCREASING Besc
# ============================================================
#
# Comparison:
#   within each Fcap, each Besc is compared with the immediately
#   preceding Besc value.
#
# Besc increments are not equal:
#   5000 -> 6561 -> 8000 -> 10000.
#
# Therefore both absolute deltas and changes per 1000 t increase
# in Besc are retained.
# ============================================================

table_marginal_Besc <- performance_ref %>%
  group_by(Fcap) %>%
  arrange(Besc,.by_group=TRUE) %>%
  mutate(
    Besc_previous=lag(Besc),
    delta_Besc=Besc-lag(Besc),
    delta_mean_catch=mean_catch-lag(mean_catch),
    delta_mean_catch_open=mean_catch_open-lag(mean_catch_open),
    delta_mean_P_closed=mean_P_closed-lag(mean_P_closed),
    delta_mean_years_closed=mean_years_closed-lag(mean_years_closed),
    delta_max_P_Blim=max_P_Blim-lag(max_P_Blim),
    delta_P_ever_Blim=P_ever_Blim-lag(P_ever_Blim),
    delta_median_mean_SSB_open=median_mean_SSB_open-lag(median_mean_SSB_open),
    delta_median_IAV=median_IAV-lag(median_IAV)
  ) %>%
  filter(!is.na(Besc_previous)) %>%
  mutate(
    delta_mean_catch_per_1000t_Besc=1000*delta_mean_catch/delta_Besc,
    delta_mean_catch_open_per_1000t_Besc=1000*delta_mean_catch_open/delta_Besc,
    delta_P_closed_per_1000t_Besc=1000*delta_mean_P_closed/delta_Besc,
    delta_closed_years_per_1000t_Besc=1000*delta_mean_years_closed/delta_Besc,
    delta_SSB_open_per_1000t_Besc=1000*delta_median_mean_SSB_open/delta_Besc
  ) %>%
  select(
    scenario_number,
    Fcap,
    Besc_previous,
    Besc,
    delta_Besc,
    mean_catch,
    delta_mean_catch,
    delta_mean_catch_per_1000t_Besc,
    mean_catch_open,
    delta_mean_catch_open,
    delta_mean_catch_open_per_1000t_Besc,
    mean_P_closed,
    delta_mean_P_closed,
    delta_P_closed_per_1000t_Besc,
    mean_years_closed,
    delta_mean_years_closed,
    delta_closed_years_per_1000t_Besc,
    max_P_Blim,
    delta_max_P_Blim,
    P_ever_Blim,
    delta_P_ever_Blim,
    median_mean_SSB_open,
    delta_median_mean_SSB_open,
    delta_SSB_open_per_1000t_Besc,
    median_IAV,
    delta_median_IAV
  ) %>%
  ungroup() %>%
  arrange(Fcap,Besc)

stopifnot(nrow(table_marginal_Besc)==length(expected_Fcap)*(length(expected_Besc)-1))
stopifnot(all(table_marginal_Besc$delta_Besc>0))

# ============================================================
# 6. REFERENCE METRICS BY Fcap
# ============================================================
#
# This table keeps the original Reference-OM metrics together with
# the immediately preceding Fcap value. It is intended for direct
# inspection of the upper end of the Fcap grid, including Fcap=2.
# ============================================================

table_Fcap_reference <- performance_ref %>%
  group_by(Besc) %>%
  arrange(Fcap,.by_group=TRUE) %>%
  mutate(
    Fcap_previous=lag(Fcap),
    delta_mean_catch=mean_catch-lag(mean_catch),
    delta_mean_P_closed=mean_P_closed-lag(mean_P_closed),
    delta_mean_years_closed=mean_years_closed-lag(mean_years_closed),
    delta_max_P_Blim=max_P_Blim-lag(max_P_Blim)
  ) %>%
  ungroup() %>%
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
    median_years_closed,
    median_mean_SSB_open,
    median_IAV,
    Fcap_previous,
    delta_mean_catch,
    delta_mean_P_closed,
    delta_mean_years_closed,
    delta_max_P_Blim
  ) %>%
  arrange(Besc,Fcap)

# ============================================================
# 7. REFERENCE METRICS BY Besc
# ============================================================
#
# This table keeps the original Reference-OM metrics together with
# the immediately preceding Besc value. It is intended for direct
# inspection of the upper end of the Besc grid, including Besc=10000.
# ============================================================

table_Besc_reference <- performance_ref %>%
  group_by(Fcap) %>%
  arrange(Besc,.by_group=TRUE) %>%
  mutate(
    Besc_previous=lag(Besc),
    delta_mean_catch=mean_catch-lag(mean_catch),
    delta_mean_P_closed=mean_P_closed-lag(mean_P_closed),
    delta_mean_years_closed=mean_years_closed-lag(mean_years_closed),
    delta_median_mean_SSB_open=median_mean_SSB_open-lag(median_mean_SSB_open),
    delta_max_P_Blim=max_P_Blim-lag(max_P_Blim)
  ) %>%
  ungroup() %>%
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
    median_years_closed,
    median_mean_SSB_open,
    median_IAV,
    Besc_previous,
    delta_mean_catch,
    delta_mean_P_closed,
    delta_mean_years_closed,
    delta_median_mean_SSB_open,
    delta_max_P_Blim
  ) %>%
  arrange(Fcap,Besc)

# ============================================================
# 8. TARGETED DIAGNOSTIC VIEWS
# ============================================================
#
# These subsets do not imply rejection.
# They isolate the two upper boundaries currently under examination:
#   Fcap = 2
#   Besc = 10000
# ============================================================

diagnostic_Fcap2 <- table_Fcap_reference %>%
  filter(Fcap==max(expected_Fcap)) %>%
  arrange(Besc)

diagnostic_Besc10000 <- table_Besc_reference %>%
  filter(Besc==max(expected_Besc)) %>%
  arrange(Fcap)

stopifnot(nrow(diagnostic_Fcap2)==length(expected_Besc))
stopifnot(nrow(diagnostic_Besc10000)==length(expected_Fcap))

# ============================================================
# 9. SAVE DIAGNOSTIC OUTPUTS
# ============================================================

write.csv(
  table_marginal_Fcap,
  file.path(output_dir,"Table24b_01_marginal_Fcap.csv"),
  row.names=FALSE
)

write.csv(
  table_marginal_Besc,
  file.path(output_dir,"Table24b_02_marginal_Besc.csv"),
  row.names=FALSE
)

write.csv(
  table_Fcap_reference,
  file.path(output_dir,"Table24b_03_Fcap_reference_metrics.csv"),
  row.names=FALSE
)

write.csv(
  table_Besc_reference,
  file.path(output_dir,"Table24b_04_Besc_reference_metrics.csv"),
  row.names=FALSE
)

diagnostic_outputs <- list(
  input_file=input_file,
  expected_Fcap=expected_Fcap,
  expected_Besc=expected_Besc,
  reference_performance=performance_ref,
  marginal_Fcap=table_marginal_Fcap,
  marginal_Besc=table_marginal_Besc,
  Fcap_reference=table_Fcap_reference,
  Besc_reference=table_Besc_reference,
  diagnostic_Fcap2=diagnostic_Fcap2,
  diagnostic_Besc10000=diagnostic_Besc10000
)

saveRDS(
  diagnostic_outputs,
  file.path(output_dir,"marginal_MP_tradeoffs_reference_OM.rds")
)

# ============================================================
# 10. CONSOLE DIAGNOSTICS
# ============================================================

cat("\n============================================================\n")
cat("MARGINAL CHANGES WITH INCREASING Fcap\n")
cat("============================================================\n")

print(
  table_marginal_Fcap %>%
    select(
      Besc,
      Fcap_previous,
      Fcap,
      delta_mean_catch,
      delta_mean_P_closed,
      delta_mean_years_closed,
      delta_max_P_Blim,
      delta_median_IAV
    ),
  n=Inf
)

cat("\n============================================================\n")
cat("UPPER Fcap BOUNDARY: Fcap = 2\n")
cat("============================================================\n")

print(
  diagnostic_Fcap2 %>%
    select(
      Besc,
      Fcap_previous,
      Fcap,
      mean_catch,
      delta_mean_catch,
      mean_P_closed,
      delta_mean_P_closed,
      mean_years_closed,
      delta_mean_years_closed,
      max_P_Blim,
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
      delta_mean_catch,
      delta_mean_P_closed,
      delta_mean_years_closed,
      delta_median_mean_SSB_open,
      delta_max_P_Blim
    ),
  n=Inf
)

cat("\n============================================================\n")
cat("UPPER Besc BOUNDARY: Besc = 10000\n")
cat("============================================================\n")

print(
  diagnostic_Besc10000 %>%
    select(
      Fcap,
      Besc_previous,
      Besc,
      mean_catch,
      delta_mean_catch,
      mean_P_closed,
      delta_mean_P_closed,
      mean_years_closed,
      delta_mean_years_closed,
      median_mean_SSB_open,
      delta_median_mean_SSB_open,
      max_P_Blim,
      delta_max_P_Blim
    ),
  n=Inf
)

# ============================================================
# 11. FINAL VALIDATION CHECKS
# ============================================================

stopifnot(file.exists(file.path(output_dir,"Table24b_01_marginal_Fcap.csv")))
stopifnot(file.exists(file.path(output_dir,"Table24b_02_marginal_Besc.csv")))
stopifnot(file.exists(file.path(output_dir,"Table24b_03_Fcap_reference_metrics.csv")))
stopifnot(file.exists(file.path(output_dir,"Table24b_04_Besc_reference_metrics.csv")))
stopifnot(file.exists(file.path(output_dir,"marginal_MP_tradeoffs_reference_OM.rds")))
stopifnot(nrow(table_marginal_Fcap)==16)
stopifnot(nrow(table_marginal_Besc)==15)
stopifnot(nrow(diagnostic_Fcap2)==4)
stopifnot(nrow(diagnostic_Besc10000)==5)

cat("\n============================================================\n")
cat("MARGINAL MP TRADE-OFF DIAGNOSTICS COMPLETED\n")
cat("============================================================\n")
cat("Candidate MPs:",expected_scenarios,"\n")
cat("Fcap marginal comparisons:",nrow(table_marginal_Fcap),"\n")
cat("Besc marginal comparisons:",nrow(table_marginal_Besc),"\n")
cat("Output directory:",output_dir,"\n")
cat("\nNo candidate MP has been automatically excluded.\n")
cat("The outputs quantify the marginal consequences needed to document\n")
cat("the candidate range carried forward to Stage 3 robustness.\n")

