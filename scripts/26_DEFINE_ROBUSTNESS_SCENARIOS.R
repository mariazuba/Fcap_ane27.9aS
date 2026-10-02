# ============================================================
# 26_DEFINE_ROBUSTNESS_SCENARIOS.R
# Formal primary robustness design for anchovy ane.27.9aS
# Reads the validated script-25 selection; creates metadata only.
# Does not modify model inputs, the HCR or the runner; runs no MSE.
# Scenario definitions follow prompt_escenarios_robustez_30sep.txt.
# ============================================================

rm(list=ls())
library(dplyr)
library(here)

# 1. PATHS ----------------------------------------------------
selection_dir <- here("outputs", "mse", "candidate_MP_selection", "robustness")
selection_file <- file.path(selection_dir, "selected_MPs_for_robustness.rds")
output_dir <- here("outputs", "mse", "robustness", "design")
stopifnot(file.exists(selection_file))

# 2. VALIDATED MP SELECTION -----------------------------------
selection <- readRDS(selection_file)
selected_MPs <- selection$selected_MPs
required_columns <- c("scenario_number", "scenario_id", "Fcap", "Besc", "carry_forward_to_robustness")
stopifnot(all(required_columns %in% names(selected_MPs)))
stopifnot(nrow(selected_MPs)==12)
stopifnot(!anyNA(selected_MPs[required_columns]))
stopifnot(all(selected_MPs$carry_forward_to_robustness))
stopifnot(n_distinct(selected_MPs$scenario_id)==12)
stopifnot(n_distinct(selected_MPs$scenario_number)==12)
stopifnot(setequal(unique(selected_MPs$Fcap), c(1.50, 1.70, 1.80)))
stopifnot(setequal(unique(selected_MPs$Besc), c(5000, 6561, 7250, 8000)))
stopifnot(!anyDuplicated(selected_MPs[c("Fcap", "Besc")]))
stopifnot(all(selected_MPs$max_P_Blim<=0.05))

# 3. PRIMARY SCENARIO DEFINITIONS ------------------------------
# sigmaR is specified by its source, not by rounded constants.
# REF: use the existing reference SR process-error input unchanged.
# R-OM: use source-OM sigmaR and the same standard-normal deviates.
# REF seasonal shares are OM-specific; NA means retain these inputs.
# C1/C2 are OM catch-allocation vectors, never MP F proportions.
# I0 is the perfect-information endpoint: both AE and forecast
# uncertainty differ from REF. I0 vs I1 isolates AE; I1 vs REF
# isolates forecast uncertainty. Thus the whole design is mainly OFAT.

scenarios <- tibble(
  robustness_id=c("REF", "R-OM", "S-C1", "S-C2", "I0", "I1"),
  scenario_type=c("reference", "recruitment_sensitivity", "seasonal_sensitivity", "seasonal_sensitivity", "information_endpoint", "information_sensitivity"),
  recruitment_PE_mode=c("reference_common", "source_OM_sigmaR", "reference_common", "reference_common", "reference_common", "reference_common"),
  OM_seasonal_share_mode=c("conditioned_OM_recent10", "conditioned_OM_recent10", "fixed_C1", "fixed_C2", "conditioned_OM_recent10", "conditioned_OM_recent10"),
  OM_share_Q1=c(NA_real_, NA_real_, 0.271, 0.235, NA_real_, NA_real_),
  OM_share_Q2=c(NA_real_, NA_real_, 0.534, 0.309, NA_real_, NA_real_),
  OM_share_Q3=c(NA_real_, NA_real_, 0.164, 0.299, NA_real_, NA_real_),
  OM_share_Q4=c(NA_real_, NA_real_, 0.031, 0.157, NA_real_, NA_real_),
  assessment_error_SSB_F=c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  recruitment_assessment_error=FALSE,
  forecast_uncertainty=c("S2", "S2", "S2", "S2", "S0", "S0"),
  comparison_baseline=c(NA_character_, "REF", "REF", "REF", "I1", "REF"),
  scientific_purpose=c(
    "Reference OM/MP configuration.",
    "Effect of source-OM recruitment process-error SD using paired standard-normal deviates.",
    "Effect of historical cluster-1 OM seasonal catch allocation.",
    "Effect of historical cluster-2 OM seasonal catch allocation.",
    "Perfect-state, deterministic-forecast endpoint; compare with I1 for assessment error and REF for combined information effect.",
    "Historical assessment error with deterministic forecast; compare with REF for forecast uncertainty."
  ),
  implementation_status=c(
    "VALIDATED_REFERENCE", "IMPLEMENTED_PENDING_FULL_MSE_VALIDATION",
    "PENDING_IMPLEMENTATION_AND_VALIDATION", "PENDING_IMPLEMENTATION_AND_VALIDATION",
    "IMPLEMENTATION_AND_EVALUATION_STATUS_TO_CHECK", "IMPLEMENTATION_AND_EVALUATION_STATUS_TO_CHECK"
  )
)

scenarios <- scenarios %>%
  mutate(projection_start=2025L, projection_end=2054L) %>%
  mutate(Blim=4721, Bpa=6561, annual_risk_threshold=0.05) %>%
  mutate(MP_propf_method="get_Fprop", MP_nFyears=5L) %>%
  mutate(n_OMs=100L, replicates_per_OM=10L, n_trajectories=1000L, n_blocks=10L) %>%
  mutate(recruitment_pairing="same_standard_normal_deviates_as_REF") %>%
  mutate(assessment_pairing="same_historical_AE_draws_where_enabled") %>%
  mutate(forecast_pairing="same_empirical_S2_CV_by_global_trajectory_where_enabled") %>%
  mutate(reuse_existing_REF=robustness_id=="REF")

# 4. VALIDATE DEFINITIONS -------------------------------------
stopifnot(identical(scenarios$robustness_id, c("REF", "R-OM", "S-C1", "S-C2", "I0", "I1")))
share_columns <- paste0("OM_share_Q", 1:4)
fixed_rows <- scenarios$robustness_id %in% c("S-C1", "S-C2")
fixed_shares <- as.matrix(scenarios[fixed_rows, share_columns])
stopifnot(all(is.finite(fixed_shares)))
stopifnot(all(fixed_shares>=0 & fixed_shares<=1))
stopifnot(all(abs(rowSums(fixed_shares)-1)<1e-12))
stopifnot(all(is.na(as.matrix(scenarios[!fixed_rows, share_columns]))))
stopifnot(all(scenarios$MP_nFyears==5L))
stopifnot(!any(scenarios$recruitment_assessment_error))

# 5. CROSS SELECTED MPs WITH THE SIX CONFIGURATIONS ------------
# Keep original REF identifiers unchanged. robustness_run_id is a
# separate design key and is not a runner filename convention.
MP_keys <- selected_MPs %>%
  select(scenario_number, scenario_id, Fcap, Besc) %>%
  rename(reference_scenario_number=scenario_number, reference_scenario_id=scenario_id)
run_grid <- merge(as.data.frame(MP_keys), as.data.frame(scenarios), by=NULL) %>%
  as_tibble() %>%
  mutate(robustness_run_id=paste(robustness_id, reference_scenario_id, sep="__")) %>%
  arrange(match(robustness_id, scenarios$robustness_id), Besc, Fcap)
stopifnot(nrow(run_grid)==72)
stopifnot(n_distinct(run_grid$robustness_run_id)==72)
stopifnot(all(table(run_grid$robustness_id)==12))
new_runs <- run_grid %>% filter(!reuse_existing_REF)
stopifnot(nrow(new_runs)==60)

# 6. SAVE FORMAL DESIGN ---------------------------------------
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(scenarios, file.path(output_dir, "Table26_01_robustness_scenarios.csv"), row.names=FALSE, na="")
write.csv(run_grid, file.path(output_dir, "Table26_02_MP_robustness_design.csv"), row.names=FALSE, na="")
write.csv(new_runs, file.path(output_dir, "Table26_03_new_MP_robustness_runs.csv"), row.names=FALSE, na="")
robustness_design <- list(
  definition_status="FORMAL_DESIGN_IMPLEMENTED_PENDING_RUNNER_VALIDATION",
  source_selection_file=selection_file,
  methodological_source="prompt_escenarios_robustez_30sep.txt",
  selected_MPs=selected_MPs,
  scenarios=scenarios,
  run_grid=run_grid,
  new_runs=new_runs,
  information_comparisons=c("I0 vs I1: assessment error", "I1 vs REF: forecast uncertainty", "I0 vs REF: combined information effect"),
  implementation_notes=c(
    "Scenario CSVs are definitions; they do not implement model changes.",
    "Retain exact reference sigmaR inputs; approximately 0.216 is descriptive.",
    "R-OM source sigmaR range approximately 0.168-0.278 is descriptive, not an imposed bound.",
    "S-C1/S-C2 overwrite only fleets.ctrl$seasonal.share during 2025-2054.",
    "Keep historical inputs through 2024 unchanged.",
    "Retain dynamic MP get_Fprop with nFyears=5; never impose C1/C2 on propf.",
    "S0 means CV=0; S2 uses the existing empirical trajectory-paired CV.",
    "Reuse existing REF results after verifying identifiers and pairing.",
    "Verify runner mappings and validate each scenario before CESGA submission."
  )
)
saveRDS(robustness_design, file.path(output_dir, "robustness_design.rds"))

# 7. SUMMARY --------------------------------------------------
cat("\nROBUSTNESS DESIGN\n")
print(scenarios %>% select(robustness_id, recruitment_PE_mode, OM_seasonal_share_mode, assessment_error_SSB_F, forecast_uncertainty), n=Inf, width=Inf)
cat("\nSelected MPs:", nrow(MP_keys), "\n")
cat("Configurations:", nrow(scenarios), "\n")
cat("MP x configuration combinations:", nrow(run_grid), "\n")
cat("Existing REF combinations to reuse: 12\n")
cat("New MP x alternative combinations:", nrow(new_runs), "\n")
cat("New block runs if using 10 blocks per combination:", sum(new_runs$n_blocks), "\n")
cat("Scenario definitions and complete cross-design: PASS\n")
cat("Runner implementation and scenario validation: PENDING\n")
cat("Outputs saved in:\n", output_dir, "\n", sep="")
cat("\nSCRIPT 26 COMPLETED SUCCESSFULLY\n")
