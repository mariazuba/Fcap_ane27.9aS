# ============================================================
# 25_SELECT_REFINED_MPS_FOR_ROBUSTNESS.R
# Scientific selection of refined candidate MPs for robustness
# European anchovy ane.27.9aS (Gulf of Cadiz)
# Shortcut MSE in FLBEIA + FLasher conditioned from SS3
#
# Workflow: 22b -> 24c -> 25 -> robustness scenario definition
#
# Purpose: document the reduction of the definitive 28-MP grid
# to 12 representative MPs. This is not a ranking or final MP
# selection. No simulations or performance metrics are recalculated.
#
# Selected Fcap: 1.50, 1.70, 1.80
# Selected Besc: 5000, 6561, 7250, 8000 t
# All 28 MPs remain in the master selection table.
# ============================================================

rm(list=ls())

library(dplyr)
library(here)

# 1. PATHS ----------------------------------------------------
performance_file <- here("data", "mse", "performance_refined", "refined_candidate_MP_performance_full.csv")
tradeoff_file <- here("outputs", "mse", "candidate_MP_selection", "refined_diagnostics", "refined_MP_tradeoffs_reference_OM.rds")
output_dir <- here("outputs", "mse", "candidate_MP_selection", "robustness")

stopifnot(file.exists(performance_file))
stopifnot(file.exists(tradeoff_file))

# 2. EXPECTED GRID AND AGREED SELECTION ------------------------
expected_Fcap <- c(1.00, 1.25, 1.50, 1.60, 1.70, 1.75, 1.80)
expected_Besc <- c(5000, 6561, 7250, 8000)
expected_n_MPs <- length(expected_Fcap)*length(expected_Besc)
risk_threshold <- 0.05

# Fcap 1.50: lower representative level of the refined region.
# Fcap 1.70: intermediate representative level of the refined region.
# Fcap 1.80: upper representative level of the refined region.
# Fcap 1.60 and 1.75 are omitted to reduce closely spaced levels;
# this does not imply rejection of these MPs.
# Fcap 1.00 and 1.25 remain in the REF evaluation but are outside
# the refined Fcap region selected for robustness testing.
#
# All four Besc levels are retained to preserve escapement contrasts:
# 5000: lower; 6561: Bpa-linked; 7250: intermediate; 8000: higher.
# Biological differences may emerge under alternative OMs.
selected_Fcap <- c(1.50, 1.70, 1.80)
selected_Besc <- c(5000, 6561, 7250, 8000)
expected_selected_MPs <- length(selected_Fcap)*length(selected_Besc)
stopifnot(expected_selected_MPs==12)

# 3. LOAD VALIDATED REFERENCE-OM PERFORMANCE -------------------
performance_ref <- read.csv(performance_file, stringsAsFactors=FALSE)
required_columns <- c(
  "scenario_number", "scenario_id", "Fcap", "Besc",
  "max_P_Blim", "n_years_P_Blim_gt_5", "P_ever_Blim",
  "max_P_Bpa", "P_ever_Bpa", "min_SSB", "mean_catch",
  "mean_catch_open", "mean_P_closed", "mean_years_closed",
  "median_years_closed", "median_mean_SSB_open", "mean_IAV", "median_IAV"
)
missing_columns <- setdiff(required_columns, names(performance_ref))
if (length(missing_columns)>0) {
  stop("Missing performance columns: ", paste(missing_columns, collapse=", "))
}
stopifnot(nrow(performance_ref)==expected_n_MPs)
stopifnot(!anyNA(performance_ref[c("scenario_number", "scenario_id", "Fcap", "Besc")]))
stopifnot(all(nzchar(trimws(performance_ref$scenario_id))))
stopifnot(n_distinct(performance_ref$scenario_number)==expected_n_MPs)
stopifnot(n_distinct(performance_ref$scenario_id)==expected_n_MPs)
stopifnot(setequal(unique(performance_ref$Fcap), expected_Fcap))
stopifnot(setequal(unique(performance_ref$Besc), expected_Besc))
stopifnot(!anyDuplicated(performance_ref[c("Fcap", "Besc")]))

# 4. VERIFY SUPPORTING TRADE-OFF DIAGNOSTIC ---------------------
tradeoff_outputs <- readRDS(tradeoff_file)
required_tradeoff_objects <- c("precautionarity_summary", "reference_performance", "marginal_Fcap", "marginal_Besc")
stopifnot(all(required_tradeoff_objects %in% names(tradeoff_outputs)))
tradeoff_ref <- tradeoff_outputs$reference_performance
stopifnot(nrow(tradeoff_ref)==expected_n_MPs)
stopifnot(setequal(performance_ref$scenario_id, tradeoff_ref$scenario_id))

# Verify that 24c used the same reference table, allowing only
# differences in row/column order and CSV numeric precision.
stopifnot(all(required_columns %in% names(tradeoff_ref)))
tradeoff_match <- match(performance_ref$scenario_id, tradeoff_ref$scenario_id)
stopifnot(!anyNA(tradeoff_match))
for (column in required_columns) {
  comparison <- all.equal(
    performance_ref[[column]],
    tradeoff_ref[[column]][tradeoff_match],
    tolerance=1e-8,
    check.attributes=FALSE
  )
  if (!isTRUE(comparison)) {
    stop("22b/24c reference mismatch in column: ", column)
  }
}

# 5. VERIFY PRECAUTIONARITY PREMISE ----------------------------
# This verifies admissibility under REF; it does not rank MPs.
stopifnot(all(is.finite(performance_ref$max_P_Blim)))
stopifnot(all(performance_ref$max_P_Blim>=0))
stopifnot(all(performance_ref$max_P_Blim<=risk_threshold))
stopifnot(all(performance_ref$n_years_P_Blim_gt_5==0))

# 6. DOCUMENT SELECTION ---------------------------------------
selection_table <- performance_ref %>%
  mutate(carry_forward_to_robustness=Fcap %in% selected_Fcap & Besc %in% selected_Besc) %>%
  mutate(Fcap_role=case_when(
    Fcap==1.50 ~ "lower_refined_Fcap",
    Fcap==1.70 ~ "intermediate_refined_Fcap",
    Fcap==1.80 ~ "upper_refined_Fcap",
    TRUE ~ "not_carried_forward"
  )) %>%
  mutate(Besc_role=case_when(
    Besc==5000 ~ "lower_escapement_contrast",
    Besc==6561 ~ "Bpa_linked_escapement",
    Besc==7250 ~ "intermediate_escapement",
    Besc==8000 ~ "higher_escapement_contrast",
    TRUE ~ NA_character_
  )) %>%
  mutate(selection_status=if_else(carry_forward_to_robustness, "SELECTED_FOR_ROBUSTNESS", "NOT_SELECTED_FOR_ROBUSTNESS")) %>%
  mutate(selection_rationale=case_when(
    carry_forward_to_robustness ~ paste0("Representative Fcap level (", Fcap_role, ") crossed with Besc contrast (", Besc_role, ")."),
    Fcap %in% c(1.60, 1.75) ~ paste0("Retained in the Reference-OM candidate set but not carried forward to robustness because adjacent representative Fcap levels provide the required contrast without implying rejection of this MP."),
    Fcap %in% c(1.00, 1.25) ~ "Retained in the Reference-OM candidate set but outside the refined Fcap region selected for robustness testing.",
    TRUE ~ "Not carried forward to robustness."
  )) %>%
  arrange(Besc, Fcap)

selected_MPs <- selection_table %>%
  filter(carry_forward_to_robustness) %>%
  arrange(Besc, Fcap)

# 7. VALIDATE THE COMPLETE 3 x 4 SELECTION ---------------------
stopifnot(nrow(selection_table)==28)
stopifnot(nrow(selected_MPs)==12)
stopifnot(n_distinct(selected_MPs$scenario_number)==12)
stopifnot(n_distinct(selected_MPs$scenario_id)==12)
stopifnot(setequal(unique(selected_MPs$Fcap), selected_Fcap))
stopifnot(setequal(unique(selected_MPs$Besc), selected_Besc))
selected_design <- selected_MPs %>% count(Fcap, Besc, name="n")
stopifnot(nrow(selected_design)==12)
stopifnot(all(selected_design$n==1))
stopifnot(all(table(selected_MPs$Fcap)==length(selected_Besc)))
stopifnot(all(table(selected_MPs$Besc)==length(selected_Fcap)))
stopifnot(all(selected_MPs$max_P_Blim<=risk_threshold))
stopifnot(all(selected_MPs$n_years_P_Blim_gt_5==0))

# 8. OUTPUT TABLES --------------------------------------------
table_all <- selection_table %>%
  select(scenario_number, scenario_id, Fcap, Besc, max_P_Blim,
         n_years_P_Blim_gt_5, P_ever_Blim, mean_catch, mean_catch_open,
         mean_P_closed, mean_years_closed, median_IAV, Fcap_role,
         Besc_role, carry_forward_to_robustness, selection_status,
         selection_rationale)

table_selected <- selected_MPs %>%
  select(scenario_number, scenario_id, Fcap, Besc, max_P_Blim,
         P_ever_Blim, mean_catch, mean_catch_open, mean_P_closed,
         mean_years_closed, median_IAV, Fcap_role, Besc_role,
         selection_rationale)

# 9. SAVE ONLY AFTER ALL VALIDATIONS PASS ----------------------
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(table_all, file.path(output_dir, "Table25_01_all_refined_MPs_selection.csv"), row.names=FALSE)
write.csv(table_selected, file.path(output_dir, "Table25_02_MPs_selected_for_robustness.csv"), row.names=FALSE)

selection_outputs <- list(
  source_performance_file=performance_file,
  source_tradeoff_file=tradeoff_file,
  risk_threshold=risk_threshold,
  selected_Fcap=selected_Fcap,
  selected_Besc=selected_Besc,
  n_refined_MPs=nrow(selection_table),
  n_selected_MPs=nrow(selected_MPs),
  all_MPs=selection_table,
  selected_MPs=selected_MPs
)
saveRDS(selection_outputs, file.path(output_dir, "selected_MPs_for_robustness.rds"))

# 10. CONSOLE SUMMARY -----------------------------------------
cat("\n============================================================\n")
cat("MP SELECTION FOR ROBUSTNESS\n")
cat("============================================================\n")
cat("Refined candidate MPs:", nrow(selection_table), "\n")
cat("MPs selected for robustness:", nrow(selected_MPs), "\n")
cat("Selected Fcap:", paste(selected_Fcap, collapse=", "), "\n")
cat("Selected Besc:", paste(selected_Besc, collapse=", "), "\n")
cat("Maximum annual P(SSB < Blim) among selected MPs:", max(selected_MPs$max_P_Blim), "\n")
cat("\nSELECTED MPs\n")
print(as_tibble(table_selected %>% select(scenario_number, scenario_id, Fcap, Besc, max_P_Blim, mean_catch, mean_P_closed, mean_years_closed)), n=Inf)
cat("\nFINAL VALIDATION\n")
cat("Expected selected MPs: 12\n")
cat("Observed selected MPs:", nrow(selected_MPs), "\n")
cat("Complete 3 x 4 Fcap-Besc design: PASS\n")
cat("22b/24c reference consistency: PASS\n")
cat("Reference-OM annual precautionarity premise: PASS\n")
cat("No ranking or performance recalculation performed: PASS\n")
cat("\nOutputs saved in:\n", output_dir, "\n", sep="")
cat("\nSCRIPT 25 COMPLETED SUCCESSFULLY\n")
