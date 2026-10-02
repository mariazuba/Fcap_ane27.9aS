# ============================================================
# 19b_REFINE_CANDIDATE_MP_GRID.R
# Refined Fcap-Besc candidate grid after Reference-OM diagnostics
#
# Stock: European anchovy ane.27.9aS (Gulf of Cadiz)
#
# Purpose:
#   Build a refined Reference-OM candidate grid while preserving the
#   original exploratory grid produced by 19_MSE_SCENARIOS.R.
#
# Exploratory grid:
#   Fcap = 1.00, 1.25, 1.50, 1.75, 2.00
#   Besc = 5000, 6561, 8000, 10000
#
# Refined grid:
#   Fcap = 1.00, 1.25, 1.50, 1.60, 1.70, 1.75, 1.80
#   Besc = 5000, 6561, 7250, 8000
#
# Expected:
#   28 refined MPs = 12 existing + 16 new.
#
# This script does NOT run FLBEIA or recalculate performance.
# ============================================================

rm(list=ls())

library(dplyr)
library(tidyr)
library(here)

# ============================================================
# 1. PATHS
# ============================================================

scenario_dir     <- here("data","mse","scenarios")
exploratory_file <- file.path(scenario_dir,"candidate_MP_grid.csv")
forecast_CV_file <- file.path(scenario_dir,"forecast_CV_trajectory_map.csv")
output_dir       <- file.path(scenario_dir,"refined")

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

stopifnot(file.exists(exploratory_file))
stopifnot(file.exists(forecast_CV_file))

# ============================================================
# 2. LOAD AND VALIDATE EXPLORATORY DESIGN
# ============================================================

exploratory_MP  <- read.csv(exploratory_file,stringsAsFactors=FALSE)
forecast_CV_map <- read.csv(forecast_CV_file,stringsAsFactors=FALSE)

required_columns <- c(
  "scenario_number","MP_id","scenario_id","Fcap","Besc","stage",
  "robustness_factor","sigmaR_scenario","sigmaR","seasonal_exploitation",
  "Fprop_Q1","Fprop_Q2","Fprop_Q3","Fprop_Q4",
  "assessment_uncertainty","assessment_error","forecast_uncertainty",
  "forecast_probability","forecast_alpha","stock","fleet","metier",
  "Blim","Bpa","first_projection_year","last_projection_year",
  "recruitment_age","recruitment_season","ssb_season","risk_threshold")

missing_columns <- setdiff(required_columns,names(exploratory_MP))

stopifnot(length(missing_columns)==0)
stopifnot(nrow(exploratory_MP)==20)
stopifnot(!anyDuplicated(exploratory_MP$scenario_number))
stopifnot(!anyDuplicated(exploratory_MP$MP_id))
stopifnot(!anyDuplicated(exploratory_MP$scenario_id))
stopifnot(!anyDuplicated(exploratory_MP[c("Fcap","Besc")]))

stopifnot(all(c("global_iter","forecast_CV") %in% names(forecast_CV_map)))
stopifnot(nrow(forecast_CV_map)==1000)
stopifnot(identical(forecast_CV_map$global_iter,seq_len(1000)))
stopifnot(!anyDuplicated(forecast_CV_map$global_iter))
stopifnot(all(is.finite(forecast_CV_map$forecast_CV)))
stopifnot(all(forecast_CV_map$forecast_CV>0))

# ============================================================
# 3. REFERENCE-OM SETTINGS
# ============================================================

reference_columns <- c(
  "stage","robustness_factor","sigmaR_scenario","sigmaR",
  "seasonal_exploitation","Fprop_Q1","Fprop_Q2","Fprop_Q3","Fprop_Q4",
  "assessment_uncertainty","assessment_error","forecast_uncertainty",
  "forecast_probability","forecast_alpha","stock","fleet","metier",
  "Blim","Bpa","first_projection_year","last_projection_year",
  "recruitment_age","recruitment_season","ssb_season","risk_threshold")

reference_settings <- exploratory_MP %>% dplyr::select(all_of(reference_columns)) %>%distinct()

stopifnot(nrow(reference_settings)==1)

# ============================================================
# 4. DEFINE REFINED GRID
# ============================================================

Fcap_values_refined <- c(1.00,1.25,1.50,1.60,1.70,1.75,1.80)
Besc_values_refined <- c(5000,6561,7250,8000)

expected_refined <- 28
expected_existing <- 12
expected_new <- 16

refined_MP <- crossing(Fcap=Fcap_values_refined,Besc=Besc_values_refined)

refined_MP <- refined_MP %>%
              mutate(refined_scenario_number=row_number()) %>%
              crossing(reference_settings) %>%
              mutate(MP_id=sprintf("Fcap_%0.2f_Besc_%g",Fcap,Besc),
                     scenario_id=sprintf("Fcap_%0.2f_Besc_%g_sigmaR_%0.3f_SE_%s_AE_%s",
                                 Fcap,Besc,sigmaR,seasonal_exploitation,assessment_uncertainty))

stopifnot(nrow(refined_MP)==expected_refined)
stopifnot(!anyDuplicated(refined_MP[c("Fcap","Besc")]))
stopifnot(!anyDuplicated(refined_MP$MP_id))
stopifnot(!anyDuplicated(refined_MP$scenario_id))

# ============================================================
# 5. IDENTIFY EXISTING AND NEW MPs
# ============================================================

exploratory_lookup <- exploratory_MP %>%
                      dplyr::select(Fcap,Besc,exploratory_scenario_number=scenario_number,
                      exploratory_MP_id=MP_id,exploratory_scenario_id=scenario_id)

refined_MP <- refined_MP %>%
              left_join(exploratory_lookup,by=c("Fcap","Besc")) %>%
              mutate(run_status=ifelse(is.na(exploratory_scenario_number),"new","existing"),
                     requires_run=run_status=="new")

stopifnot(sum(refined_MP$run_status=="existing")==expected_existing)
stopifnot(sum(refined_MP$run_status=="new")==expected_new)

# ============================================================
# 6. VERIFY EXPECTED EXISTING AND NEW COMBINATIONS
# ============================================================

expected_existing_grid <- crossing(
                          Fcap=c(1.00,1.25,1.50,1.75),
                          Besc=c(5000,6561,8000)) %>%
                          arrange(Fcap,Besc)

existing_grid <- refined_MP %>%
                  filter(run_status=="existing") %>%
                  dplyr::select(Fcap,Besc) %>%
                  arrange(Fcap,Besc)

stopifnot(identical(existing_grid,expected_existing_grid))

new_from_Fcap <- refined_MP %>% filter(Fcap %in% c(1.60,1.70,1.80),Besc %in% c(5000,6561,8000))

new_from_Besc <- refined_MP %>% filter(Besc==7250)

stopifnot(nrow(new_from_Fcap)==9)
stopifnot(nrow(new_from_Besc)==7)
stopifnot(nrow(bind_rows(new_from_Fcap,new_from_Besc) %>% distinct(Fcap,Besc))==expected_new)

# ============================================================
# 7. BUILD RUNNER TABLE FOR THE 16 NEW MPs
# ============================================================
#
# Script 20 selects rows by scenario_number. The runner table therefore
# uses a compact scenario_number 1:16. refined_scenario_number retains
# the position of each MP in the complete 28-row refined grid.
# ============================================================

new_MP_RUN <- refined_MP %>%
              filter(requires_run) %>%
              arrange(refined_scenario_number) %>%
              mutate(run_scenario_number=row_number(),
                     scenario_number=run_scenario_number)

stopifnot(nrow(new_MP_RUN)==expected_new)
stopifnot(identical(new_MP_RUN$scenario_number,seq_len(expected_new)))
stopifnot(!anyDuplicated(new_MP_RUN$scenario_id))

# ============================================================
# 8. CHECK COMPATIBILITY WITH SCRIPT 20
# ============================================================

runner_required_columns <- c(
  "scenario_number","scenario_id","Fcap","Besc","sigmaR_scenario","sigmaR",
  "seasonal_exploitation","Fprop_Q1","Fprop_Q2","Fprop_Q3","Fprop_Q4",
  "assessment_uncertainty","assessment_error","forecast_uncertainty",
  "forecast_probability","forecast_alpha")

runner_missing_columns <- setdiff(runner_required_columns,names(new_MP_RUN))

stopifnot(length(runner_missing_columns)==0)
stopifnot(all(new_MP_RUN$Fcap>0))
stopifnot(all(new_MP_RUN$Besc>0))
stopifnot(all(new_MP_RUN$sigmaR>0))
stopifnot(all(new_MP_RUN$assessment_error))
stopifnot(all(new_MP_RUN$forecast_probability))
stopifnot(all(new_MP_RUN$forecast_alpha>0 & new_MP_RUN$forecast_alpha<1))

seasonal_total <- new_MP_RUN$Fprop_Q1+new_MP_RUN$Fprop_Q2+new_MP_RUN$Fprop_Q3+new_MP_RUN$Fprop_Q4

stopifnot(all(abs(seasonal_total-1)<1e-8))

# ============================================================
# 9. ORGANIZE OUTPUT TABLES
# ============================================================

refined_MP <- refined_MP %>%
              dplyr::select(
              refined_scenario_number,exploratory_scenario_number,run_status,requires_run,
              MP_id,scenario_id,Fcap,Besc,all_of(reference_columns),
              exploratory_MP_id,exploratory_scenario_id)

existing_MP <- refined_MP %>% filter(run_status=="existing") %>% arrange(refined_scenario_number)

new_MP <- refined_MP %>% filter(run_status=="new") %>% arrange(refined_scenario_number)

new_MP_RUN <- new_MP_RUN %>%
              dplyr::select(
              scenario_number,run_scenario_number,refined_scenario_number,
              exploratory_scenario_number,run_status,requires_run,
              MP_id,scenario_id,Fcap,Besc,all_of(reference_columns),
              exploratory_MP_id,exploratory_scenario_id)

# ============================================================
# 10. SAVE REFINED DESIGN
# ============================================================

refined_file  <- file.path(output_dir,"candidate_MP_grid_refined.csv")
existing_file <- file.path(output_dir,"candidate_MP_grid_refined_existing.csv")
new_file      <- file.path(output_dir,"candidate_MP_grid_refined_new.csv")
new_RUN_file  <- file.path(output_dir,"candidate_MP_grid_refined_new_RUN.csv")
summary_file  <- file.path(output_dir,"refined_MP_grid_summary.rds")

write.csv(refined_MP,refined_file,row.names=FALSE)
write.csv(existing_MP,existing_file,row.names=FALSE)
write.csv(new_MP,new_file,row.names=FALSE)
write.csv(new_MP_RUN,new_RUN_file,row.names=FALSE)

refined_summary <- list(
                    exploratory_file=exploratory_file,
                    forecast_CV_file=forecast_CV_file,
                    Fcap_values_refined=Fcap_values_refined,
                    Besc_values_refined=Besc_values_refined,
                    n_refined=expected_refined,
                    n_existing=expected_existing,
                    n_new=expected_new,
                    refined_MP=refined_MP,
                    existing_MP=existing_MP,
                    new_MP=new_MP,
                    new_MP_RUN=new_MP_RUN)

saveRDS(refined_summary,summary_file)

# ============================================================
# 11. REPORT AND FINAL VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("REFINED CANDIDATE MP GRID\n")
cat("============================================================\n")
cat("Fcap values:",paste(Fcap_values_refined,collapse=", "),"\n")
cat("Besc values:",paste(Besc_values_refined,collapse=", "),"\n")
cat("Refined candidate MPs:",nrow(refined_MP),"\n")
cat("Existing combinations:",nrow(existing_MP),"\n")
cat("New combinations:",nrow(new_MP),"\n")
cat("\nNew combinations requiring MSE execution:\n")

print(new_MP_RUN %>% dplyr::select(scenario_number,refined_scenario_number,Fcap,Besc,scenario_id),n=Inf)

stopifnot(file.exists(refined_file))
stopifnot(file.exists(existing_file))
stopifnot(file.exists(new_file))
stopifnot(file.exists(new_RUN_file))
stopifnot(file.exists(summary_file))
stopifnot(nrow(refined_MP)==28)
stopifnot(nrow(existing_MP)==12)
stopifnot(nrow(new_MP)==16)
stopifnot(nrow(new_MP_RUN)==16)
stopifnot(!any(refined_MP$Fcap==2.00))
stopifnot(!any(refined_MP$Besc==10000))
stopifnot(all(c(1.60,1.70,1.80) %in% refined_MP$Fcap))
stopifnot(7250 %in% refined_MP$Besc)

cat("\n============================================================\n")
cat("REFINED MP GRID COMPLETED\n")
cat("============================================================\n")
cat("Exploratory grid preserved:",exploratory_file,"\n")
cat("Refined grid:",refined_file,"\n")
cat("Runner table for new MPs:",new_RUN_file,"\n")
cat("Forecast CV map preserved:",forecast_CV_file,"\n")
cat("\nNo FLBEIA simulations were run.\n")
cat("Only the 16 rows in candidate_MP_grid_refined_new_RUN.csv require new MSE runs.\n")
