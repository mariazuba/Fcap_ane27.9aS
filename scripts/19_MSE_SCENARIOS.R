# ============================================================
# 19_MSE_SCENARIOS.R
# Experimental design for Fcap-Besc MSE
#
# STAGE 1
#   Candidate MPs:
#     Fcap x Besc
#
# STAGE 2
#   Performance evaluation:
#     Precautionarity
#     Yield
#     Stability
#
# STAGE 3
#   Robustness tests applied only to selected MPs:
#     Recruitment variability (sigmaR)
#     Seasonal exploitation
#     Assessment uncertainty
#
# Forecast uncertainty is NOT a Stage 3 robustness factor.
# Empirical forecast uncertainty (S2) is an integral component
# of the probabilistic reference HCR evaluated in Stage 2.
# ============================================================

rm(list=ls())

library(dplyr)
library(tidyr)
library(here)

# ============================================================
# 1. FIXED STOCK SETTINGS
# ============================================================

stock <- "ANE"
fleet <- "SEINE"
metier <- "ALL"

Blim <- 4721
Bpa <- 6561

first_projection_year <- 2025
last_projection_year <- 2054

recruitment_age <- "0"
recruitment_season <- "3"
ssb_season <- "2"

risk_threshold <- 0.05

# ============================================================
# 2. CANDIDATE MP PARAMETERS
# ============================================================
#
# Fcap:
#   Range selected from the exploitation diagnostic.
#   Fcap = 1.0 remains close to historical experience.
#   Transition occurs around Fcap = 1.25-1.50.
#   Fcap = 1.75 explores moderately higher exploitation.
#   Fcap = 2.00 defines the upper candidate bound retained after
#   the exploitation diagnostic.
#
# Besc:
#   5000  = close to Blim
#   6561  = Bpa
#   8000  = moderately above Bpa
#   10000 = high escapement
# ============================================================

Fcap_values <- c(1.00,1.25,1.50,1.75,2.00)

Besc_values <- c(5000,6561,8000,10000)

# ============================================================
# 3. REFERENCE OPERATING MODEL
# ============================================================
#
# All candidate MPs are evaluated first under exactly the same
# reference OM assumptions.
#
# Recruitment process error:
#   Historical BH residual variability.
#
# Seasonal exploitation:
#   Current reference seasonal allocation.
#
# Assessment uncertainty:
#   Historical assessment error.
# ============================================================

# Recruitment process error:
#   Reference sigmaR = median across the 100 bootstrap-conditioned OMs.

sigmaR_reference <- 0.216

seasonal_reference <- c(
  Q1=0.138,
  Q2=0.417,
  Q3=0.334,
  Q4=0.111)

assessment_reference <- "historical"

# Forecast uncertainty:
#   S2 empirical contemporary uncertainty from 100 bootstrap-conditioned
#   SS3 assessments. Forecast uncertainty is part of the reference HCR,
#   not a robustness dimension.
#
#   One CV is assigned to each global MSE trajectory and retained through
#   the projection. The same trajectory-CV mapping is used for all MPs.
#
#   Probabilistic criterion:
#     P(SSBforecast < Besc) <= alpha
#   equivalently:
#     Q_alpha(SSBforecast) >= Besc

forecast_uncertainty_reference <- "bootstrap_empirical"
forecast_probability <- TRUE
forecast_alpha <- risk_threshold
forecast_uncertainty_file <- here("data","Rdata","forecast_uncertainty_characterisation.rds")
forecast_CV_seed <- 20260923
n_trajectories <- 1000

stopifnot(file.exists(forecast_uncertainty_file))
forecast_uncertainty <- readRDS(forecast_uncertainty_file)
CV_S2 <- as.numeric(forecast_uncertainty$contemporary$SSB_2025_CV)
CV_S2 <- CV_S2[is.finite(CV_S2) & CV_S2>0]
stopifnot(length(CV_S2)==100)

set.seed(forecast_CV_seed)
forecast_CV_map <- tibble(global_iter=seq_len(n_trajectories),forecast_CV=sample(CV_S2,size=n_trajectories,replace=TRUE))

stopifnot(nrow(forecast_CV_map)==n_trajectories)
stopifnot(!anyDuplicated(forecast_CV_map$global_iter))
stopifnot(all(is.finite(forecast_CV_map$forecast_CV)))
stopifnot(all(forecast_CV_map$forecast_CV>0))

stopifnot(abs(sum(seasonal_reference)-1)<1e-8)
stopifnot(all(seasonal_reference>=0))

# ============================================================
# 4. STAGE 1 - CANDIDATE MP GRID
# ============================================================
#
# 5 Fcap x 4 Besc = 20 candidate MPs.
#
# No robustness dimensions are varied at this stage.
# ============================================================

candidate_MP <- crossing(
  Fcap=Fcap_values,
  Besc=Besc_values)

candidate_MP <- candidate_MP %>% mutate(
  stage="candidate_MP",
  robustness_factor="reference",
  sigmaR_scenario="reference",
  sigmaR=sigmaR_reference,
  seasonal_exploitation="reference",
  Fprop_Q1=seasonal_reference["Q1"],
  Fprop_Q2=seasonal_reference["Q2"],
  Fprop_Q3=seasonal_reference["Q3"],
  Fprop_Q4=seasonal_reference["Q4"],
  assessment_uncertainty=assessment_reference,
  assessment_error=TRUE,
  forecast_uncertainty=forecast_uncertainty_reference,
  forecast_probability=forecast_probability,
  forecast_alpha=forecast_alpha)

# ============================================================
# 5. ADD FIXED STOCK INFORMATION
# ============================================================

candidate_MP <- candidate_MP %>% mutate(
  stock=stock,
  fleet=fleet,
  metier=metier,
  Blim=Blim,
  Bpa=Bpa,
  first_projection_year=first_projection_year,
  last_projection_year=last_projection_year,
  recruitment_age=recruitment_age,
  recruitment_season=recruitment_season,
  ssb_season=ssb_season,
  risk_threshold=risk_threshold)

# ============================================================
# 6. CANDIDATE MP IDENTIFIERS
# ============================================================

candidate_MP <- candidate_MP %>% mutate(
  MP_id=sprintf("Fcap_%0.2f_Besc_%g",Fcap,Besc),
  scenario_id=sprintf("Fcap_%0.2f_Besc_%g_sigmaR_%0.3f_SE_%s_AE_%s",Fcap,Besc,sigmaR,seasonal_exploitation,assessment_uncertainty))

candidate_MP <- candidate_MP %>% mutate(scenario_number=row_number())

print(candidate_MP %>% dplyr::select(scenario_number,MP_id,Fcap,Besc))

stopifnot(identical(candidate_MP$scenario_number,seq_len(nrow(candidate_MP))))

# ============================================================
# 7. CHECK CANDIDATE MP GRID
# ============================================================

stopifnot(nrow(candidate_MP)==length(Fcap_values)*length(Besc_values))
stopifnot(nrow(candidate_MP)==20)
stopifnot(all(candidate_MP$Fcap>=1 & candidate_MP$Fcap<=2))
stopifnot(all(candidate_MP$Besc>Blim))
stopifnot(!anyDuplicated(candidate_MP$MP_id))
stopifnot(!anyDuplicated(candidate_MP$scenario_id))
stopifnot(n_distinct(candidate_MP$sigmaR)==1)
stopifnot(n_distinct(candidate_MP$seasonal_exploitation)==1)
stopifnot(n_distinct(candidate_MP$assessment_uncertainty)==1)
stopifnot(unique(candidate_MP$sigmaR)==sigmaR_reference)
stopifnot(unique(candidate_MP$assessment_uncertainty)==assessment_reference)
stopifnot(n_distinct(candidate_MP$forecast_uncertainty)==1)
stopifnot(all(candidate_MP$forecast_probability))
stopifnot(all(candidate_MP$forecast_alpha==forecast_alpha))


cat("\n============================================\n")
cat("STAGE 1 - CANDIDATE MP EXPERIMENT\n")
cat("============================================\n")

cat("Fcap values:",paste(Fcap_values,collapse=", "),"\n")
cat("Besc values:",paste(Besc_values,collapse=", "),"\n")
cat("Candidate MPs:",nrow(candidate_MP),"\n")
cat("Reference sigmaR:",sigmaR_reference,"\n")
cat("Reference seasonal exploitation:",paste(seasonal_reference,collapse=", "),"\n")
cat("Reference assessment uncertainty:",assessment_reference,"\n")
cat("Reference forecast uncertainty:",forecast_uncertainty_reference,"\n")
cat("Forecast probability alpha:",forecast_alpha,"\n")
cat("Forecast CV source values:",length(CV_S2),"\n")
cat("Forecast CV mapped trajectories:",nrow(forecast_CV_map),"\n")
cat("Forecast CV median:",round(median(forecast_CV_map$forecast_CV),4),"\n")
cat("Forecast CV P05-P95:",round(quantile(forecast_CV_map$forecast_CV,.05),4),"-",round(quantile(forecast_CV_map$forecast_CV,.95),4),"\n\n")

print(candidate_MP %>% dplyr::select(MP_id,Fcap,Besc))

# ============================================================
# 8. ROBUSTNESS - RECRUITMENT VARIABILITY
# ============================================================
#
# Robustness scenarios are defined here but are NOT crossed
# with the complete candidate MP grid.
#
# After Stage 2, they will be applied only to selected MPs.
#
# historical = BH residual variability
# medium     = increased recruitment variability
# high       = high recruitment variability
# ============================================================

sigmaR_scenarios <- tibble(
  sigmaR_scenario=c("reference","medium","high"),
  sigmaR=c(0.216,0.50,0.70),
  reference=c(TRUE,FALSE,FALSE))

# ============================================================
# 9. ROBUSTNESS - SEASONAL EXPLOITATION
# ============================================================
#
# Historical analysis identified three recurrent seasonal
# exploitation patterns.
#
# The exact Q1-Q4 mean proportions for clusters 1-3 should be
# inserted here from the seasonal clustering analysis.
#
# For now only the current reference pattern is activated.
# ============================================================

seasonal_scenarios <- tibble(
  seasonal_exploitation="reference",
  Fprop_Q1=0.138,
  Fprop_Q2=0.417,
  Fprop_Q3=0.334,
  Fprop_Q4=0.111,
  reference=TRUE)

# ============================================================
# 10. ROBUSTNESS - ASSESSMENT UNCERTAINTY
# ============================================================
#
# historical:
#   Historical assessment error is applied to SSB and F.
#
# perfect:
#   No assessment error is applied.
#
# Historical assessment uncertainty is the reference scenario.
# ============================================================

assessment_scenarios <- tibble(
  assessment_uncertainty=c("historical","perfect"),
  assessment_error=c(TRUE,FALSE),
  reference=c(TRUE,FALSE))
# ============================================================
# 11. CHECK ROBUSTNESS DEFINITIONS
# ============================================================

seasonal_total <- seasonal_scenarios %>% transmute(total=Fprop_Q1+Fprop_Q2+Fprop_Q3+Fprop_Q4)

stopifnot(all(abs(seasonal_total$total-1)<1e-8))
stopifnot(all(sigmaR_scenarios$sigmaR>0))
stopifnot(sum(sigmaR_scenarios$reference)==1)
stopifnot(sum(seasonal_scenarios$reference)==1)
stopifnot(sum(assessment_scenarios$reference)==1)
stopifnot(!anyDuplicated(sigmaR_scenarios$sigmaR_scenario))
stopifnot(!anyDuplicated(seasonal_scenarios$seasonal_exploitation))
stopifnot(!anyDuplicated(assessment_scenarios$assessment_uncertainty))
stopifnot(identical(candidate_MP$scenario_number,seq_len(nrow(candidate_MP))))
stopifnot(n_distinct(candidate_MP$sigmaR)==1)
stopifnot(n_distinct(candidate_MP$seasonal_exploitation)==1)
stopifnot(n_distinct(candidate_MP$assessment_uncertainty)==1)
# ============================================================
# 12. ROBUSTNESS CATALOGUE
# ============================================================
#
# These tables define available robustness tests.
# They do NOT generate additional MSE runs.
# ============================================================

cat("\n============================================\n")
cat("STAGE 3 - ROBUSTNESS CATALOGUE\n")
cat("============================================\n")

cat("\nRecruitment variability:\n")
print(sigmaR_scenarios)

cat("\nSeasonal exploitation:\n")
print(seasonal_scenarios)

cat("\nAssessment uncertainty:\n")
print(assessment_scenarios)

# ============================================================
# 13. SAVE EXPERIMENTAL DESIGN
# ============================================================

scenario_dir <- here("data","mse","scenarios")
dir.create(scenario_dir,recursive=TRUE,showWarnings=FALSE)

forecast_CV_map_file <- file.path(scenario_dir,"forecast_CV_trajectory_map.csv")
write.csv(forecast_CV_map,forecast_CV_map_file,row.names=FALSE)
write.csv(candidate_MP,file.path(scenario_dir,"candidate_MP_grid.csv"),row.names=FALSE)
write.csv(sigmaR_scenarios,file.path(scenario_dir,"robustness_sigmaR.csv"),row.names=FALSE)
write.csv(seasonal_scenarios,file.path(scenario_dir,"robustness_seasonal_exploitation.csv"),row.names=FALSE)
write.csv(assessment_scenarios,file.path(scenario_dir,"robustness_assessment_uncertainty.csv"),row.names=FALSE)

experimental_design <- list(
  candidate_MP=candidate_MP,
  forecast_CV_map=forecast_CV_map,
  forecast_CV_seed=forecast_CV_seed,
  forecast_uncertainty_reference=forecast_uncertainty_reference,
  forecast_alpha=forecast_alpha,
  sigmaR_scenarios=sigmaR_scenarios,
  seasonal_scenarios=seasonal_scenarios,
  assessment_scenarios=assessment_scenarios)

saveRDS(experimental_design,file.path(scenario_dir,"MSE_experimental_design.rds"))

# ============================================================
# 14. SUMMARY
# ============================================================

cat("\n============================================\n")
cat("MSE EXPERIMENTAL DESIGN\n")
cat("============================================\n")
cat("Stage 1 candidate MPs:",nrow(candidate_MP),"\n")
cat("Reference forecast uncertainty:",forecast_uncertainty_reference,"\n")
cat("Forecast CV trajectory map:",forecast_CV_map_file,"\n")
cat("Stage 2: evaluate precautionarity, yield and stability\n")
cat("Stage 3: apply targeted robustness tests to selected MPs\n")
cat("\nExperimental design saved in:",scenario_dir,"\n")

