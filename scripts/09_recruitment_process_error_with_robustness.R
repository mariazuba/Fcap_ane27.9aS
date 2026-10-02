# ============================================================
# 09_recruitment_process_error.R
# Generate future recruitment process error
#
# Stock / analysis:
#   European anchovy, ane.27.9aS (Gulf of Cádiz)
#
# Role in the workflow:
#   This script follows the expansion of the 100 conditioned OMs
#   to 1000 trajectories in 08_expand_OM_100_to_1000.R.
#   It introduces stochastic recruitment process error during the
#   projection period while preserving the conditioned OM-specific
#   Beverton-Holt relationships.
#
# Purpose:
#   - Load the 1000-iteration conditioned stock-recruitment object.
#   - Recover the historical sigmaR distribution across the 100
#     conditioned OMs.
#   - Implement a reference process-error scenario using the median
#     sigmaR across OMs.
#   - Implement a robustness scenario in which each trajectory inherits
#     the sigmaR estimated for its source OM.
#   - Generate paired bias-corrected lognormal recruitment multipliers
#     using the same standard-normal deviates in both scenarios.
#   - Assign recruitment process error to Q3 in FLSRsim during the
#     projection period.
#   - Validate the generated deviations, multipliers and their assignment.
#   - Save diagnostic tables and figures used in the Working Document.
#
# Inputs:
#   data/Rdata/SRs_reference_conditioned_1000iter.rds
#   data/Rdata/iteration_map_1000.csv
#
# Outputs — reference scenario:
#   data/Rdata/SRs_reference_1000iter_process_error.rds
#   data/Rdata/recruitment_deviations_reference_1000iter.rds
#   data/Rdata/recruitment_multipliers_reference_1000iter.rds
#
# Outputs — OM-specific sigmaR robustness scenario:
#   data/Rdata/SRs_om_specific_1000iter_process_error.rds
#   data/Rdata/recruitment_deviations_om_specific_1000iter.rds
#   data/Rdata/recruitment_multipliers_om_specific_1000iter.rds
#
# Common stochastic component:
#   data/Rdata/recruitment_standard_normal_deviates_1000iter.rds
#
# Working Document outputs:
#   outputs/recruitment/tables/sigmaR_by_OM.csv
#   outputs/recruitment/tables/sigmaR_historical_summary.csv
#   outputs/recruitment/tables/recruitment_process_error_scenarios_summary.csv
#   outputs/recruitment/tables/sigmaR_om_specific_assignment_check.csv
#   outputs/recruitment/figures/recruitment_process_error_scenarios.png
#   outputs/recruitment/figures/recruitment_process_error_scenarios.pdf
#
# Projection period:
#   2025-2054
#
# Simulation structure:
#   Historical OMs:       100
#   Replicates per OM:     10
#   Total trajectories:  1000
#
# Recruitment process error:
#   Recruitment season: Q3
#   Distribution: bias-corrected lognormal
#   Reference scenario: median historical sigmaR across the 100 OMs
#   Robustness scenario: OM-specific historical sigmaR
#   Random seed: 1234
#
# Important:
#   The common-median sigmaR formulation is retained as the reference
#   process-error scenario.
#   The OM-specific sigmaR formulation is implemented as a robustness
#   scenario and does not replace the reference formulation.
#   Both scenarios use the same matrix of standard-normal deviates,
#   allowing paired comparison of the alternative sigmaR parameterisations.
#   Recruitment deviations are generated independently across projection
#   years and trajectories before scaling by the corresponding sigmaR.
#   The lognormal bias correction gives a theoretical expected multiplier
#   of 1; the realised finite-sample mean need not equal exactly 1.
#   Process-error multipliers are assigned only to Q3 during 2025-2054;
#   uncertainty multipliers in the remaining seasons are fixed at 1.
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================

rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)
library(ggplot2)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

fig_dir <- file.path(here(),"outputs","recruitment","figures")
tab_dir <- file.path(here(),"outputs","recruitment","tables")

dir.create(fig_dir,recursive = TRUE,showWarnings = FALSE)
dir.create(tab_dir,recursive = TRUE,showWarnings = FALSE)

rds_boot_dir <- file.path(here(),"data","Rdata")

first.proj.yr <- 2025
last.proj.yr <- 2054
ss.rec <- 3
n_om <- 100
n_it <- 1000
seed_R <- 1234

# ============================================================
# 2. LOAD INPUT OBJECTS
# ============================================================
# Load the expanded stock-recruitment object and the explicit mapping
# between the 1000 trajectories and their 100 source OMs.

SRs <- readRDS(file.path(rds_boot_dir,"SRs_reference_conditioned_1000iter.rds"))
expansion_map <- read.csv(file.path(rds_boot_dir,"iteration_map_1000.csv"))

# ============================================================
# 3. VALIDATE INPUT STRUCTURE
# ============================================================
# Confirm the expected stock, mapping fields, number of trajectories
# and number of source OMs before generating stochastic deviations.

stopifnot("ANE" %in% names(SRs))
stopifnot(all(c("iter1000","om","sigmaR") %in% names(expansion_map)))
stopifnot(nrow(expansion_map) == n_it)
stopifnot(dplyr::n_distinct(expansion_map$iter1000) == n_it)
stopifnot(dplyr::n_distinct(expansion_map$om) == n_om)

# ============================================================
# 4. DERIVE REFERENCE RECRUITMENT PROCESS ERROR
# ============================================================
# Recover one historical residual SD for each source OM and use the
# median across the 100 conditioned OMs as the common projection sigmaR.

sigmaR_boot <- expansion_map %>%
  dplyr::distinct(om,sigmaR)

stopifnot(nrow(sigmaR_boot) == n_om)
stopifnot(all(is.finite(sigmaR_boot$sigmaR)))
stopifnot(all(sigmaR_boot$sigmaR > 0))

sigmaR_ref <- median(sigmaR_boot$sigmaR)

cat("\n============================================\n")
cat("REFERENCE RECRUITMENT PROCESS ERROR\n")
cat("============================================\n")
cat("Bootstrap OMs =",nrow(sigmaR_boot),"\n")
cat("sigmaR range =",range(sigmaR_boot$sigmaR),"\n")
cat("sigmaR median =",sigmaR_ref,"\n")
cat("============================================\n")

# ============================================================
# 5. DEFINE sigmaR FOR BOTH SCENARIOS
# ============================================================
# Reference: common median sigmaR.
# Robustness: each trajectory inherits sigmaR from its source OM.

sigmaR_iter_reference <- rep(sigmaR_ref,n_it)
sigmaR_iter_om_specific <- expansion_map$sigmaR

stopifnot(length(sigmaR_iter_reference) == n_it)
stopifnot(length(sigmaR_iter_om_specific) == n_it)
stopifnot(all(is.finite(sigmaR_iter_om_specific)))
stopifnot(all(sigmaR_iter_om_specific > 0))

# ============================================================
# 6. GENERATE PAIRED STANDARD-NORMAL DEVIATES
# ============================================================
# The same standard-normal matrix is used for both scenarios so that
# differences between them arise from sigmaR parameterisation.

set.seed(seed_R)

proj_yrs <- as.character(first.proj.yr:last.proj.yr)

z_R <- matrix(
  rnorm(length(proj_yrs) * n_it),
  nrow = length(proj_yrs),
  ncol = n_it,
  dimnames = list(year = proj_yrs,iter = as.character(seq_len(n_it)))
)

stopifnot(all(is.finite(z_R)))

# ============================================================
# 7. GENERATE REFERENCE SCENARIO — COMMON MEDIAN sigmaR
# ============================================================

eps_R_reference <- sweep(z_R,2,sigmaR_iter_reference,"*")
log_mult_R_reference <- sweep(eps_R_reference,2,0.5 * sigmaR_iter_reference^2,"-")
mult_R_reference <- exp(log_mult_R_reference)

# ============================================================
# 8. GENERATE ROBUSTNESS SCENARIO — OM-SPECIFIC sigmaR
# ============================================================

eps_R_om_specific <- sweep(z_R,2,sigmaR_iter_om_specific,"*")
log_mult_R_om_specific <- sweep(eps_R_om_specific,2,0.5 * sigmaR_iter_om_specific^2,"-")
mult_R_om_specific <- exp(log_mult_R_om_specific)

# ============================================================
# 9. VALIDATE GENERATED PROCESS ERROR
# ============================================================

sigmaR_realised_reference <- sd(as.numeric(eps_R_reference))
sigmaR_realised_by_iter <- apply(eps_R_om_specific,2,sd)

stopifnot(all(is.finite(eps_R_reference)))
stopifnot(all(is.finite(mult_R_reference)))
stopifnot(all(mult_R_reference > 0))
stopifnot(abs(sigmaR_realised_reference - sigmaR_ref) < 0.01)
stopifnot(abs(sd(log(as.numeric(mult_R_reference))) - sigmaR_realised_reference) < 1e-10)

stopifnot(all(is.finite(eps_R_om_specific)))
stopifnot(all(is.finite(mult_R_om_specific)))
stopifnot(all(mult_R_om_specific > 0))
stopifnot(all(is.finite(sigmaR_realised_by_iter)))

cat("\nReference scenario:\n")
cat("Target sigmaR =",sigmaR_ref,"\n")
cat("Realised sigmaR =",sigmaR_realised_reference,"\n")
cat("Mean multiplier =",mean(mult_R_reference),"\n")
cat("Median multiplier =",median(mult_R_reference),"\n")

cat("\nOM-specific robustness scenario:\n")
cat("Target sigmaR range =",range(sigmaR_iter_om_specific),"\n")
cat("Realised per-trajectory sigmaR range =",range(sigmaR_realised_by_iter),"\n")
cat("Mean multiplier =",mean(mult_R_om_specific),"\n")
cat("Median multiplier =",median(mult_R_om_specific),"\n")

# ============================================================
# 10. ASSIGN PROCESS ERROR TO FLSRsim
# ============================================================

SRs_reference <- SRs
SRs_om_specific <- SRs

SRs_reference$ANE@uncertainty[] <- 1
SRs_om_specific$ANE@uncertainty[] <- 1

for(i in seq_len(n_it)) {
  SRs_reference$ANE@uncertainty[,proj_yrs,,ac(ss.rec),,i] <- mult_R_reference[,i]
}

for(i in seq_len(n_it)) {
  SRs_om_specific$ANE@uncertainty[,proj_yrs,,ac(ss.rec),,i] <- mult_R_om_specific[,i]
}

# ============================================================
# 11. VALIDATE FLSRsim ASSIGNMENT
# ============================================================

unc_Q3_reference <- as.numeric(SRs_reference$ANE@uncertainty[,proj_yrs,,ac(ss.rec),,])
unc_Q3_om_specific <- as.numeric(SRs_om_specific$ANE@uncertainty[,proj_yrs,,ac(ss.rec),,])

stopifnot(isTRUE(all.equal(unc_Q3_reference,as.numeric(mult_R_reference),tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(unc_Q3_om_specific,as.numeric(mult_R_om_specific),tolerance = 1e-12)))

other_seasons <- setdiff(seq_len(4),ss.rec)

for(s in other_seasons) {
  stopifnot(all(as.numeric(SRs_reference$ANE@uncertainty[,proj_yrs,,ac(s),,]) == 1))
}

for(s in other_seasons) {
  stopifnot(all(as.numeric(SRs_om_specific$ANE@uncertainty[,proj_yrs,,ac(s),,]) == 1))
}

cat("\nFLSRsim process-error assignment: OK for both scenarios\n")

# ============================================================
# 12. WD TABLES
# ============================================================

sigmaR_om_table <- sigmaR_boot %>%
  dplyr::rename(sigmaR_historical = sigmaR)

sigmaR_summary <- data.frame(
  n_om = n_om,
  minimum = min(sigmaR_boot$sigmaR),
  p05 = as.numeric(quantile(sigmaR_boot$sigmaR,0.05)),
  median = median(sigmaR_boot$sigmaR),
  mean = mean(sigmaR_boot$sigmaR),
  p95 = as.numeric(quantile(sigmaR_boot$sigmaR,0.95)),
  maximum = max(sigmaR_boot$sigmaR)
)

process_error_scenario_summary <- data.frame(
  scenario = c("common_median","om_specific"),
  n_trajectories = c(n_it,n_it),
  sigmaR_parameterisation = c("Median across 100 OMs","OM-specific"),
  sigmaR_min = c(sigmaR_ref,min(sigmaR_iter_om_specific)),
  sigmaR_median = c(sigmaR_ref,median(sigmaR_iter_om_specific)),
  sigmaR_max = c(sigmaR_ref,max(sigmaR_iter_om_specific)),
  mean_multiplier = c(mean(mult_R_reference),mean(mult_R_om_specific)),
  median_multiplier = c(median(mult_R_reference),median(mult_R_om_specific)),
  multiplier_min = c(min(mult_R_reference),min(mult_R_om_specific)),
  multiplier_max = c(max(mult_R_reference),max(mult_R_om_specific))
)

sigmaR_assignment_check <- data.frame(
  iter1000 = expansion_map$iter1000,
  om = expansion_map$om,
  sigmaR_target = sigmaR_iter_om_specific,
  sigmaR_realised = sigmaR_realised_by_iter
)

write.csv(sigmaR_om_table,file.path(tab_dir,"sigmaR_by_OM.csv"),row.names = FALSE)
write.csv(sigmaR_summary,file.path(tab_dir,"sigmaR_historical_summary.csv"),row.names = FALSE)
write.csv(process_error_scenario_summary,file.path(tab_dir,"recruitment_process_error_scenarios_summary.csv"),row.names = FALSE)
write.csv(sigmaR_assignment_check,file.path(tab_dir,"sigmaR_om_specific_assignment_check.csv"),row.names = FALSE)

# ============================================================
# 13. WD FIGURE
# ============================================================

mult_R_df <- rbind(
  data.frame(scenario = "Common median sigmaR",multiplier = as.numeric(mult_R_reference)),
  data.frame(scenario = "OM-specific sigmaR",multiplier = as.numeric(mult_R_om_specific))
)

p_process_error <- ggplot(mult_R_df,aes(x = multiplier)) +
  geom_histogram(bins = 60,boundary = 1) +
  geom_vline(xintercept = 1,linetype = "dashed",linewidth = 0.8) +
  facet_wrap(~scenario,ncol = 1) +
  labs(
    x = "Recruitment process-error multiplier",
    y = "Frequency",
    title = "Simulated future recruitment process error"
  ) +
  theme_bw()

ggsave(file.path(fig_dir,"recruitment_process_error_scenarios.png"),p_process_error,width = 8,height = 8,dpi = 300)
ggsave(file.path(fig_dir,"recruitment_process_error_scenarios.pdf"),p_process_error,width = 8,height = 8)

# ============================================================
# 14. FINAL IMPLEMENTATION SUMMARY
# ============================================================

cat("\n============================================\n")
cat("RECRUITMENT PROCESS ERROR COMPLETE\n")
cat("============================================\n")
cat("Projection years:",first.proj.yr,"-",last.proj.yr,"\n")
cat("Iterations:",n_it,"\n")
cat("Reference sigmaR:",round(sigmaR_ref,6),"\n")
cat("Reference realised sigmaR:",round(sigmaR_realised_reference,6),"\n")
cat("Reference mean multiplier:",round(mean(mult_R_reference),6),"\n")
cat("OM-specific sigmaR range:",range(sigmaR_iter_om_specific),"\n")
cat("OM-specific mean multiplier:",round(mean(mult_R_om_specific),6),"\n")
cat("Recruitment season:",ss.rec,"\n")
cat("============================================\n")

# ============================================================
# 15. SAVE OUTPUT OBJECTS
# ============================================================

saveRDS(SRs_reference,file.path(rds_boot_dir,"SRs_reference_1000iter_process_error.rds"))
saveRDS(eps_R_reference,file.path(rds_boot_dir,"recruitment_deviations_reference_1000iter.rds"))
saveRDS(mult_R_reference,file.path(rds_boot_dir,"recruitment_multipliers_reference_1000iter.rds"))

saveRDS(SRs_om_specific,file.path(rds_boot_dir,"SRs_om_specific_1000iter_process_error.rds"))
saveRDS(eps_R_om_specific,file.path(rds_boot_dir,"recruitment_deviations_om_specific_1000iter.rds"))
saveRDS(mult_R_om_specific,file.path(rds_boot_dir,"recruitment_multipliers_om_specific_1000iter.rds"))

saveRDS(z_R,file.path(rds_boot_dir,"recruitment_standard_normal_deviates_1000iter.rds"))
