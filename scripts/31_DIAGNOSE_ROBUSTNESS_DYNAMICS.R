# ============================================================
# 31_DIAGNOSE_ROBUSTNESS_DYNAMICS.R
# Diagnose variation hidden by aggregate robustness performance
# Anchovy ane.27.9aS: R-OM, S-C1, S-C2, I0 or I1 versus REF
#
# Uses paired trajectory metrics and annual ensemble metrics from
# script 30, and implementation flags from script 28. No MSE reruns.
# Descriptive diagnostics: distributions of paired differences,
# mean differences per conditioned OM (10 replicates), largest
# trajectory changes, annual risk/closure/catch changes, and TAC flags.
# These tables locate differences; they do not establish their cause.
# No raw seasonal extraction, recruitment mechanism analysis,
# significance tests, ranking or figures are performed here.
#
# Usage (R console, project root):
#   source("scripts/31_DIAGNOSE_ROBUSTNESS_DYNAMICS.R") # default R-OM
#   Sys.setenv(MSE_DIAGNOSE_ROBUSTNESS_ID="S-C1")
#   source("scripts/31_DIAGNOSE_ROBUSTNESS_DYNAMICS.R")
# Inputs: script-30 comparison RDS and script-28 validation RDS.
# Outputs: outputs/mse/robustness/diagnostics/<scenario>/
#          CSV tables and 31_robustness_diagnostics.rds.
# All differences are robustness minus REF. Probability differences
# in percentage-point columns are multiplied by 100.
# Quantiles describe paired differences, not confidence intervals.
# Re-running replaces only these diagnostic outputs.
# ============================================================

rm(list=ls())
library(dplyr)
library(tidyr)
library(here)

robustness_id <- Sys.getenv("MSE_DIAGNOSE_ROBUSTNESS_ID", "R-OM")
stopifnot(robustness_id %in% c("R-OM", "S-C1", "S-C2", "I0", "I1"))
comparison_file <- here("outputs", "mse", "robustness", "comparison", robustness_id, "30_comparison_with_REF.rds")
validation_file <- here("outputs", "mse", "robustness", "validation", "post_run", robustness_id, "21c_refined_new_validation_summary.rds")
output_dir <- here("outputs", "mse", "robustness", "diagnostics", robustness_id)
stopifnot(file.exists(comparison_file), file.exists(validation_file))
comparison <- readRDS(comparison_file)
validation <- readRDS(validation_file)
stopifnot(identical(comparison$robustness_id, robustness_id))
stopifnot(identical(validation$settings$robustness_id, robustness_id))
stopifnot(all(validation$scenario_validation$structural_PASS))
trajectory <- comparison$trajectory
annual <- comparison$annual
stopifnot(nrow(trajectory)==48000, nrow(annual)==720)
stopifnot(!anyDuplicated(trajectory[c("reference_scenario_id","horizon","iter1000")]))

# Compare distributions of paired changes within each horizon.
metrics <- c("mean_catch", "mean_SSB", "min_SSB", "IAV", "n_years_closed")
delta_cols <- paste0(metrics, "_delta")
stopifnot(all(delta_cols %in% names(trajectory)))
long <- trajectory %>%
  dplyr::select(robustness_id, reference_scenario_number, reference_scenario_id,
                Fcap, Besc, horizon, iter1000, om, replicate, all_of(delta_cols)) %>%
  pivot_longer(all_of(delta_cols), names_to="metric", values_to="delta") %>%
  mutate(metric=sub("_delta$", "", metric))
stopifnot(all(is.finite(long$delta)))
trajectory_summary <- long %>%
  group_by(robustness_id, reference_scenario_number, reference_scenario_id, Fcap, Besc, horizon, metric) %>%
  summarise(n_trajectories=n(), mean_delta=mean(delta), median_delta=median(delta),
            p05_delta=quantile(delta,0.05), p95_delta=quantile(delta,0.95),
            min_delta=min(delta), max_delta=max(delta),
            fraction_positive=mean(delta>0), fraction_negative=mean(delta<0),
            fraction_zero=mean(delta==0), .groups="drop")
stopifnot(all(trajectory_summary$n_trajectories==1000))

# Replicates share their historical OM; report their mean change per OM.
om_summary <- long %>%
  group_by(robustness_id, reference_scenario_number, reference_scenario_id, Fcap, Besc, horizon, metric, om) %>%
  summarise(n_replicates=n(), mean_delta=mean(delta), median_delta=median(delta),
            min_delta=min(delta), max_delta=max(delta), .groups="drop")
stopifnot(all(om_summary$n_replicates==10))

# Five largest absolute changes per MP, horizon and metric.
extreme_trajectories <- long %>%
  group_by(reference_scenario_id, horizon, metric) %>%
  slice_max(order_by=abs(delta), n=5, with_ties=FALSE) %>% ungroup()

# Annual ensemble differences: Short/Medium/Long duplicate Full years,
# so use only Full here to avoid counting an annual event twice.
annual_diagnostics <- annual %>% filter(horizon=="Full") %>%
  mutate(delta_P_Blim_pp=100*P_Blim_delta,
         delta_P_Bpa_pp=100*P_Bpa_delta,
         delta_P_closed_pp=100*P_closed_delta,
         REF_risk_gt_5=P_Blim_REF>0.05,
         robustness_risk_gt_5=P_Blim_rob>0.05)
stopifnot(nrow(annual_diagnostics)==360)
annual_summary <- annual_diagnostics %>%
  group_by(robustness_id, reference_scenario_number, reference_scenario_id, Fcap, Besc) %>%
  summarise(max_abs_delta_P_Blim_pp=max(abs(delta_P_Blim_pp)),
            max_abs_delta_P_closed_pp=max(abs(delta_P_closed_pp)),
            max_abs_delta_mean_catch=max(abs(mean_catch_delta)),
            year_largest_catch_change=year[which.max(abs(mean_catch_delta))],
            n_years_REF_risk_gt_5=sum(REF_risk_gt_5),
            n_years_robustness_risk_gt_5=sum(robustness_risk_gt_5), .groups="drop")

# Review whether ever-below-Blim status changes for a paired trajectory.
full_trajectory <- trajectory %>% filter(horizon=="Full")
stopifnot(all(c("min_SSB_REF", "min_SSB_rob") %in% names(full_trajectory)))
rob_eval <- readRDS(comparison$robustness_file)
thresholds <- rob_eval$performance_full %>% dplyr::select(reference_scenario_id, Blim)
stopifnot(nrow(thresholds)==12, !anyDuplicated(thresholds$reference_scenario_id))
risk_switches <- full_trajectory %>% left_join(thresholds, by="reference_scenario_id") %>%
  mutate(ever_below_REF=min_SSB_REF<Blim, ever_below_rob=min_SSB_rob<Blim) %>%
  filter(ever_below_REF!=ever_below_rob) %>%
  dplyr::select(robustness_id, reference_scenario_number, reference_scenario_id, Fcap, Besc,
                iter1000, om, replicate, min_SSB_REF, min_SSB_rob, ever_below_REF, ever_below_rob)

TAC_flags <- validation$TAC_flags
TAC_summary <- validation$TAC_summary
stopifnot(!is.null(TAC_flags), !is.null(TAC_summary))
if (nrow(TAC_flags)>0) {
  TAC_flags <- TAC_flags %>% mutate(shortfall_pct=100*(1-Catch_TAC_ratio))
}

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
tables <- list(trajectory_difference_summary=trajectory_summary,
               OM_difference_summary=om_summary, largest_trajectory_changes=extreme_trajectories,
               annual_diagnostics=annual_diagnostics, annual_difference_summary=annual_summary,
               ever_Blim_status_changes=risk_switches, TAC_flags=TAC_flags, TAC_summary=TAC_summary)
for (name in names(tables)) write.csv(tables[[name]], file.path(output_dir,paste0("31_",name,".csv")), row.names=FALSE)
saveRDS(list(robustness_id=robustness_id, source_comparison=comparison_file,
             source_validation=validation_file, tables=tables), file.path(output_dir,"31_robustness_diagnostics.rds"))
cat("\nPAIRED TRAJECTORY DIAGNOSTICS - FULL HORIZON\n")
print(trajectory_summary %>% filter(horizon=="Full", metric=="mean_catch") %>%
        dplyr::select(reference_scenario_number, Fcap, Besc, mean_delta, median_delta,
                      p05_delta, p95_delta, fraction_positive, fraction_negative), n=Inf, width=Inf)
cat("\nAnnual diagnostic summary:\n")
print(annual_summary, n=Inf, width=Inf)
cat("\nEver-below-Blim status changes:", nrow(risk_switches), "MP x trajectory cases\n")
cat("Catch/TAC flags:", nrow(TAC_flags), "\n")
cat("\nDIAGNOSTICS COMPLETED -", robustness_id, "\nOutputs saved in:", output_dir, "\n")

