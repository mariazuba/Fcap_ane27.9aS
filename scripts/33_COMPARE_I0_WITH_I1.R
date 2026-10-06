# ============================================================
# 33_COMPARE_I0_WITH_I1.R
# Comparacion directa de informacion: diferencias I0 menos I1.
# Ambos tienen forecast determinista; I0 retira el assessment error.
# Lee las evaluaciones del script 29. No vuelve a simular.
# Empareja MP, horizonte, iter1000, OM y replica.
# Probabilidades en fracciones; 0.01 equivale a 1 punto porcentual.
# Las diferencias de medianas agregadas no son medianas de diferencias.
# Ejecutar desde la raiz: source("scripts/33_COMPARE_I0_WITH_I1.R")
# Reejecutar sobrescribe solamente las salidas de esta comparacion.
# ============================================================

library(dplyr)
library(here)

file_I0 <- here("data", "mse", "robustness_performance", "I0", "robustness_MP_performance.rds")
file_I1 <- here("data", "mse", "robustness_performance", "I1", "robustness_MP_performance.rds")
stopifnot(file.exists(file_I0), file.exists(file_I1))
I0 <- readRDS(file_I0)
I1 <- readRDS(file_I1)
stopifnot(identical(I0$robustness_id, "I0"), identical(I1$robustness_id, "I1"))
stopifnot(identical(I0$horizons, I1$horizons))
stopifnot(identical(I0$risk_threshold, I1$risk_threshold))

MP_keys <- c("reference_scenario_number", "reference_scenario_id", "Fcap", "Besc")
MP_settings <- I0$scenario_grid %>% dplyr::select(all_of(MP_keys))
stopifnot(nrow(MP_settings)==12, !anyDuplicated(MP_settings$reference_scenario_id))
stopifnot(nrow(anti_join(MP_settings, I1$scenario_grid, by=MP_keys))==0)
stopifnot(nrow(I1$scenario_grid)==12)

# Comprobar identicos umbrales y limites temporales por MP y horizonte.
settings <- c("Fcap", "Besc", "Blim", "Bpa", "first_year", "last_year")
setting_keys <- c("reference_scenario_id", "horizon", settings)
stopifnot(nrow(I0$performance)==48, nrow(I1$performance)==48)
stopifnot(nrow(anti_join(I0$performance, I1$performance, by=setting_keys))==0)

compare_I0_I1 <- function(table_I0, table_I1, extra_keys, expected_rows) {
  join_keys <- c("reference_scenario_id", extra_keys)
  stopifnot(all(join_keys %in% names(table_I0)), all(join_keys %in% names(table_I1)))
  stopifnot(nrow(table_I0)==expected_rows, nrow(table_I1)==expected_rows)
  stopifnot(!anyDuplicated(table_I0[join_keys]), !anyDuplicated(table_I1[join_keys]))
  stopifnot(nrow(anti_join(table_I0, table_I1, by=join_keys))==0)
  stopifnot(nrow(anti_join(table_I1, table_I0, by=join_keys))==0)
  excluded <- c(join_keys, "scenario_number", "reference_scenario_number", "Fcap", "Besc", "Blim", "Bpa", "first_year", "last_year", "n")
  common <- setdiff(intersect(names(table_I0), names(table_I1)), excluded)
  numeric_metrics <- common[vapply(common, function(x) is.numeric(table_I0[[x]]) && is.numeric(table_I1[[x]]), logical(1))]
  logical_metrics <- common[vapply(common, function(x) is.logical(table_I0[[x]]) && is.logical(table_I1[[x]]), logical(1))]
  metrics <- c(numeric_metrics, logical_metrics)
  stopifnot(length(numeric_metrics)>0)
  left <- table_I0 %>% dplyr::select(all_of(join_keys), all_of(metrics))
  right <- table_I1 %>% dplyr::select(all_of(join_keys), all_of(metrics))
  out <- left_join(left, right, by=join_keys, suffix=c("_I0", "_I1"))
  for (metric in numeric_metrics) {
    out[[paste0(metric, "_delta")]] <- out[[paste0(metric, "_I0")]] - out[[paste0(metric, "_I1")]]
  }
  out <- left_join(out, MP_settings, by="reference_scenario_id")
  stopifnot(nrow(out)==expected_rows)
  out
}

performance_I0_I1 <- compare_I0_I1(I0$performance, I1$performance, "horizon", 48)
annual_I0_I1 <- compare_I0_I1(I0$annual_metrics, I1$annual_metrics, c("horizon", "year"), 720)
trajectory_I0_I1 <- compare_I0_I1(I0$trajectory_metrics, I1$trajectory_metrics, c("horizon", "iter1000", "om", "replicate"), 48000)
full_I0_I1 <- performance_I0_I1 %>% filter(horizon=="Full") %>% arrange(Fcap, Besc)
trajectory_full_I0_I1 <- trajectory_I0_I1 %>% filter(horizon=="Full")
stopifnot(nrow(trajectory_full_I0_I1)==12000)
stopifnot(!anyNA(trajectory_full_I0_I1$ever_below_Blim_I0), !anyNA(trajectory_full_I0_I1$ever_below_Blim_I1))

# Tabla completa, incluidos casos que conservan su estado.
risk_status_I0_I1 <- trajectory_full_I0_I1 %>%
  count(ever_below_Blim_I1, ever_below_Blim_I0, name="n_cases")
risk_switches_I0_I1 <- trajectory_full_I0_I1 %>%
  filter(ever_below_Blim_I0 != ever_below_Blim_I1)
paired_catch_I0_I1 <- trajectory_full_I0_I1 %>%
  group_by(reference_scenario_number, Fcap, Besc) %>%
  summarise(mean_delta=mean(mean_catch_delta),
            median_delta=median(mean_catch_delta),
            p05_delta=as.numeric(quantile(mean_catch_delta, 0.05)),
            p95_delta=as.numeric(quantile(mean_catch_delta, 0.95)),
            fraction_positive=mean(mean_catch_delta>0),
            fraction_negative=mean(mean_catch_delta<0),
            fraction_zero=mean(mean_catch_delta==0),
            .groups="drop")

summary_I0_I1 <- full_I0_I1 %>% dplyr::select(
  reference_scenario_number, Fcap, Besc,
  max_P_Blim_I0, max_P_Blim_I1, max_P_Blim_delta,
  mean_catch_I0, mean_catch_I1, mean_catch_delta,
  mean_P_closed_I0, mean_P_closed_I1, mean_P_closed_delta,
  median_IAV_I0, median_IAV_I1, median_IAV_delta)

output_I0_I1 <- here("outputs", "mse", "robustness", "comparison", "I0_vs_I1")
dir.create(output_I0_I1, recursive=TRUE, showWarnings=FALSE)
write.csv(performance_I0_I1, file.path(output_I0_I1, "33_performance_comparison.csv"), row.names=FALSE)
write.csv(full_I0_I1, file.path(output_I0_I1, "33_full_horizon_comparison.csv"), row.names=FALSE)
write.csv(annual_I0_I1, file.path(output_I0_I1, "33_annual_comparison.csv"), row.names=FALSE)
write.csv(trajectory_I0_I1, file.path(output_I0_I1, "33_trajectory_comparison.csv"), row.names=FALSE)
write.csv(risk_status_I0_I1, file.path(output_I0_I1, "33_risk_status_counts.csv"), row.names=FALSE)
write.csv(risk_switches_I0_I1, file.path(output_I0_I1, "33_risk_switches_full.csv"), row.names=FALSE)
write.csv(paired_catch_I0_I1, file.path(output_I0_I1, "33_paired_catch_summary_full.csv"), row.names=FALSE)
saveRDS(list(difference_definition="I0 minus I1", source_I0=file_I0, source_I1=file_I1,
             performance=performance_I0_I1, performance_full=full_I0_I1,
             annual=annual_I0_I1, trajectory=trajectory_I0_I1,
             risk_status=risk_status_I0_I1, risk_switches=risk_switches_I0_I1,
             paired_catch=paired_catch_I0_I1),
        file.path(output_I0_I1, "33_comparison_I0_with_I1.rds"))
print(summary_I0_I1, n=Inf, width=Inf)
print(paired_catch_I0_I1, n=Inf, width=Inf)
print(risk_status_I0_I1)
cat("\nCOMPARACION COMPLETADA: I0 menos I1\n")
cat("MPs: 12 | Horizontes: 4 | Filas emparejadas: 48000\n")
cat("Resultados guardados en:", output_I0_I1, "\n")
