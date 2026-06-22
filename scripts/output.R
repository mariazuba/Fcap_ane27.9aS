################################################################################
# output.R
# Extract MSE outputs and calculate performance statistics for Fcap x Besc grid
# Anchovy 9aS - capelin-style performance statistics
################################################################################

rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(FLBEIA)
  library(FLCore)
  library(FLFishery)
  library(dplyr)
  library(readr)
  library(here)
  library(stringr)
  library(tidyr)
  library(purrr)
})

setwd(here())

#===============================================================================
# 1. User options
#===============================================================================
# Definición de paŕametros de análisis
experiment_name <- "Fcap_Besc_grid"
stock_name      <- "ANE"
ssb_season      <- 2
rec_season      <- 3
risk_threshold  <- 0.05 # riesgo aceptable = 5%
open_threshold  <- 0.05 # sólo se consideran años "abiertos" si la pesquería abre en al menos 5% de iteraciones
closure_threshold <- 1e-8 # umbral numérico para considerar TAC/Catch igual a cero

# If TRUE, keeps only one file per recruitment_scenario x Fcap x Bescapement,
# prioritising the largest scenario_id. Useful when old and new runs are mixed.
keep_unique_grid <- TRUE # si es TRUE deja solo un archivo por combinación

#===============================================================================
# 2. Paths
#===============================================================================
# Se define donde estan los resultados y dóde guarda salidas
res_dir  <- here("outputs", "mse", "res")
meta_dir <- here("outputs", "mse", "summary")
out_dir  <- here("outputs", "Fcap_results", experiment_name)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_meta_file <- function(res_file) {
  file.path(meta_dir,paste0(tools::file_path_sans_ext(basename(res_file)), "_metadata.rds"))
}

#===============================================================================
# 3. Select files
#===============================================================================
 # Busca archivos .rds con nombre tipo: sc100_SR_bh_extreme_Fcap_2p00_Besc_10000.rds
# extrae desde el nombre: scenario_id, SR_model, recruitment_scenario, Fcap, Bescapement

mse_files_raw <- list.files(res_dir, 
                            pattern = "^sc[0-9]+_SR_bh_(hist|SS3|s05|s07|extreme)_Fcap_.*_Besc_.*_CP_cluster_[0-9]+\\.rds$",
                            full.names = TRUE) |> sort()

if (length(mse_files_raw) == 0) {stop("No MSE result files found in: ", res_dir)}

file_grid <- tibble(
  file_path            = mse_files_raw,
  file                 = basename(file_path),
  scenario_id          = str_extract(file, "^sc[0-9]+") |> str_remove("sc") |> as.integer(),
  SR_model             = str_extract(file, "(?<=_SR_)[A-Za-z0-9]+"),
  recruitment_scenario = str_extract(file, "(?<=_bh_)[A-Za-z0-9]+"),
  Fcap                 = str_extract(file, "(?<=Fcap_)[0-9]+p[0-9]+") |> str_replace("p", ".") |> as.numeric(),
  Bescapement          = str_extract(file, "(?<=Besc_)[0-9]+") |> as.numeric(),
  CP_cluster           = str_extract(file, "CP_cluster_[0-9]+"),
  meta_file            = map_chr(file_path, make_meta_file),
  has_metadata         = file.exists(meta_file) # verifica que exista el archivo de metadata asociado.
)

if (any(!file_grid$has_metadata)) {warning("Some files have no metadata. They will be removed from extraction.")}

file_grid <- file_grid |> filter(has_metadata)

# si es TRUE deja solo un archivo por combinación recruitment_scenarioxBescapement x Fcap 
# priorizando el scenario_id más alto
# (parámetro está definido arriba)
if (keep_unique_grid) {files_used <- file_grid |>
    arrange(recruitment_scenario, Bescapement, Fcap, CP_cluster, desc(scenario_id)) |>
    distinct(recruitment_scenario, Bescapement, Fcap, CP_cluster, .keep_all = TRUE)
} else { files_used <- file_grid}

mse_files <- files_used$file_path

message("Files found: ", length(mse_files_raw))
message("Files with metadata: ", nrow(file_grid))
message("Files used: ", length(mse_files))

write_csv(file_grid,  file.path(out_dir, "files_found.csv"))
write_csv(files_used, file.path(out_dir, "files_used.csv"))

#===============================================================================
# 4. Helper functions ----
#===============================================================================

#-------------------------------------------------------------------------------
## FUNCIÓN to_df_value() ----
#-------------------------------------------------------------------------------
# la Función to_df_value(): convierte objetos FLQuant a data.frame
# Por ejemplo, transforma SSB, reclutamiento, TAC o catch en tablas con 
# columnas year, iter, SSB/Recruitment/TAC/Catch
# Es una función auxiliar para evitar repetir código

to_df_value <- function(x, value_name) {df <- as.data.frame(x)
  
                  if ("data" %in% names(df)) {
                      names(df)[names(df) == "data"] <- value_name
                  } else {value_col <- setdiff(names(df), c("age", "year", "unit", "season", "area", "iter"))[1]
                      names(df)[names(df) == value_col] <- value_name
                          }
  
              df |> mutate(
              year = as.numeric(as.character(year)),
              iter = as.integer(as.character(iter)),
              !!value_name := as.numeric(.data[[value_name]]))
}
#-------------------------------------------------------------------------------
## FUNCIÓN calc_aav() ----
#-------------------------------------------------------------------------------
# la Función calc_aav(): Calcula el Average Annual Variation, es decir, 
# la variación relativa media entre años consecutivos: abs((valor_actual - valor_anterior)/valor_anterior)
# Se usa para TAC y Catch.  Sirve para evaluar estabilidad interanual de la regla

calc_aav <- function(df, value_col) {
              df |> arrange(iter, year) |> group_by(iter) |>
              mutate(prev_value = lag(.data[[value_col]]),
                    rel_change = ifelse(
              is.na(prev_value) | prev_value == 0,NA_real_,
              abs((.data[[value_col]] - prev_value) / prev_value))) |>
              summarise(AAV = mean(rel_change, na.rm = TRUE), .groups = "drop")
}
#-------------------------------------------------------------------------------
## FUNCIÓN extract_scenario_info() ----
#-------------------------------------------------------------------------------
# Extrae información del escenario desde: 1) el nombre del archivo y 2) metadata, si existe
# devuelve variables como: scenario_id, SR_model_recruitment_scenario, Fcap, Bescapement, Blim, Bpa
# esto permite que todas las tablas finales mantengan identificada cada combinación de HCR

extract_scenario_info <- function(res_file, meta = NULL) {file_name <- basename(res_file)
  
                          from_file <- tibble(scenario_file        = file_name,
                                              scenario_id_file     = str_extract(file_name, "^sc[0-9]+") |> str_remove("sc") |>as.integer(),
                                              SR_model_file        = str_extract(file_name, "(?<=_SR_)[A-Za-z0-9]+"),
                                              recruitment_scenario = str_extract(file_name, "(?<=_bh_)[A-Za-z0-9]+"),
                                              Fcap_file            = str_extract(file_name, "(?<=Fcap_)[0-9]+p[0-9]+") |> str_replace("p", ".") |> as.numeric(),
                                              Bescapement_file     = str_extract(file_name, "(?<=Besc_)[0-9]+") |> as.numeric(),
                                              CP_cluster_file      = str_extract(file_name, "(CP_cluster_)[0-9]+"))
  
  if (!is.null(meta) && "scenario" %in% names(meta)) {sc <- as_tibble(meta$scenario)
    
    if ("SR" %in% names(sc)) sc <- rename(sc, SR_model = SR)
    if ("obs" %in% names(sc)) sc <- rename(sc, obs_error = obs)
    if ("Besc" %in% names(sc)) sc <- rename(sc, Bescapement = Besc)
    
  } else {sc <- tibble()}
  
  if (!"scenario_id"          %in% names(sc)) sc$scenario_id          <- from_file$scenario_id_file
  if (!"SR_model"             %in% names(sc)) sc$SR_model             <- from_file$SR_model_file
  if (!"recruitment_scenario" %in% names(sc)) sc$recruitment_scenario <- from_file$recruitment_scenario
  if (!"CP_cluster"           %in% names(sc)) sc$CP_cluster           <- from_file$CP_cluster_file
  if (!"obs_error"            %in% names(sc)) sc$obs_error            <- NA_character_
  if (!"assess"               %in% names(sc)) sc$assess               <- NA_character_
  if (!"advice"               %in% names(sc)) sc$advice               <- NA_character_
  if (!"Fcap"                 %in% names(sc)) sc$Fcap                 <- from_file$Fcap_file
  if (!"Bescapement"          %in% names(sc)) sc$Bescapement          <- from_file$Bescapement_file
  if (!"catch_prop"           %in% names(sc)) sc$catch_prop           <- NA_character_
  if (!"ni"                   %in% names(sc)) sc$ni                   <- NA_integer_
  if (!"proj_nyr"             %in% names(sc)) sc$proj_nyr             <- NA_integer_
  
  sc |> mutate(
      scenario_file        = file_name,
      recruitment_scenario = from_file$recruitment_scenario,
      CP_cluster           = from_file$CP_cluster_file,
      catch_pattern        = catch_prop,
      implementation_error = NA_character_,
      Blim = if (!is.null(meta) && "Blim" %in% names(meta)) meta$Blim else NA_real_,
      Bpa  = if (!is.null(meta) && "Bpa"  %in% names(meta)) meta$Bpa  else NA_real_) |>
      relocate(any_of(c("scenario_file", "scenario_id", "CP_cluster", "SR_model", "recruitment_scenario",
        "obs_error", "assess", "advice", "catch_pattern", "implementation_error",
        "Fcap", "Bescapement", "Blim", "Bpa", "catch_prop", "ni", "proj_nyr")))
}

#===============================================================================
# 5. Extract one result file ----
#===============================================================================

#-------------------------------------------------------------------------------
## FUNCIÓN extract_mse_outputs() ----
#-------------------------------------------------------------------------------

extract_mse_outputs <- function(res_file) {meta_file <- make_meta_file(res_file)
                        if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)
  
                        # primero lee un archivo de resultados
                        res  <- readRDS(res_file)
                        meta <- readRDS(meta_file)
                        stk  <- res$stocks[[stock_name]]
                        # define años de proyección
                        proj_yrs      <- as.numeric(meta$proj.yrs)
                        proj_yrs_eval <- proj_yrs[-length(proj_yrs)] 
                        all_stk_yrs   <- as.numeric(dimnames(stk)$year)
                        hist_yrs      <- all_stk_yrs[all_stk_yrs < min(proj_yrs)]
                        plot_yrs      <- c(hist_yrs, proj_yrs_eval)
                        yrs_eval      <- as.character(proj_yrs_eval)
                        yrs_plot      <- as.character(plot_yrs)
                        sc_info       <- extract_scenario_info(res_file, meta)
                        #-------------------------------------------------------
                        # Historical SSB for plotting complete time series
                        #-------------------------------------------------------
                        # Extrae la SSB histórica, desde el primer año disponible hasta el año anterior
                        # al inicio de la proyección. Se guarda separada de las variables de desempeño,
                        # porque TAC, catch y recruitment se evalúan sólo en los años proyectados.
                        
                        yrs_hist <- as.character(hist_yrs)
                        
                        ssb_hist_df <- to_df_value(
                          ssb(stk)[, ac(yrs_hist), , ssb_season],
                          "SSB"
                        ) |>
                          select(year, iter, SSB) |>
                          bind_cols(sc_info) |>
                          group_by(across(all_of(names(sc_info))), year) |>
                          summarise(
                            median_SSB = median(SSB, na.rm = TRUE),
                            p05_SSB    = quantile(SSB, 0.05, na.rm = TRUE),
                            p95_SSB    = quantile(SSB, 0.95, na.rm = TRUE),
                            .groups = "drop"
                          ) |>
                          mutate(period = "Historical")
                        # Extrae por año de iteración
                        ssb_df   <- to_df_value(ssb(stk)[, ac(yrs_eval), , ssb_season], "SSB") |> select(year, iter, SSB)
                        rec_df   <- to_df_value(rec(stk)[, ac(yrs_eval), , rec_season], "Recruitment") |> select(year, iter, Recruitment)
                        tac_df   <- to_df_value(res$advice$TAC[stock_name, ac(yrs_eval)], "TAC") |> select(year, iter, TAC)
                        catch_df <- to_df_value(catch(stk)[, ac(yrs_eval)], "Catch") |> select(year, iter, Catch)
                        # Luego se une todo en una tabla anual y crea indicadores 
                        # below_Blim, below_Bescapement, fishery_closed, fishery_open
                        # es decir, para cada año e iteración dice sí: 
                        # - la SSB cae bajo Blim
                        # - la SSB cae bajo Bescapement
                        # - la pesquería queda cerrada
                        # - la pesquería queda abierta
                        annual_df <- ssb_df |> 
                                    left_join(rec_df,   by = c("year", "iter")) |>
                                    left_join(tac_df,   by = c("year", "iter")) |>
                                    left_join(catch_df, by = c("year", "iter")) |>
                                    bind_cols(sc_info) |>
                                    mutate(below_Blim        = SSB < Blim,
                                           below_Bescapement = SSB < Bescapement,
                                           tac_zero          = TAC <= closure_threshold,
                                           catch_zero        = Catch <= closure_threshold,
                                           fishery_closed    = tac_zero,
                                           fishery_open      = !fishery_closed)
                        #-------------------------------------------------------
                        #  Performance statistics:
                        #-------------------------------------------------------
                        # Primero calcula medianas por iteración a través de los años
                        iter_perf_all_years <- annual_df |> group_by(across(all_of(names(sc_info))), iter) |>
                                                summarise(median_SSB_iter         = median(SSB, na.rm = TRUE),
                                                          median_TAC_iter         = median(TAC, na.rm = TRUE),
                                                          median_Catch_iter       = median(Catch, na.rm = TRUE),
                                                          median_Recruitment_iter = median(Recruitment, na.rm = TRUE),
                                                          n_years_closed_iter     = sum(fishery_closed, na.rm = TRUE), .groups = "drop")
  
                        iter_perf_open_years <- annual_df |> filter(fishery_open) |> group_by(across(all_of(names(sc_info))), iter) |>
                                                summarise(median_SSB_open_iter   = median(SSB, na.rm = TRUE),
                                                          median_TAC_open_iter   = median(TAC, na.rm = TRUE),
                                                          median_Catch_open_iter = median(Catch, na.rm = TRUE), .groups = "drop")
  
                        tac_aav   <- calc_aav(annual_df, "TAC") |> rename(AAV_TAC = AAV)
                        catch_aav <- calc_aav(annual_df, "Catch") |>rename(AAV_Catch = AAV)
                        aav_df    <- tac_aav |>left_join(catch_aav, by = "iter")
                        # Luego resume entre iteraciones (_all_years) y 
                        # también calcula estadísticas solo para años donde la pesquería estuvo abierta (_open_years)
                        # Esto es útil porque algunas reglas pueden cerrar la pesquería muchos años, y 
                        # eso cambia la interpretación del catch medio.
                        summary_df <- iter_perf_all_years |> left_join(iter_perf_open_years, by = c(names(sc_info), "iter")) |>
                                      group_by(across(all_of(names(sc_info)))) |>
                                      summarise(median_SSB_all_years = median(median_SSB_iter, na.rm = TRUE),
                                                p05_SSB_all_years    = quantile(median_SSB_iter, 0.05, na.rm = TRUE),
                                                p95_SSB_all_years    = quantile(median_SSB_iter, 0.95, na.rm = TRUE),
      
                                                median_TAC_all_years = median(median_TAC_iter, na.rm = TRUE),
                                                p05_TAC_all_years    = quantile(median_TAC_iter, 0.05, na.rm = TRUE),
                                                p95_TAC_all_years    = quantile(median_TAC_iter, 0.95, na.rm = TRUE),
      
                                                median_Catch_all_years = median(median_Catch_iter, na.rm = TRUE),
                                                p05_Catch_all_years    = quantile(median_Catch_iter, 0.05, na.rm = TRUE),
                                                p95_Catch_all_years    = quantile(median_Catch_iter, 0.95, na.rm = TRUE),
      
                                                median_SSB_open_years  = median(median_SSB_open_iter, na.rm = TRUE),
                                                p05_SSB_open_years     = quantile(median_SSB_open_iter, 0.05, na.rm = TRUE),
                                                p95_SSB_open_years     = quantile(median_SSB_open_iter, 0.95, na.rm = TRUE),
      
                                                median_TAC_open_years  = median(median_TAC_open_iter, na.rm = TRUE),
                                                p05_TAC_open_years     = quantile(median_TAC_open_iter, 0.05, na.rm = TRUE),
                                                p95_TAC_open_years     = quantile(median_TAC_open_iter, 0.95, na.rm = TRUE),
      
                                                median_Catch_open_years = median(median_Catch_open_iter, na.rm = TRUE),
                                                p05_Catch_open_years    = quantile(median_Catch_open_iter, 0.05, na.rm = TRUE),
                                                p95_Catch_open_years    = quantile(median_Catch_open_iter, 0.95, na.rm = TRUE),
      
                                                median_Recruitment_all_years = median(median_Recruitment_iter, na.rm = TRUE),
      
                                                n_years_closed_mean   = mean(n_years_closed_iter, na.rm = TRUE),
                                                n_years_closed_median = median(n_years_closed_iter, na.rm = TRUE),
                                                n_years_closed_p05    = quantile(n_years_closed_iter, 0.05, na.rm = TRUE),
                                                n_years_closed_p95    = quantile(n_years_closed_iter, 0.95, na.rm = TRUE), .groups = "drop") |>
                                        bind_cols(aav_df |>
                                        summarise(mean_AAV_TAC = mean(AAV_TAC, na.rm = TRUE),
                                                  mean_AAV_Catch = mean(AAV_Catch, na.rm = TRUE)))
  
                        list(
                          annual   = annual_df,
                          summary  = summary_df,
                          ssb_hist = ssb_hist_df
                        )
}

#################################################################################################################################################
#################################################################################################################################################

#===============================================================================
# 6. Extract all files ----
#===============================================================================
# aplica la función a todos los escenarios
all_outputs  <- lapply(mse_files, extract_mse_outputs)
# después une todo:
# annual_all queda como la gran tabla base con todos los escenarios, años e iteraciones
annual_all   <- bind_rows(lapply(all_outputs, `[[`, "annual"))
summary_all  <- bind_rows(lapply(all_outputs, `[[`, "summary"))
ssb_hist_all <- bind_rows(lapply(all_outputs, `[[`, "ssb_hist"))
#===============================================================================
# HORIZONTES TEMPORALES ----
#===============================================================================
# divide la proyección en: short, medium, long, según el tercio del período proyectado
# esto permite comparar desempeño a corto, medio y largo plazo
proj_years   <- sort(unique(annual_all$year))
n_proj       <- length(proj_years)

annual_all <- annual_all |> mutate(year_index = match(year, proj_years),
                                    horizon   = case_when(year_index <= ceiling(n_proj / 3) ~ "Short",
                                                 year_index <= ceiling(2 * n_proj / 3) ~ "Medium", TRUE ~ "Long"),
                                    horizon = factor(horizon, levels = c("Short", "Medium", "Long")))


scenario_cols <- c(
  "scenario_file", "scenario_id","CP_cluster", "SR_model", "recruitment_scenario",
  "catch_pattern", "obs_error", "implementation_error",
  "Fcap", "Bescapement", "Blim", "Bpa"
)
scenario_cols <- intersect(scenario_cols, names(annual_all))



#===============================================================================
# 7. Annual probabilities ----
#===============================================================================

#-------------------------------------------------------
# Riesgo anual ----
#-------------------------------------------------------
# Calcula, por año y escenario y marca si es riesgo anual cumple con el umbral del 5%
risk_by_year <- annual_all |> group_by(across(all_of(scenario_cols)), year) |>
                summarise(n_iter                     = n(),
                          prob_fishery_open          = mean(fishery_open, na.rm = TRUE),
                          prob_fishery_closed        = mean(fishery_closed, na.rm = TRUE),
                          prob_SSB_below_Blim        = mean(SSB < Blim, na.rm = TRUE),
                          prob_SSB_below_Bescapement = mean(SSB < Bescapement, na.rm = TRUE), 
                          .groups = "drop") |>
                mutate(risk_ok_annual_Blim           = prob_SSB_below_Blim <= risk_threshold,
                       risk_ok_annual_Bescapement    = prob_SSB_below_Bescapement <= risk_threshold)
#-------------------------------------------------------
# Riesgo tipo ICES categoría 3
#-------------------------------------------------------
# Years where the fishery opened in at least 5% of iterations under any HCR,
# calculated separately for each recruitment scenario.

# Primero identifia los años donde alguna HCR abre la pesquería al menos 5%
open_years_any_hcr <- risk_by_year |> group_by(recruitment_scenario, year) |>
                      summarise(prob_open_any_hcr = max(prob_fishery_open, na.rm = TRUE), .groups = "drop") |>
                      filter(prob_open_any_hcr >= open_threshold)

# Luego calcula el máximo riesgo en esos años
# la regla candidata cumple si max_risk3_Blim <=0.05. Es decir, si el riesgo máximo de caer bajo Blim no supera 5%
risk3_by_hcr      <- risk_by_year |> semi_join(open_years_any_hcr, by = c("recruitment_scenario", "year")) |>
                     group_by(across(all_of(scenario_cols))) |>
                     summarise(max_risk3_Blim                 = max(prob_SSB_below_Blim, na.rm = TRUE),
                               n_years_risk3_Blim_gt_5        = sum(prob_SSB_below_Blim > risk_threshold, na.rm = TRUE),
                               max_risk3_Bescapement          = max(prob_SSB_below_Bescapement, na.rm = TRUE),
                               n_years_risk3_Bescapement_gt_5 = sum(prob_SSB_below_Bescapement > risk_threshold, na.rm = TRUE),
                               risk_ok_ices_type3             = max_risk3_Blim <= risk_threshold, .groups = "drop")

# Also keep risks across all years for diagnostics.
risk_all_years_by_hcr <- risk_by_year |> group_by(across(all_of(scenario_cols))) |>
                         summarise(max_risk_all_years_Blim        = max(prob_SSB_below_Blim, na.rm = TRUE),
                                   n_years_all_Blim_gt_5          = sum(prob_SSB_below_Blim > risk_threshold, na.rm = TRUE),
                                   max_risk_all_years_Bescapement = max(prob_SSB_below_Bescapement, na.rm = TRUE),
                                   n_years_all_Bescapement_gt_5   = sum(prob_SSB_below_Bescapement > risk_threshold, na.rm = TRUE), .groups = "drop")

#===============================================================================
# 8. Annual summaries for time series plots ----
#===============================================================================

annual_summary <- annual_all |> group_by(across(all_of(scenario_cols)), year) |>
                  summarise(median_SSB = median(SSB, na.rm = TRUE),
                            p05_SSB    = quantile(SSB, 0.05, na.rm = TRUE),
                            p25_SSB    = quantile(SSB, 0.25, na.rm = TRUE),
                            p75_SSB    = quantile(SSB, 0.75, na.rm = TRUE),
                            p95_SSB    = quantile(SSB, 0.95, na.rm = TRUE),
    
                            median_Recruitment = median(Recruitment, na.rm = TRUE),
                            p05_Recruitment    = quantile(Recruitment, 0.05, na.rm = TRUE),
                            p95_Recruitment    = quantile(Recruitment, 0.95, na.rm = TRUE),
    
                            median_TAC = median(TAC, na.rm = TRUE),
                            p05_TAC = quantile(TAC, 0.05, na.rm = TRUE),
                            p95_TAC = quantile(TAC, 0.95, na.rm = TRUE),
    
                            median_Catch = median(Catch, na.rm = TRUE),
                            p05_Catch = quantile(Catch, 0.05, na.rm = TRUE),
                            p95_Catch = quantile(Catch, 0.95, na.rm = TRUE), .groups = "drop")

#===============================================================================
# 9. Final performance table ----
#===============================================================================

# Une estadísticas biológicas, pesqueras y de riesgo
performance_table <- summary_all |> left_join(risk3_by_hcr |>
                     select(all_of(scenario_cols),
                            max_risk3_Blim,
                            n_years_risk3_Blim_gt_5,
                            max_risk3_Bescapement,
                            n_years_risk3_Bescapement_gt_5,
                            risk_ok_ices_type3), by = scenario_cols) |>
                    left_join(risk_all_years_by_hcr, by = scenario_cols) |>
                    arrange(recruitment_scenario, Bescapement, Fcap)
# Luego identifica reglas candidatas
# y las ordena priorizando: 1. que cumple riesgo, 2. mayor catch en años abiertos, 3. menor variabilidad del TAC
candidate_hcrs <- performance_table |>filter(risk_ok_ices_type3) |>
                  arrange(recruitment_scenario, desc(median_Catch_open_years), mean_AAV_TAC)

# Distancia a Blim ----
# Calcula qué tan cerca queda cada escenario del límite biológico. 
# Esto sirve para ver cuáles reglas están más cerca de una situación peligrosa
distance_to_Blim <- annual_all |> group_by(across(all_of(scenario_cols))) |>
                    summarise(min_SSB           = min(SSB, na.rm = TRUE),
                              min_SSB_over_Blim = min_SSB / unique(Blim),
                              safety_margin     = min_SSB - unique(Blim), .groups = "drop") |>
                    arrange(min_SSB_over_Blim)

# Performance horizon ----
performance_horizons <- annual_all |> group_by(across(all_of(scenario_cols)), horizon, iter) |>
                        summarise(median_SSB_all_years   = median(SSB, na.rm = TRUE),
                                  median_Catch_all_years = median(Catch, na.rm = TRUE),
                                  median_TAC_all_years   = median(TAC, na.rm = TRUE),
                                  n_years_closed         = sum(fishery_closed, na.rm = TRUE), .groups = "drop") |>
                        group_by(across(all_of(scenario_cols)), horizon) |>
                        summarise(median_SSB        = median(median_SSB_all_years, na.rm = TRUE),
                                  median_Catch      = median(median_Catch_all_years, na.rm = TRUE),
                                  median_TAC        = median(median_TAC_all_years, na.rm = TRUE),
                                  mean_years_closed = mean(n_years_closed, na.rm = TRUE),.groups = "drop")

write_csv(performance_horizons, file.path(out_dir, "performance_horizons.csv"))
#===============================================================================
# 10. Save outputs
#===============================================================================

report_outputs <-  list(annual_all            = annual_all,
                        annual_summary        = annual_summary,
                        ssb_hist_all          = ssb_hist_all,
                        risk_by_year          = risk_by_year,
                        open_years_any_hcr    = open_years_any_hcr,
                        risk3_by_hcr          = risk3_by_hcr,
                        risk_all_years_by_hcr = risk_all_years_by_hcr,
                        performance_table     = performance_table,
                        candidate_hcrs        = candidate_hcrs,
                        distance_to_Blim      = distance_to_Blim,
                        scenario_cols         = scenario_cols,
                        files_found           = file_grid,
                        files_used            = files_used,
                        settings              = list(experiment_name = experiment_name,
                                                     stock_name      = stock_name,
                                                     ssb_season      = ssb_season,
                                                     rec_season      = rec_season,
                                                     risk_threshold  = risk_threshold,
                                                     open_threshold  = open_threshold))

saveRDS(report_outputs, file.path(out_dir, "report_outputs.rds"))
saveRDS(report_outputs, file.path(out_dir, paste0(experiment_name, "_report_outputs.rds")))

write_csv(performance_table, file.path(out_dir, "performance_table.csv"))
write_csv(annual_summary, file.path(out_dir, "annual_summary.csv"))
write_csv(risk_by_year, file.path(out_dir, "risk_by_year.csv"))
write_csv(risk3_by_hcr, file.path(out_dir, "risk3_by_hcr.csv"))
write_csv(candidate_hcrs, file.path(out_dir, "candidate_hcrs.csv"))
write_csv(distance_to_Blim, file.path(out_dir, "distance_to_Blim.csv"))

message("Done. Outputs saved in: ", out_dir)

#==============================================================================
# Diagnostic: does the HCR close the fishery?
#==============================================================================

diagnose_closure <- function(res_file, meta_file = NULL) {
  
  res <- readRDS(res_file)
  
  if (is.null(meta_file)) {
    meta_file <- gsub("outputs/mse/res", "outputs/mse/summary", res_file)
    meta_file <- gsub("\\.rds$", "_metadata.rds", meta_file)
  }
  
  meta <- readRDS(meta_file)
  stk  <- res$stocks$ANE
  
  yrs <- meta$proj.yrs
  yrs_eval <- yrs[-length(yrs)]
  
  tac_df <- as.data.frame(res$advice$TAC["ANE", ac(yrs_eval)])
  names(tac_df)[names(tac_df) == "data"] <- "TAC"
  
  catch_df <- as.data.frame(catch(stk)[, ac(yrs_eval)])
  names(catch_df)[names(catch_df) == "data"] <- "Catch"
  
  ssb_df <- as.data.frame(ssb(stk)[, ac(yrs_eval), , 2])
  names(ssb_df)[names(ssb_df) == "data"] <- "SSB"
  
  out <- ssb_df |> 
    dplyr::select(year, iter, SSB) |>
    dplyr::left_join(tac_df   |> dplyr::select(year, iter, TAC),   by = c("year", "iter")) |>
    dplyr::left_join(catch_df |> dplyr::select(year, iter, Catch), by = c("year", "iter")) |> 
    dplyr::mutate(
      year              = as.numeric(as.character(year)),
      iter              = as.integer(as.character(iter)),
      SSB               = as.numeric(SSB),
      TAC               = as.numeric(TAC),
      Catch             = as.numeric(Catch),
      Bescapement       = meta$Bescapement,
      Blim              = meta$Blim,
      TAC_zero          = TAC <= closure_threshold,
      Catch_zero        = Catch <= closure_threshold,
      below_Bescapement = SSB < Bescapement,
      below_Blim        = SSB < Blim)
  
  summary_total <- out |> dplyr::summarise(
      n                      = dplyr::n(),
      min_TAC                = min(TAC, na.rm = TRUE),
      p05_TAC                = quantile(TAC, 0.05, na.rm = TRUE),
      median_TAC             = median(TAC, na.rm = TRUE),
      max_TAC                = max(TAC, na.rm = TRUE),
      n_TAC_zero             = sum(TAC_zero, na.rm = TRUE),
      prop_TAC_zero          = mean(TAC_zero, na.rm = TRUE),
      min_Catch              = min(Catch, na.rm = TRUE),
      n_Catch_zero           = sum(Catch_zero, na.rm = TRUE),
      prop_Catch_zero        = mean(Catch_zero, na.rm = TRUE),
      prop_below_Bescapement = mean(below_Bescapement, na.rm = TRUE),
      prop_below_Blim        = mean(below_Blim, na.rm = TRUE))
  
  summary_year <- out |> dplyr::group_by(year) |> dplyr::summarise(
      n_iter                 = dplyr::n(),
      min_TAC                = min(TAC, na.rm = TRUE),
      prop_TAC_zero          = mean(TAC_zero, na.rm = TRUE),
      prop_Catch_zero        = mean(Catch_zero, na.rm = TRUE),
      prop_below_Bescapement = mean(below_Bescapement, na.rm = TRUE),
      prop_below_Blim        = mean(below_Blim, na.rm = TRUE), .groups = "drop")
  
  list(
    res_file      = basename(res_file),
    meta          = meta,
    raw           = out,
    summary_total = summary_total,
    summary_year  = summary_year)
}

# 
# res_file <- "outputs/mse/res/sc100_SR_bh_extreme_Fcap_2p00_Besc_10000.rds"
# 
# diag <- diagnose_closure(res_file)
# 
# diag$summary_total
# diag$summary_year
# 

