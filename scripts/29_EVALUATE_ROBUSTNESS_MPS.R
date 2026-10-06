# ============================================================
# 29_EVALUATE_ROBUSTNESS_MPS.R
# Performance evaluation of robustness MPs for anchovy ane.27.9aS
#
# Adapted from 22b_EVALUATE_REFINED_CANDIDATE_MPS.R.
# All original metric formulas and horizon definitions are retained.
# Reads 12 selected MPs x 10 blocks x 100 trajectories; no MSE reruns.
# Requires a successful script-28 validation for the selected scenario.
#
# Variables: recruitment Q3, SSB Q2, realised annual landings, TAC, Fadv.
# Horizons: Short 2025-2034; Medium 2035-2044; Long 2045-2054;
#           Full 2025-2054.
# Metrics: annual/ever biomass risk, closures (Fadv<=1e-8), realised
# catch, conditional open-year catch/SSB, and interannual variability.
# IAV uses |C(t)-C(t-1)|/[(C(t)+C(t-1))/2]; consecutive zero catches
# contribute zero. The first year of each horizon has no IAV step.
# Open-year metrics retain the original conditional definitions.
# Catch-TAC summaries are implementation diagnostics.
# No comparison, ranking or figures are generated here.
#
# R console, from project root:
#   source("scripts/29_EVALUATE_ROBUSTNESS_MPS.R")  # default: R-OM
#   Sys.setenv(MSE_EVALUATE_ROBUSTNESS_ID="S-C1")
#   source("scripts/29_EVALUATE_ROBUSTNESS_MPS.R")
#
# Inputs: script-26 design, script-28 validation, robustness block RDS,
#         block maps and metadata.
# Outputs: data/mse/robustness_performance/<robustness_id>/
#          robustness_MP_performance.rds and six CSV metric tables.
# REF performance files remain unchanged; reference keys are retained
# for comparison in script 30. Re-running replaces evaluation outputs.
# ============================================================

rm(list = ls())

library(FLCore)
library(FLBEIA)
library(dplyr)
library(tidyr)
library(purrr)
library(here)

# ============================================================
# DIRECTORIES AND SETTINGS
# ============================================================
robustness_id <- Sys.getenv("MSE_EVALUATE_ROBUSTNESS_ID", "R-OM")

stopifnot(robustness_id %in% c("R-OM", "S-C1", "S-C2", "I0", "I1"))

scenario_file <- here("outputs", "mse", "robustness", "design", "Table26_02_MP_robustness_design.csv")
scenario_root <- here("data", "mse", "robustness_runs")
performance_dir <- here("data", "mse", "robustness_performance", robustness_id)
validation_file <- here("outputs", "mse", "robustness", "validation", "post_run", robustness_id, "21c_refined_new_validation_summary.rds")

stopifnot(file.exists(scenario_file), file.exists(validation_file))

design <- read.csv(scenario_file, stringsAsFactors=FALSE)

required <- c("robustness_id", "reference_scenario_number", "reference_scenario_id", "robustness_run_id",
              "Fcap", "Besc", "Blim", "Bpa", "projection_start", "projection_end", "annual_risk_threshold")

stopifnot(all(required %in% names(design)))

scenario_grid <- design[design$robustness_id==robustness_id, , drop=FALSE]

stopifnot(nrow(scenario_grid)==12, !anyDuplicated(scenario_grid$reference_scenario_id))

stopifnot(setequal(scenario_grid$reference_scenario_number, c(9,10,11,12,17,18,19,20,25,26,27,28)))

scenario_grid$scenario_number <- scenario_grid$reference_scenario_number
scenario_grid$scenario_id <- scenario_grid$robustness_run_id
scenario_grid$first_projection_year <- scenario_grid$projection_start
scenario_grid$last_projection_year <- scenario_grid$projection_end
scenario_grid <- scenario_grid %>% arrange(scenario_number)

stopifnot(all(scenario_grid$projection_start==2025), all(scenario_grid$projection_end==2054))
stopifnot(all(scenario_grid$annual_risk_threshold==0.05))

post_run_validation <- readRDS(validation_file)

stopifnot(identical(post_run_validation$settings$robustness_id, robustness_id))
stopifnot(nrow(post_run_validation$scenario_validation)==12, nrow(post_run_validation$block_validation)==120)
stopifnot(all(post_run_validation$scenario_validation$structural_PASS))
stopifnot(all(post_run_validation$block_validation$structural_PASS))
stopifnot(setequal(post_run_validation$scenario_validation$scenario_id, scenario_grid$scenario_id))

dir.create(performance_dir, recursive=TRUE, showWarnings=FALSE)
block_ids <- 1:10
risk_threshold <- 0.05
horizons <- list(Short=2025:2034, Medium=2035:2044, Long=2045:2054, Full=2025:2054)

print(scenario_grid %>% dplyr::select(scenario_number, scenario_id, Fcap, Besc))

cat("\nRobustness configuration:", robustness_id, "| Selected MPs:", nrow(scenario_grid), "\n")

# ============================================================
# HELPERS
# ============================================================

add_global_map <- function(dat,map_b) {
  out <- dat %>% mutate(iter_local=as.integer(as.character(iter))) %>% 
    left_join(map_b %>% dplyr::select(iter_local,iter1000,om,replicate),by="iter_local")
  stopifnot(!anyNA(out$iter1000),!anyNA(out$om),!anyNA(out$replicate))
  out
}

standardise_dynamic <- function(x,map_b,proj_years) {
  as.data.frame(x) %>% rename(value=data) %>% 
    mutate(year=as.integer(as.character(year))) %>% 
    add_global_map(map_b) %>% filter(year %in% proj_years) %>% 
    transmute(year,iter1000,om,replicate,value)
}

validate_MP_output <- function(x,variable_name) {
  x %>% summarise(variable=variable_name,
                  n=dplyr::n(),
                  n_OM=n_distinct(om),
                  n_trajectories=n_distinct(iter1000),
                  first_year=min(year),
                  last_year=max(year),n_NA=sum(is.na(value)),
                  n_nonfinite=sum(!is.finite(value)),
                  n_negative=sum(value<0,na.rm=TRUE),
                  min=min(value,na.rm=TRUE),
                  max=max(value,na.rm=TRUE))
}

# ============================================================
# EVALUATE ONE CANDIDATE MP
# ============================================================

evaluate_candidate_MP <- function(scenario_row) {
  
  scenario_number <- scenario_row$scenario_number[[1]]
  scenario_id <- scenario_row$scenario_id[[1]]
  Fcap <- scenario_row$Fcap[[1]]
  Besc <- scenario_row$Besc[[1]]
  Blim <- scenario_row$Blim[[1]]
  Bpa <- scenario_row$Bpa[[1]]
  first_proj <- scenario_row$first_projection_year[[1]]
  last_proj <- scenario_row$last_projection_year[[1]]
  proj_years <- first_proj:last_proj
  
  cat("\n============================================================\n")
  cat("SCENARIO:",scenario_number,"\n")
  cat("ID:",scenario_id,"\n")
  cat("Fcap:",Fcap,"| Besc:",Besc,"\n")
  cat("============================================================\n")
  
  # ==========================================================
  # FILES
  # ==========================================================
  
  mp_dir <- file.path(scenario_root,scenario_id)
  
  mp_files <- file.path(mp_dir,sprintf("block_%02d.rds",block_ids))
  map_files <- file.path(mp_dir,sprintf("block_%02d_map.csv",block_ids))
  metadata_files <- file.path(mp_dir,sprintf("block_%02d_metadata.csv",block_ids))
  
  if(!dir.exists(mp_dir)) stop("Scenario directory not found: ",mp_dir)
  if(!all(file.exists(mp_files))) stop("Missing RDS blocks for scenario ",scenario_number)
  if(!all(file.exists(map_files))) stop("Missing map files for scenario ",scenario_number)
  if(!all(file.exists(metadata_files))) stop("Missing metadata files for scenario ",scenario_number)
  
  # ==========================================================
  # CHECK 1 - COMPLETE MP ENSEMBLE
  # ==========================================================
  
  iteration_maps <- map_dfr(block_ids,~read.csv(map_files[.x]) %>% mutate(block=.x))
  metadata_all <- map_dfr(block_ids,~read.csv(metadata_files[.x]) %>% mutate(block=.x))
  
  stopifnot(nrow(iteration_maps)==1000)
  stopifnot(n_distinct(iteration_maps$iter1000)==1000)
  stopifnot(n_distinct(iteration_maps$om)==100)
  stopifnot(all(table(iteration_maps$om)==10))
  stopifnot(all(sort(unique(iteration_maps$replicate))==1:10))
  stopifnot(all(metadata_all$Fcap==Fcap))
  stopifnot(all(metadata_all$Besc==Besc))
  stopifnot(all(metadata_all$Blim==Blim))
  stopifnot(all(metadata_all$first_projection_year==first_proj))
  stopifnot(all(metadata_all$last_projection_year==last_proj))
  
  cat("CHECK 1 PASSED - Complete ensemble: 100 OMs x 10 replicates = 1000 trajectories\n")
  
  # ==========================================================
  # EXTRACT ONE BLOCK
  # ==========================================================
  
  extract_MP_block <- function(block_id) {
    
    OM_b <- readRDS(mp_files[block_id])
    map_b <- read.csv(map_files[block_id])
    
    biol <- OM_b$biols$ANE
    cobj <- OM_b$fleets$SEINE@metiers$ALL@catches$ANE
    
    R_b <- standardise_dynamic(biol@n["0",,,"3",,],map_b,proj_years)
    Mat <- predict(biol@mat)
    SSB_q2 <- quantSums(biol@n[,,,"2",,]*biol@wt[,,,"2",,]*Mat[,,,"2",,])
    B_b <- standardise_dynamic(SSB_q2,map_b,proj_years)
    C_ann <- seasonSums(quantSums(cobj@landings))
    C_b <- standardise_dynamic(C_ann,map_b,proj_years)
    TAC_b <- standardise_dynamic(OM_b$advice$TAC["ANE",,,,,],map_b,proj_years)
    Fadv_b <- standardise_dynamic(OM_b$advice$Fadv["ANE",,,,,],map_b,proj_years)
    
    list(R=R_b,SSB=B_b,Catch=C_b,TAC=TAC_b,Fadv=Fadv_b)
  }
  
  # ==========================================================
  # EXTRACT COMPLETE MP ENSEMBLE
  # ==========================================================
  
  MP_blocks <- map(block_ids,extract_MP_block)
  
  R_MP <- map_dfr(MP_blocks,"R")
  B_MP <- map_dfr(MP_blocks,"SSB")
  C_MP <- map_dfr(MP_blocks,"Catch")
  TAC_MP <- map_dfr(MP_blocks,"TAC")
  Fadv_MP <- map_dfr(MP_blocks,"Fadv")
  
  rm(MP_blocks)
  gc()
  
  # ==========================================================
  # CHECK 2 - COMPLETE EXTRACTED ENSEMBLE
  # ==========================================================
  
  check_complete_MP <- bind_rows(validate_MP_output(R_MP,"Recruitment"),
                                 validate_MP_output(B_MP,"SSB"),
                                 validate_MP_output(C_MP,"Realised catch"),
                                 validate_MP_output(TAC_MP,"TAC advice"),
                                 validate_MP_output(Fadv_MP,"F advice"))
  
  stopifnot(all(check_complete_MP$n==30000))
  stopifnot(all(check_complete_MP$n_OM==100))
  stopifnot(all(check_complete_MP$n_trajectories==1000))
  stopifnot(all(check_complete_MP$first_year==first_proj))
  stopifnot(all(check_complete_MP$last_year==last_proj))
  stopifnot(all(check_complete_MP$n_NA==0))
  stopifnot(all(check_complete_MP$n_nonfinite==0))
  stopifnot(all(check_complete_MP$n_negative==0))
  for (dat in list(R_MP, B_MP, C_MP, TAC_MP, Fadv_MP)) {
    stopifnot(!anyDuplicated(dat[c("year", "iter1000")]))
    stopifnot(all(table(dat$iter1000)==30))
  }
  stopifnot(all(Fadv_MP$value<=Fcap+1e-8))

  
  cat("CHECK 2 PASSED - Complete 1000-trajectory MP ensemble extracted\n")
  
  # ==========================================================
  # IMPLEMENTATION DIAGNOSTICS
  # ==========================================================
  
  implementation_MP <- C_MP %>% rename(Catch=value) %>% 
                        left_join(TAC_MP %>% rename(TAC=value),by=c("year","iter1000","om","replicate")) %>% 
                        mutate(Catch_minus_TAC=Catch-TAC,
                               Catch_TAC_ratio=if_else(TAC>0,Catch/TAC,NA_real_))
  
  implementation_summary <- implementation_MP %>% summarise(Catch_TAC_diff_median=median(Catch_minus_TAC,na.rm=TRUE),
                                                            Catch_TAC_diff_p05=quantile(Catch_minus_TAC,0.05,na.rm=TRUE),
                                                            Catch_TAC_diff_p95=quantile(Catch_minus_TAC,0.95,na.rm=TRUE),
                                                            Catch_TAC_ratio_median=median(Catch_TAC_ratio,na.rm=TRUE),
                                                            Catch_TAC_ratio_p05=quantile(Catch_TAC_ratio,0.05,na.rm=TRUE),
                                                            Catch_TAC_ratio_p95=quantile(Catch_TAC_ratio,0.95,na.rm=TRUE))
  
  # ==========================================================
  # PERFORMANCE BY HORIZON
  # ==========================================================
  
  evaluate_horizon <- function(horizon_name,horizon_years) {
    
    yrs <- intersect(horizon_years,proj_years)
    
    B_h <- B_MP %>% filter(year %in% yrs)
    C_h <- C_MP %>% filter(year %in% yrs)
    Fadv_h <- Fadv_MP %>% filter(year %in% yrs)
    
    # --------------------------------------------------------
    # BIOLOGICAL RISK
    # --------------------------------------------------------
    
    risk_h <- B_h %>% mutate(below_Blim=value<Blim,below_Bpa=value<Bpa)
    
    risk_annual_h <- risk_h %>% group_by(year) %>% 
                      summarise(n=dplyr::n(),
                                P_Blim=mean(below_Blim),
                                P_Bpa=mean(below_Bpa),
                                SSB_median=median(value),
                                SSB_p05=quantile(value,0.05),
                                SSB_p95=quantile(value,0.95),.groups="drop")
    
    risk_trajectory_h <- risk_h %>% group_by(iter1000,om,replicate) %>% 
                          summarise(ever_below_Blim=any(below_Blim),
                                    ever_below_Bpa=any(below_Bpa),
                                    min_SSB=min(value),mean_SSB=mean(value),
                                    median_SSB=median(value),.groups="drop")
    
    risk_summary_h <- risk_trajectory_h %>% 
                        summarise(P_ever_Blim=mean(ever_below_Blim),
                                  P_ever_Bpa=mean(ever_below_Bpa),
                                  min_SSB=min(min_SSB))
    
    risk_performance_h <- risk_annual_h %>% 
                            summarise(max_P_Blim=max(P_Blim),
                                      n_years_P_Blim_gt_5=sum(P_Blim>risk_threshold),
                                      max_P_Bpa=max(P_Bpa),
                                      n_years_P_Bpa_gt_5=sum(P_Bpa>risk_threshold))
    
    # --------------------------------------------------------
    # FISHERY CLOSURES
    # --------------------------------------------------------
    
    closure_h <- Fadv_h %>% mutate(closed=value<=1e-8)
    
    closure_annual_h <- closure_h %>% group_by(year) %>% 
                        summarise(n=dplyr::n(),
                                  n_closed=sum(closed),
                                  P_closed=mean(closed),.groups="drop")
    
    closure_trajectory_h <- closure_h %>% group_by(iter1000,om,replicate) %>% 
                            summarise(n_years_closed=sum(closed),
                                      ever_closed=any(closed),.groups="drop")
    
    closure_summary_h <- closure_trajectory_h %>% 
                          summarise(P_ever_closed=mean(ever_closed),
                                    mean_years_closed=mean(n_years_closed),
                                    median_years_closed=median(n_years_closed),
                                    max_years_closed=max(n_years_closed))
    
    closure_performance_h <- closure_annual_h %>% 
                              summarise(mean_P_closed=mean(P_closed),
                                        max_P_closed=max(P_closed),
                                        year_max_P_closed=year[which.max(P_closed)])
    
    # --------------------------------------------------------
    # YIELD
    # --------------------------------------------------------
    
    yield_trajectory_h <- C_h %>% group_by(iter1000,om,replicate) %>% 
                          summarise(mean_catch=mean(value),
                                    median_catch=median(value),.groups="drop")
    
    yield_summary_h <- yield_trajectory_h %>% 
                      summarise(mean_mean_catch=mean(mean_catch),
                                median_mean_catch=median(mean_catch),
                                p05_mean_catch=quantile(mean_catch,0.05),
                                p95_mean_catch=quantile(mean_catch,0.95),
                                min_mean_catch=min(mean_catch),
                                max_mean_catch=max(mean_catch))
    
    yield_annual_h <- C_h %>% group_by(year) %>% 
                      summarise(mean_catch=mean(value),
                                median_catch=median(value),
                                p05_catch=quantile(value,0.05),
                                p95_catch=quantile(value,0.95),.groups="drop")
    
    # --------------------------------------------------------
    # YIELD IN OPEN YEARS
    # --------------------------------------------------------
    
    yield_open_h <- C_h %>% left_join(closure_h %>% dplyr::select(year,iter1000,om,replicate,closed),
                                      by=c("year","iter1000","om","replicate")) %>% filter(!closed)
    
    yield_open_trajectory_h <- yield_open_h %>% group_by(iter1000,om,replicate) %>% 
                                summarise(mean_catch_open=mean(value),
                                          median_catch_open=median(value),
                                          n_years_open=dplyr::n(),.groups="drop")
    
    yield_open_summary_h <- yield_open_trajectory_h %>%
      summarise(n_trajectories_open=n_distinct(iter1000),
                P_trajectories_open=n_distinct(iter1000)/1000,
                mean_mean_catch_open=mean(mean_catch_open),
                median_mean_catch_open=median(mean_catch_open),
                p05_mean_catch_open=quantile(mean_catch_open,0.05),
                p95_mean_catch_open=quantile(mean_catch_open,0.95),
                median_median_catch_open=median(median_catch_open),
                p05_median_catch_open=quantile(median_catch_open,0.05),
                p95_median_catch_open=quantile(median_catch_open,0.95))
    
    # --------------------------------------------------------
    # SSB IN OPEN YEARS
    # --------------------------------------------------------
    
    SSB_open_h <- B_h %>% left_join(closure_h %>% dplyr::select(year,iter1000,om,replicate,closed),
                                    by=c("year","iter1000","om","replicate")) %>% filter(!closed)
    
    SSB_open_trajectory_h <- SSB_open_h %>% group_by(iter1000,om,replicate) %>% 
                              summarise(mean_SSB_open=mean(value),
                                        median_SSB_open=median(value),.groups="drop")
    
    SSB_open_summary_h <- SSB_open_trajectory_h %>%
      summarise(n_trajectories_SSB_open=n_distinct(iter1000),
                P_trajectories_SSB_open=n_distinct(iter1000)/1000,
                median_mean_SSB_open=median(mean_SSB_open),
                p05_mean_SSB_open=quantile(mean_SSB_open,0.05),
                p95_mean_SSB_open=quantile(mean_SSB_open,0.95))
    
    # --------------------------------------------------------
    # INTERANNUAL CATCH VARIABILITY
    # --------------------------------------------------------
    
    IAV_trajectory_h <- C_h %>% arrange(iter1000,year) %>% group_by(iter1000,om,replicate) %>% 
                        mutate(Catch_prev=lag(value),
                               IAV_step=if_else(value==0 & Catch_prev==0,0,abs(value-Catch_prev)/((value+Catch_prev)/2))) %>% 
                                summarise(IAV=mean(IAV_step,na.rm=TRUE),.groups="drop")
    
    IAV_summary_h <- IAV_trajectory_h %>% summarise(mean_IAV=mean(IAV),
                                                    median_IAV=median(IAV),
                                                    p05_IAV=quantile(IAV,0.05),
                                                    p95_IAV=quantile(IAV,0.95),
                                                    min_IAV=min(IAV),
                                                    max_IAV=max(IAV))
    
    # --------------------------------------------------------
    # PERFORMANCE SUMMARY
    # --------------------------------------------------------
    
    performance_h <- tibble(scenario_number=scenario_number,
                            scenario_id=scenario_id,
                            Fcap=Fcap,
                            Besc=Besc,
                            Blim=Blim,
                            Bpa=Bpa,
                            horizon=horizon_name,
                            first_year=min(yrs),
                            last_year=max(yrs),
                            max_P_Blim=risk_performance_h$max_P_Blim,
                            n_years_P_Blim_gt_5=risk_performance_h$n_years_P_Blim_gt_5,
                            P_ever_Blim=risk_summary_h$P_ever_Blim,
                            max_P_Bpa=risk_performance_h$max_P_Bpa,
                            n_years_P_Bpa_gt_5=risk_performance_h$n_years_P_Bpa_gt_5,
                            P_ever_Bpa=risk_summary_h$P_ever_Bpa,
                            min_SSB=risk_summary_h$min_SSB,
                            mean_P_closed=closure_performance_h$mean_P_closed,
                            max_P_closed=closure_performance_h$max_P_closed,
                            P_ever_closed=closure_summary_h$P_ever_closed,
                            mean_years_closed=closure_summary_h$mean_years_closed,
                            median_years_closed=closure_summary_h$median_years_closed,
                            mean_catch=yield_summary_h$mean_mean_catch,
                            median_catch=yield_summary_h$median_mean_catch,
                            p05_catch=yield_summary_h$p05_mean_catch,
                            p95_catch=yield_summary_h$p95_mean_catch,
                            mean_catch_open=yield_open_summary_h$mean_mean_catch_open,
                            median_mean_catch_open=yield_open_summary_h$median_mean_catch_open,
                            p05_mean_catch_open=yield_open_summary_h$p05_mean_catch_open,
                            p95_mean_catch_open=yield_open_summary_h$p95_mean_catch_open,
                            median_mean_SSB_open=SSB_open_summary_h$median_mean_SSB_open,
                            p05_mean_SSB_open=SSB_open_summary_h$p05_mean_SSB_open,
                            p95_mean_SSB_open=SSB_open_summary_h$p95_mean_SSB_open,
                            mean_IAV=IAV_summary_h$mean_IAV,
                            median_IAV=IAV_summary_h$median_IAV,
                            p05_IAV=IAV_summary_h$p05_IAV,
                            p95_IAV=IAV_summary_h$p95_IAV)
    
    trajectory_h <- yield_trajectory_h %>% 
      left_join(yield_open_trajectory_h,by=c("iter1000","om","replicate")) %>% 
      left_join(IAV_trajectory_h,by=c("iter1000","om","replicate")) %>% 
      left_join(risk_trajectory_h,by=c("iter1000","om","replicate")) %>% 
      left_join(SSB_open_trajectory_h,by=c("iter1000","om","replicate")) %>%
      left_join(closure_trajectory_h,by=c("iter1000","om","replicate")) %>% 
      mutate(n_years_open=replace_na(n_years_open,0L)) %>%
      mutate(scenario_number=scenario_number,
             scenario_id=scenario_id,
             Fcap=Fcap,Besc=Besc,
             horizon=horizon_name,.before=1)
    
    annual_h <- risk_annual_h %>% 
                left_join(closure_annual_h %>% dplyr::select(year,P_closed),by="year") %>% 
                left_join(yield_annual_h,by="year") %>% 
                mutate(scenario_number=scenario_number,scenario_id=scenario_id,Fcap=Fcap,Besc=Besc,horizon=horizon_name,.before=1)
    
    list(performance=performance_h,trajectory=trajectory_h,annual=annual_h)
  }
  
  horizon_results <- imap(horizons,~evaluate_horizon(.y,.x))
  
  performance_MP <- map_dfr(horizon_results,"performance")
  trajectory_MP <- map_dfr(horizon_results,"trajectory")
  annual_MP <- map_dfr(horizon_results,"annual")
  
  # ==========================================================
  # STRUCTURAL FINAL CHECKS
  # ==========================================================
  
  stopifnot(nrow(performance_MP)==length(horizons))
  stopifnot(all(performance_MP$max_P_Blim>=0 & performance_MP$max_P_Blim<=1))
  stopifnot(all(performance_MP$P_ever_Blim>=0 & performance_MP$P_ever_Blim<=1))
  stopifnot(all(performance_MP$max_P_Bpa>=0 & performance_MP$max_P_Bpa<=1))
  stopifnot(all(performance_MP$P_ever_Bpa>=0 & performance_MP$P_ever_Bpa<=1))
  stopifnot(all(performance_MP$mean_catch>=0))
  stopifnot(all(performance_MP$mean_IAV>=0 & performance_MP$mean_IAV<=2))
  stopifnot(all(is.finite(trajectory_MP$mean_SSB)))
  stopifnot(all(is.finite(trajectory_MP$median_SSB)))
  stopifnot(all(trajectory_MP$mean_SSB>=0))
  stopifnot(all(trajectory_MP$median_SSB>=0))
  stopifnot(all(trajectory_MP$n_years_open+trajectory_MP$n_years_closed==case_when(trajectory_MP$horizon=="Full"~30L,TRUE~10L)))
  cat("PERFORMANCE EVALUATION COMPLETED\n")
  print(performance_MP %>% dplyr::select(horizon,max_P_Blim,P_ever_Blim,max_P_Bpa,P_ever_Bpa,mean_catch,median_IAV,P_ever_closed))
  
  list(performance=performance_MP,
       trajectory=trajectory_MP,
       annual=annual_MP,
       implementation=implementation_summary,
       validation=check_complete_MP)
}

# ============================================================
# EVALUATE THE 12 SELECTED MPs
# ============================================================

candidate_results <- map(seq_len(nrow(scenario_grid)),~evaluate_candidate_MP(scenario_grid[.x,,drop=FALSE]))

# ============================================================
# COMBINE RESULTS
# ============================================================

candidate_performance <- map_dfr(candidate_results,"performance")
candidate_trajectory_metrics <- map_dfr(candidate_results,"trajectory")
candidate_annual_metrics <- map_dfr(candidate_results,"annual")

candidate_implementation <- map2_dfr(candidate_results,seq_len(nrow(scenario_grid)),~.x$implementation %>% 
                                       mutate(scenario_number=scenario_grid$scenario_number[.y],
                                              scenario_id=scenario_grid$scenario_id[.y],
                                              Fcap=scenario_grid$Fcap[.y],
                                              Besc=scenario_grid$Besc[.y],.before=1))

candidate_validation <- map2_dfr(candidate_results,seq_len(nrow(scenario_grid)),~.x$validation %>% 
                                   mutate(scenario_number=scenario_grid$scenario_number[.y],
                                          scenario_id=scenario_grid$scenario_id[.y],
                                          Fcap=scenario_grid$Fcap[.y],
                                          Besc=scenario_grid$Besc[.y],.before=1))


# Preserve REF identities for subsequent paired comparisons.
add_robustness_keys <- function(dat) {
  dat %>% left_join(scenario_grid %>% dplyr::select(scenario_number, robustness_id,
                    reference_scenario_number, reference_scenario_id), by="scenario_number")
}
candidate_performance <- add_robustness_keys(candidate_performance)
candidate_trajectory_metrics <- add_robustness_keys(candidate_trajectory_metrics)
candidate_annual_metrics <- add_robustness_keys(candidate_annual_metrics)
candidate_implementation <- add_robustness_keys(candidate_implementation)
candidate_validation <- add_robustness_keys(candidate_validation)

# ============================================================
# FINAL CHECKS
# ============================================================

n_scenarios <- nrow(scenario_grid)

stopifnot(nrow(candidate_performance)==n_scenarios*length(horizons))
stopifnot(n_distinct(candidate_performance$scenario_number)==n_scenarios)
stopifnot(all(table(candidate_performance$scenario_number)==length(horizons)))
stopifnot(all(sort(unique(candidate_performance$horizon))==sort(names(horizons))))
stopifnot(n_distinct(candidate_trajectory_metrics$scenario_number)==n_scenarios)
stopifnot(n_distinct(candidate_annual_metrics$scenario_number)==n_scenarios)

cat("\nFINAL CHECKS PASSED\n")
cat("Candidate MPs evaluated:",n_scenarios,"\n")
cat("Performance rows:",nrow(candidate_performance),"\n")

# ============================================================
# SUMMARY - FULL HORIZON
# ============================================================

candidate_performance_full <- candidate_performance %>% filter(horizon=="Full") %>% arrange(Fcap,Besc)

print(candidate_performance_full %>% 
        dplyr::select(scenario_number,Fcap,Besc,max_P_Blim,P_ever_Blim,max_P_Bpa,P_ever_Bpa,mean_catch,median_IAV,P_ever_closed))

# ============================================================
# SAVE OUTPUTS
# ============================================================

performance_outputs <- list(robustness_id=robustness_id, horizons=horizons, risk_threshold=risk_threshold, source_validation_file=validation_file, scenario_grid=scenario_grid,performance=candidate_performance,performance_full=candidate_performance_full,trajectory_metrics=candidate_trajectory_metrics,annual_metrics=candidate_annual_metrics,implementation=candidate_implementation,validation=candidate_validation)

saveRDS(performance_outputs,file.path(performance_dir,"robustness_MP_performance.rds"))
write.csv(candidate_performance,file.path(performance_dir,"robustness_MP_performance_by_horizon.csv"),row.names=FALSE)
write.csv(candidate_performance_full,file.path(performance_dir,"robustness_MP_performance_full.csv"),row.names=FALSE)
write.csv(candidate_trajectory_metrics,file.path(performance_dir,"robustness_MP_trajectory_metrics.csv"),row.names=FALSE)
write.csv(candidate_annual_metrics,file.path(performance_dir,"robustness_MP_annual_metrics.csv"),row.names=FALSE)
write.csv(candidate_implementation,file.path(performance_dir,"robustness_MP_implementation_diagnostics.csv"),row.names=FALSE)
write.csv(candidate_validation,file.path(performance_dir,"robustness_MP_validation.csv"),row.names=FALSE)

cat("\nOutputs saved in:",performance_dir,"\n")

