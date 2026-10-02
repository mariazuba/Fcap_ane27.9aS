# ============================================================
# 21_CHECK_MP_SCENARIO_PROBABILISTIC.R
# Global post-run validation of candidate probabilistic
# Fcap + Besc management procedures
#
# Stock / analysis:
#   European anchovy ane.27.9aS, Gulf of Cádiz
#   Shortcut MSE implemented in FLBEIA + FLasher
#
# Role in the workflow:
#
#   20_RUN_MP_SCENARIO_PROBABILISTIC.R
#                     |
#                     v
#   21_CHECK_MP_SCENARIO_PROBABILISTIC.R
#                     |
#                     v
#   22_EVALUATE_CANDIDATE_MPS.R
#
#   This script validates the structural and operational
#   consistency of the completed candidate-MP simulations
#   before management performance is evaluated.
#
# Purpose:
#   Validate all candidate MPs defined in candidate_MP_grid.csv.
#
#   For each scenario and block, the script checks:
#
#     - existence of run, map and metadata files;
#     - block size and iteration structure;
#     - trajectory / forecast-CV pairing;
#     - output iteration dimensions;
#     - F advice;
#     - TAC;
#     - realised catch;
#     - TAC implementation;
#     - HCR response: closure, reduced F or Fcap.
#
#   Across blocks, for each candidate MP, it checks:
#
#     - exactly 1000 global trajectories;
#     - iter1000 = 1:1000;
#     - unique global iterations;
#     - 100 conditioned OMs;
#     - 10 stochastic replicates per OM;
#     - unique OM x replicate combinations;
#     - exact pairing with the master forecast-CV map.
#
# Inputs:
#   data/mse/scenarios/candidate_MP_grid.csv
#   data/mse/scenarios/forecast_CV_trajectory_map.csv
#
#   data/mse/MP_runs/<scenario_id>/
#     block_XX.rds
#     block_XX_map.csv
#     block_XX_metadata.csv
#
# Outputs:
#   outputs/mse/validation/
#     21_MP_validation_by_block.csv
#     21_MP_validation_by_scenario.csv
#     21_HCR_branch_summary.csv
#     21_TAC_implementation_summary.csv
#     21_TAC_implementation_flags.csv
#     21_largest_TAC_Catch_differences.csv
#     21_MP_validation_summary.rds
#
# Working Document outputs:
#   None directly.
#
#   The validation outputs document the implementation checks
#   performed before candidate-MP performance evaluation.
#
# Important:
#   This script DOES NOT rerun FLBEIA.
#
#   Catch/TAC differences are implementation diagnostics.
#   A Catch/TAC difference > 5% is flagged but does not by
#   itself define structural failure of an MSE run.
#
# Expected production design:
#   20 candidate MPs
#   10 blocks per MP
#   100 trajectories per block
#   1000 trajectories per MP
#   100 conditioned OMs
#   10 replicates per OM
#   Projection period: 2025-2054
# ============================================================


# ============================================================
# 0. CLEAN WORKSPACE AND LOAD PACKAGES
# ============================================================

rm(list=ls())

library(FLCore)
library(FLBEIA)
library(dplyr)
library(tidyr)
library(here)


# ============================================================
# 1. SETTINGS AND EXPECTED PRODUCTION DESIGN
# ============================================================

proj.yr <- 2025
last.yr <- 2054

expected_blocks       <- 10
expected_block_size   <- 100
expected_trajectories <- 1000
expected_OMs          <- 100
expected_replicates   <- 10

tol_F         <- 1e-8
tol_TAC_ratio <- 0.05
tol_CV        <- 1e-12


# ============================================================
# 2. DIRECTORIES AND INPUT FILES
# ============================================================

scenario_dir <- here("data","mse","scenarios")
run_root     <- here("data","mse","MP_runs")
output_dir   <- here("outputs","mse","validation")

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

scenario_file     <- file.path(scenario_dir, "candidate_MP_grid.csv")
forecast_map_file <- file.path(scenario_dir,"forecast_CV_trajectory_map.csv")

stopifnot(file.exists(scenario_file))
stopifnot(file.exists(forecast_map_file))


# ============================================================
# 3. LOAD EXPERIMENTAL DESIGN
# ============================================================

scenario_table <- read.csv(scenario_file,stringsAsFactors=FALSE)
forecast_map   <- read.csv(forecast_map_file,stringsAsFactors=FALSE)

required_scenario_cols <- c("scenario_number","scenario_id","Fcap","Besc","forecast_alpha")

stopifnot(all(required_scenario_cols %in% names(scenario_table)))


# ============================================================
# 4. VALIDATE MASTER FORECAST-CV MAP
# ============================================================
#
# The forecast-CV map defines the empirical S2 forecast
# uncertainty assigned to each global trajectory.
#
# The same global trajectory-to-CV pairing must be retained
# across all candidate MPs.
# ============================================================

if("global_iter" %in% names(forecast_map)){forecast_map <- forecast_map %>%rename(iter1000=global_iter)}

required_forecast_cols <- c("iter1000","forecast_CV")

stopifnot(all(required_forecast_cols %in% names(forecast_map)))

forecast_map <- forecast_map %>%arrange(iter1000)

stopifnot(nrow(forecast_map)==expected_trajectories)
stopifnot(!anyDuplicated(forecast_map$iter1000))
stopifnot(identical(forecast_map$iter1000,seq_len(expected_trajectories)))
stopifnot(all(is.finite(forecast_map$forecast_CV)))
stopifnot(all(forecast_map$forecast_CV>0))

# ============================================================
# 5. INITIALISE VALIDATION STORAGE
# ============================================================

block_validation_list    <- list()
scenario_validation_list <- list()
branch_summary_list      <- list()
TAC_summary_list         <- list()
TAC_flags_list           <- list()
largest_TAC_differences_list <- list()

block_counter <- 0

cat("\n")
cat("============================================================\n")
cat("GLOBAL POST-RUN PROBABILISTIC MP VALIDATION\n")
cat("============================================================\n")
cat("Candidate MPs:",nrow(scenario_table),"\n")
cat("Expected blocks per MP:",expected_blocks,"\n")
cat("Expected trajectories per block:",expected_block_size,"\n")
cat("Expected trajectories per MP:",expected_trajectories,"\n")
cat("Projection years:",proj.yr,"-",last.yr,"\n")
cat("============================================================\n")


# ============================================================
# 6. LOOP OVER CANDIDATE MPs
# ============================================================

for(i in seq_len(nrow(scenario_table))){
  
  sc <- scenario_table[i,,drop=FALSE]
  
  scenario_number <- sc$scenario_number
  scenario_id <- sc$scenario_id
  Fcap <- sc$Fcap
  Besc <- sc$Besc
  forecast_alpha <- sc$forecast_alpha
  
  cat("\n")
  cat("============================================================\n")
  cat("SCENARIO",scenario_number,"OF",nrow(scenario_table),"\n")
  cat("Scenario ID:",scenario_id,"\n")
  cat("Fcap:",Fcap,"\n")
  cat("Besc:",Besc,"t\n")
  cat("Forecast alpha:",forecast_alpha,"\n")
  cat("============================================================\n")
  
  run_dir <- file.path(run_root,scenario_id)
  
  if(!dir.exists(run_dir)){stop("Run directory does not exist for scenario ",scenario_number,": ",run_dir)}
  
  scenario_map_list <- list()
  scenario_results_list <- list()
  
  
  # ==========================================================
  # 6.1 LOOP OVER PRODUCTION BLOCKS
  # ==========================================================
  
  for(block in seq_len(expected_blocks)){
  
    cat("Checking block ",sprintf("%02d",block),"... ")
    
    OM_file <- file.path(run_dir,sprintf("block_%02d.rds",block))
    map_file <- file.path(run_dir,sprintf("block_%02d_map.csv",block))
    metadata_file <- file.path(run_dir,sprintf("block_%02d_metadata.csv",block))
    
    if(!file.exists(OM_file)){stop("Missing OM file: ", OM_file)}
    if(!file.exists(map_file)){stop("Missing map file: ", map_file)}
    if(!file.exists(metadata_file)){stop("Missing metadata file: ",metadata_file)}
    
    
    # ========================================================
    # 6.2 LOAD BLOCK OUTPUTS
    # ========================================================
    
    OM <- readRDS(OM_file)
    run_map <- read.csv(map_file,stringsAsFactors=FALSE)
    metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
    required_map_cols <- c("iter1000","om","replicate","iter_local","forecast_CV")
    
    stopifnot(all(required_map_cols %in% names(run_map)))
    
    nb <- nrow(run_map)
    
    
    # ========================================================
    # 6.3 VALIDATE BLOCK METADATA
    # ========================================================
    #
    # Metadata are checked only for fields that are explicitly
    # present in the metadata file. The check confirms that
    # identifiers and MP settings stored by the production run
    # correspond to the candidate scenario being validated.
    # ========================================================
  
    stopifnot(nrow(metadata)>=1)
    if("scenario_number" %in% names(metadata)){stopifnot(all(metadata$scenario_number==scenario_number))}
    if("scenario_id" %in% names(metadata)){stopifnot(all(metadata$scenario_id==scenario_id))}
    if("block" %in% names(metadata)){stopifnot(all(metadata$block==block))}
    if("Fcap" %in% names(metadata)){stopifnot(all(abs(metadata$Fcap-Fcap)<=tol_F))}
    if("Besc" %in% names(metadata)){stopifnot(all(metadata$Besc==Besc))}
    if("forecast_alpha" %in% names(metadata)){stopifnot(all(metadata$forecast_alpha==forecast_alpha))}
    
    
    # ========================================================
    # 6.4 VALIDATE BLOCK ITERATION MAP
    # ========================================================
    
    stopifnot(nb==expected_block_size)
    stopifnot(n_distinct(run_map$iter_local)==nb)
    stopifnot(all(run_map$iter_local==seq_len(nb)))
    stopifnot(!anyDuplicated(run_map$iter1000))
    stopifnot(!anyDuplicated(run_map[,c("om","replicate")]))
    stopifnot(all(is.finite(run_map$forecast_CV)))
    stopifnot(all(run_map$forecast_CV>0))
    
    # ========================================================
    # 6.5 VALIDATE BLOCK FORECAST-CV PAIRING
    # ========================================================
    
    master_CV_block <- forecast_map %>%
      filter(iter1000 %in% run_map$iter1000) %>%
      select(iter1000,forecast_CV_master=forecast_CV)
    
    CV_check <- run_map %>%
      select(iter1000,forecast_CV) %>%
      left_join(master_CV_block,by="iter1000") %>%
      mutate(CV_difference=forecast_CV-forecast_CV_master)
    
    stopifnot(nrow(CV_check)==nb)
    stopifnot(!any(is.na(CV_check$forecast_CV_master)))
    stopifnot(all(abs(CV_check$CV_difference)<=tol_CV))
    
    
    # ========================================================
    # 6.6 VALIDATE OUTPUT ITERATION DIMENSION
    # ========================================================
    
    ni_output <- dim(OM$biols$ANE@n)[6]
    stopifnot(ni_output==nb)
    
    # ========================================================
    # 6.7 EXTRACT F ADVICE
    # ========================================================
    
    years <- as.character(proj.yr:last.yr)
    Fadv_array <- OM$advice$Fadv["ANE",years,,,,]
    Fadv <- expand.grid(year=proj.yr:last.yr,iter_local=seq_len(nb),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
    Fadv$Fadv <- as.numeric(Fadv_array)
    
    stopifnot(nrow(Fadv)==length(years)*nb)
    stopifnot(all(is.finite(Fadv$Fadv)))
    stopifnot(all(Fadv$Fadv>=0))
    stopifnot(all(Fadv$Fadv<=Fcap+tol_F))
    
    
    # ========================================================
    # 6.8 EXTRACT TAC
    # ========================================================
    
    TAC_array <- OM$advice$TAC["ANE",years,,,,]
    
    TAC <- expand.grid(
      year=proj.yr:last.yr,
      iter_local=seq_len(nb),
      KEEP.OUT.ATTRS=FALSE,
      stringsAsFactors=FALSE)
    
    TAC$TAC <- as.numeric(TAC_array)
    
    stopifnot(nrow(TAC)==length(years)*nb)
    stopifnot(all(is.finite(TAC$TAC)))
    stopifnot(all(TAC$TAC>=0))
    
    
    # ========================================================
    # 6.9 EXTRACT REALISED CATCH
    # ========================================================
    
    cobj <- OM$fleets$SEINE@metiers$ALL@catches$ANE
    
    Catch_annual <- seasonSums(quantSums(cobj@landings))
    Catch_array <- Catch_annual[,years,,,,]
    
    Catch <- expand.grid(
      year=proj.yr:last.yr,
      iter_local=seq_len(nb),
      KEEP.OUT.ATTRS=FALSE,
      stringsAsFactors=FALSE)
    
    Catch$Catch <- as.numeric(Catch_array)
    
    stopifnot(nrow(Catch)==length(years)*nb)
    stopifnot(all(is.finite(Catch$Catch)))
    stopifnot(all(Catch$Catch>=0))
    
    # ========================================================
    # 6.10 COMBINE YEAR × TRAJECTORY RESULTS
    # ========================================================
    
    results <- Fadv %>%
      left_join(TAC,by=c("year","iter_local")) %>%
      left_join(Catch,by=c("year","iter_local")) %>%
      left_join(run_map %>%select(iter_local,iter1000,om,replicate,forecast_CV),by="iter_local")
    
    stopifnot(nrow(results)==length(years)*nb)
    stopifnot(!anyDuplicated(results[,c("year","iter_local")]))
    
    # ========================================================
    # 6.11 DIAGNOSE TAC IMPLEMENTATION
    # ========================================================
    #
    # Catch/TAC differences are implementation diagnostics.
    # The 5% threshold is used to flag observations for review;
    # it is not used as a structural acceptance criterion.
    # ========================================================
    
    results <- results %>%
      mutate(
        Catch_minus_TAC=Catch-TAC,
        Catch_TAC_ratio=ifelse(TAC>0,Catch/TAC,NA_real_))
    
    
    # ========================================================
    # 6.12 CLASSIFY HCR RESPONSE
    # ========================================================
    
    results <- results %>%
      mutate(
        HCR_branch=case_when(
          abs(Fadv)<=tol_F ~ "closure",
          abs(Fadv-Fcap)<=tol_F ~ "Fcap",
          Fadv>tol_F & Fadv<Fcap-tol_F ~ "reduced",
          TRUE ~ "other"))
    
    stopifnot(!any(results$HCR_branch=="other") )
    
    # ========================================================
    # 6.13 STRUCTURAL OUTPUT CHECKS
    # ========================================================
    
    stopifnot(all(is.finite(results$Fadv)))
    stopifnot(all(results$Fadv>=0))
    stopifnot(all(results$Fadv<=Fcap+tol_F))
    stopifnot(all(is.finite(results$TAC)))
    stopifnot(all(results$TAC>=0))
    stopifnot(all(is.finite(results$Catch)))
    stopifnot(all(results$Catch>=0))
    
    # ========================================================
    # 6.14 IDENTIFY TAC IMPLEMENTATION FLAGS
    # ========================================================
    
    TAC_flags <- results %>%filter(TAC>0) %>%
      filter(abs(Catch_TAC_ratio-1)>tol_TAC_ratio) %>%
      mutate(
        scenario_number=scenario_number,
        scenario_id=scenario_id,
        block=block) %>%
      select(scenario_number,scenario_id,block,everything())
    
    
    # ========================================================
    # 6.15 IDENTIFY LARGEST TAC-CATCH DIFFERENCES
    # ========================================================
    
    largest_TAC_differences <- results %>%filter(TAC>0) %>%
      arrange(desc(abs(Catch_minus_TAC))) %>%
      slice_head(n=10) %>%
      mutate(
        scenario_number=scenario_number,
        scenario_id=scenario_id,
        block=block) %>%
      select(
        scenario_number,
        scenario_id,
        block,
        year,
        iter_local,
        iter1000,
        om,
        replicate,
        forecast_CV,
        HCR_branch,
        Fadv,
        TAC,
        Catch,
        Catch_minus_TAC,
        Catch_TAC_ratio)
    
    
    # ========================================================
    # 6.16 BUILD BLOCK VALIDATION SUMMARY
    # ========================================================
    
    block_validation <- tibble(
      scenario_number=scenario_number,
      scenario_id=scenario_id,
      block=block,
      n_trajectories=nb,
      n_OM=n_distinct(run_map$om),
      min_iter1000=min(run_map$iter1000),
      max_iter1000=max(run_map$iter1000),
      n_iter1000=n_distinct(run_map$iter1000),
      min_forecast_CV=min(run_map$forecast_CV),
      median_forecast_CV=median(run_map$forecast_CV),
      max_forecast_CV=max(run_map$forecast_CV),
      max_abs_CV_difference=max(abs(CV_check$CV_difference)),
      min_Fadv=min(results$Fadv),
      max_Fadv=max(results$Fadv),
      min_TAC=min(results$TAC),
      max_TAC=max(results$TAC),
      min_Catch=min(results$Catch),
      max_Catch=max(results$Catch),
      n_closure=sum(results$HCR_branch=="closure"),
      n_reduced=sum(results$HCR_branch=="reduced"),
      n_Fcap=sum(results$HCR_branch=="Fcap"),
      n_TAC_flags=nrow(TAC_flags),
      structural_PASS=TRUE
    )
    
    
    # ========================================================
    # 6.17 STORE BLOCK RESULTS
    # ========================================================
    
    block_counter <- block_counter+1
    block_validation_list[[block_counter]] <- block_validation
    scenario_map_list[[block]] <- run_map %>%mutate(block=block)
    scenario_results_list[[block]] <- results %>%mutate(block=block)
    
    if(nrow(TAC_flags)>0){TAC_flags_list[[length(TAC_flags_list)+1]] <- TAC_flags}
    
    largest_TAC_differences_list[[length(largest_TAC_differences_list)+1]] <- largest_TAC_differences
    
    cat("PASS\n")
    
    rm(
      OM,
      run_map,
      metadata,
      Fadv_array,
      TAC_array,
      Catch_array,
      Catch_annual,
      Fadv,
      TAC,
      Catch,
      results,
      CV_check,
      master_CV_block,
      TAC_flags,
      largest_TAC_differences)
    
    gc(verbose=FALSE)
  }
  
  
  # ==========================================================
  # 7. COMBINE BLOCKS WITHIN EACH CANDIDATE MP
  # ==========================================================
  
  scenario_map <- bind_rows(scenario_map_list)
  scenario_results <- bind_rows(scenario_results_list)
  
  # ==========================================================
  # 8. VALIDATE GLOBAL TRAJECTORY STRUCTURE
  # ==========================================================
  
  stopifnot(nrow(scenario_map)==expected_trajectories)
  stopifnot(n_distinct(scenario_map$iter1000)==expected_trajectories)
  stopifnot(!anyDuplicated(scenario_map$iter1000))
  stopifnot(identical(sort(scenario_map$iter1000),seq_len(expected_trajectories)))
  
  stopifnot(n_distinct(scenario_map$om)==expected_OMs)
  
  OM_replication <- scenario_map %>%count(om,name="n_replicates")
  
  stopifnot(nrow(OM_replication)==expected_OMs)
  stopifnot(all(OM_replication$n_replicates==expected_replicates))
  stopifnot(!anyDuplicated(scenario_map[,c("om","replicate")]))
  
  OM_replicate_check <- scenario_map %>%count(om,replicate,name="n")
  
  stopifnot(nrow(OM_replicate_check)==expected_trajectories)
  stopifnot(all(OM_replicate_check$n==1))
  
  
  # ==========================================================
  # 9. VALIDATE GLOBAL FORECAST-CV PAIRING
  # ==========================================================
  
  scenario_CV_check <- scenario_map %>%
    select(iter1000,forecast_CV) %>%
    arrange(iter1000) %>%
    left_join(forecast_map %>% select(iter1000,forecast_CV_master=forecast_CV),by="iter1000") %>%
    mutate(CV_difference=forecast_CV-forecast_CV_master)
  
  stopifnot(nrow(scenario_CV_check)==expected_trajectories)
  stopifnot(!any(is.na(scenario_CV_check$forecast_CV_master)))
  stopifnot(all(abs(scenario_CV_check$CV_difference)<=tol_CV))
  
  
  # ==========================================================
  # 10. VALIDATE YEAR × TRAJECTORY STRUCTURE
  # ==========================================================
  
  expected_year_trajectory_rows <- expected_trajectories*length(proj.yr:last.yr)
  
  stopifnot(nrow(scenario_results)==expected_year_trajectory_rows)
  
  trajectory_year_check <- scenario_results %>%count(iter1000,name="n_years")
  
  stopifnot(nrow(trajectory_year_check)==expected_trajectories)
  stopifnot(all(trajectory_year_check$n_years==length(proj.yr:last.yr)))
  
  
  # ==========================================================
  # 11. SUMMARISE HCR RESPONSE
  # ==========================================================
  #
  # This is a descriptive implementation diagnostic.
  # Management performance is evaluated in script 22.
  # ==========================================================
  
  branch_summary <- scenario_results %>%
    count(HCR_branch,name="n") %>%
    mutate(
      percent=100*n/sum(n),
      scenario_number=scenario_number,
      scenario_id=scenario_id) %>%
    select(scenario_number,scenario_id, HCR_branch,n,percent)
  
  branch_summary_list[[i]] <- branch_summary
  
  
  # ==========================================================
  # 12. SUMMARISE TAC IMPLEMENTATION
  # ==========================================================
  
  TAC_summary <- scenario_results %>%
    summarise(
      median_difference=median(Catch_minus_TAC),
      min_difference=min(Catch_minus_TAC),
      max_difference=max(Catch_minus_TAC),
      median_ratio=median(Catch_TAC_ratio,na.rm=TRUE),
      min_ratio=min(Catch_TAC_ratio,na.rm=TRUE),
      max_ratio=max(Catch_TAC_ratio,na.rm=TRUE),
      n_TAC_positive=sum(TAC>0),
      n_TAC_flags=sum(TAC>0 &abs(Catch_TAC_ratio-1)>tol_TAC_ratio,na.rm=TRUE)) %>%
    mutate(scenario_number=scenario_number,
           scenario_id=scenario_id) %>%
    select(scenario_number,scenario_id,everything())
  
  TAC_summary_list[[i]] <- TAC_summary
  
  
  # ==========================================================
  # 13. BUILD SCENARIO VALIDATION SUMMARY
  # ==========================================================
  
  scenario_validation <- tibble(
    scenario_number=scenario_number,
    scenario_id=scenario_id,
    Fcap=Fcap,
    Besc=Besc,
    forecast_alpha=forecast_alpha,
    n_blocks=n_distinct(scenario_map$block),
    n_trajectories=nrow(scenario_map),
    n_iter1000=n_distinct(scenario_map$iter1000),
    n_OM=n_distinct(scenario_map$om),
    min_replicates_per_OM=min(OM_replication$n_replicates),
    max_replicates_per_OM=max(OM_replication$n_replicates),
    n_year_trajectory=nrow(scenario_results),
    min_forecast_CV=min(scenario_map$forecast_CV),
    median_forecast_CV=median(scenario_map$forecast_CV),
    max_forecast_CV=max(scenario_map$forecast_CV),
    max_abs_CV_difference=max(abs(scenario_CV_check$CV_difference)),
    min_Fadv=min(scenario_results$Fadv),
    max_Fadv=max(scenario_results$Fadv),
    n_closure=sum(scenario_results$HCR_branch=="closure"),
    n_reduced=sum(scenario_results$HCR_branch=="reduced"),
    n_Fcap=sum(scenario_results$HCR_branch=="Fcap"),
    n_TAC_flags=TAC_summary$n_TAC_flags,
    structural_PASS=TRUE)
  
  scenario_validation_list[[i]] <- scenario_validation
  
  cat("\n")
  cat("Scenario structural validation: PASS\n")
  cat("Trajectories:",nrow(scenario_map),"\n")
  cat("OMs:",n_distinct(scenario_map$om),"\n")
  cat("Replicates per OM:",unique(OM_replication$n_replicates),"\n")
  cat("Year x trajectory observations:",nrow(scenario_results),"\n")
  cat("Closure observations:",scenario_validation$n_closure,"\n")
  cat("Reduced-F observations:",scenario_validation$n_reduced,"\n")
  cat("Fcap observations:",scenario_validation$n_Fcap,"\n")
  cat("Catch/TAC flags >5%:",scenario_validation$n_TAC_flags,"\n")
  
  rm(
    scenario_map,
    scenario_results,
    scenario_map_list,
    scenario_results_list,
    OM_replication,
    OM_replicate_check,
    scenario_CV_check,
    trajectory_year_check,
    branch_summary,
    TAC_summary,
    scenario_validation)
  
  gc(verbose=FALSE)
}


# ============================================================
# 14. COMBINE GLOBAL VALIDATION OUTPUTS
# ============================================================

block_validation <- bind_rows(block_validation_list)
scenario_validation <- bind_rows(scenario_validation_list)
branch_summary <- bind_rows(branch_summary_list)
TAC_summary <- bind_rows(TAC_summary_list)

if(length(TAC_flags_list)>0){
  TAC_flags <- bind_rows(TAC_flags_list)
}else{
  TAC_flags <- tibble()
}

largest_TAC_differences <- bind_rows(largest_TAC_differences_list)


# ============================================================
# 15. FINAL GLOBAL VALIDATION CHECKS
# ============================================================
stopifnot(nrow(scenario_validation)==nrow(scenario_table))
stopifnot(nrow(block_validation)==nrow(scenario_table)*expected_blocks)
stopifnot(all(block_validation$structural_PASS))
stopifnot(all(scenario_validation$structural_PASS))
stopifnot(all(block_validation$n_trajectories==expected_block_size))
stopifnot(all(scenario_validation$n_trajectories==expected_trajectories))
stopifnot(all(scenario_validation$n_OM==expected_OMs))
stopifnot(all(scenario_validation$min_replicates_per_OM==expected_replicates))
stopifnot(all(scenario_validation$max_replicates_per_OM==expected_replicates))
stopifnot(all(scenario_validation$max_Fadv<=scenario_validation$Fcap+tol_F))


# ============================================================
# 16. SAVE VALIDATION OUTPUTS
# ============================================================

write.csv(block_validation,file.path(output_dir,"21_MP_validation_by_block.csv"),row.names=FALSE)
write.csv(scenario_validation,file.path(output_dir,"21_MP_validation_by_scenario.csv"),row.names=FALSE)
write.csv(branch_summary,file.path(output_dir,"21_HCR_branch_summary.csv"),row.names=FALSE)
write.csv(TAC_summary,file.path(output_dir,"21_TAC_implementation_summary.csv"),row.names=FALSE)
write.csv(TAC_flags,file.path(output_dir,"21_TAC_implementation_flags.csv"),row.names=FALSE)
write.csv(largest_TAC_differences,file.path(output_dir,"21_largest_TAC_Catch_differences.csv"),row.names=FALSE)

validation_object <- list(
  settings=list(
    proj.yr=proj.yr,
    last.yr=last.yr,
    expected_blocks=expected_blocks,
    expected_block_size=expected_block_size,
    expected_trajectories=expected_trajectories,
    expected_OMs=expected_OMs,
    expected_replicates=expected_replicates,
    tol_F=tol_F,
    tol_TAC_ratio=tol_TAC_ratio,
    tol_CV=tol_CV),
  block_validation=block_validation,
  scenario_validation=scenario_validation,
  branch_summary=branch_summary,
  TAC_summary=TAC_summary,
  TAC_flags=TAC_flags,
  largest_TAC_differences=largest_TAC_differences)

saveRDS(validation_object,file.path(output_dir,"21_MP_validation_summary.rds"))


# ============================================================
# 17. FINAL VALIDATION SUMMARY
# ============================================================
cat("\n")
cat("============================================================\n")
cat("GLOBAL POST-RUN VALIDATION COMPLETED\n")
cat("============================================================\n")
cat("Candidate MPs checked:",nrow(scenario_validation),"\n")
cat("Blocks checked:", nrow(block_validation),"\n")
cat("Trajectories per MP:", expected_trajectories,"\n")
cat("Projection years:",length(proj.yr:last.yr),"\n")
cat("Year x trajectory observations per MP:",expected_trajectories*length(proj.yr:last.yr),"\n")
cat("Total MP x year x trajectory observations:",nrow(scenario_validation)*expected_trajectories*length(proj.yr:last.yr),"\n")
cat("Structural failures: 0\n")
cat("All structural checks: PASS\n")
cat("============================================================\n")
cat("\n")
cat("Catch/TAC diagnostic flags (>5%): ")
cat(nrow(TAC_flags),"\n")
cat("\n")
cat("Validation outputs saved in:\n")
cat(output_dir,"\n")
cat("============================================================\n")


