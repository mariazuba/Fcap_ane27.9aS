# ============================================================
# 28_VALIDATE_ROBUSTNESS_RUNS.R
# Post-run validation of selected MPs under robustness scenarios
#
# Uses the structural and TAC diagnostics from scripts 21/21c.
# Adds REF map pairing, mandatory scenario metadata, saved input
# checks, and finite/non-negative Q2 SSB and Q3 recruitment checks.
# No simulations, performance ranking or input changes are made.
# Saved input_checks document runner checks; they do not independently
# prove identical AE draws or standard-normal recruitment deviates.
#
# Usage from the project root (R console):
#   source("scripts/28_VALIDATE_ROBUSTNESS_RUNS.R") # R-OM by default
#   Sys.setenv(MSE_VALIDATE_ROBUSTNESS_ID="S-C1")
#   source("scripts/28_VALIDATE_ROBUSTNESS_RUNS.R")
# Requires all 12 MPs, 10 production blocks each, and REF block maps.
# Outputs: outputs/mse/robustness/validation/post_run/<scenario>/
# Catch/TAC deviations >5% are review flags, not structural failures.
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

robustness_id <- Sys.getenv("MSE_VALIDATE_ROBUSTNESS_ID", "R-OM")

stopifnot(robustness_id %in% c("R-OM", "S-C1", "S-C2", "I0", "I1"))
run_root <- here("data", "mse", "robustness_runs")
reference_root <- here("data", "mse", "MP_runs")
output_dir <- here("outputs", "mse", "robustness", "validation", "post_run", robustness_id)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
scenario_file <- here("outputs", "mse", "robustness", "design", "Table26_02_MP_robustness_design.csv")
forecast_map_file <- here("data", "mse", "scenarios", "forecast_CV_trajectory_map.csv")
stopifnot(file.exists(scenario_file), file.exists(forecast_map_file))
design <- read.csv(scenario_file, stringsAsFactors=FALSE)
required <- c("robustness_id", "reference_scenario_number", "reference_scenario_id",
              "robustness_run_id", "Fcap", "Besc", "annual_risk_threshold",
              "recruitment_PE_mode", "OM_seasonal_share_mode", "assessment_error_SSB_F", "forecast_uncertainty")
stopifnot(all(required %in% names(design)))
scenario_table <- design[design$robustness_id==robustness_id, , drop=FALSE]
stopifnot(nrow(scenario_table)==12, !anyDuplicated(scenario_table$reference_scenario_id))
stopifnot(setequal(scenario_table$reference_scenario_number, c(9,10,11,12,17,18,19,20,25,26,27,28)))
scenario_table$forecast_alpha <- scenario_table$annual_risk_threshold
scenario_table$scenario_number <- scenario_table$reference_scenario_number
scenario_table$scenario_id <- scenario_table$robustness_run_id
forecast_map <- read.csv(forecast_map_file, stringsAsFactors=FALSE)

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
if (robustness_id %in% c("I0", "I1")) forecast_map$forecast_CV <- 0

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
    
    ref_file <- file.path(reference_root, sc$reference_scenario_id, sprintf("block_%02d_map.csv", block))
    input_file <- file.path(run_dir, sprintf("block_%02d_input_validation.rds", block))
    stopifnot(file.exists(ref_file), file.exists(input_file))
    ref_map <- read.csv(ref_file, stringsAsFactors=FALSE)
    paired_cols <- c("iter1000", "om", "replicate", "bootstrap", "a", "b", "sigmaR", "iter_local")
    stopifnot(all(c(paired_cols, "forecast_CV_reference", "sigmaR_used") %in% names(run_map)))
    stopifnot(all(c(paired_cols, "forecast_CV") %in% names(ref_map)))
    stopifnot(nrow(ref_map)==nrow(run_map))
    stopifnot(identical(run_map$iter1000, ref_map$iter1000))
    for (column in paired_cols) {
      stopifnot(isTRUE(all.equal(run_map[[column]], ref_map[[column]], tolerance=tol_CV)))
    }
    stopifnot(all(abs(run_map$forecast_CV_reference-ref_map$forecast_CV)<=tol_CV))
    stopifnot(all(is.finite(run_map$sigmaR_used)), all(run_map$sigmaR_used>=0))
    if (robustness_id=="R-OM") stopifnot(all(abs(run_map$sigmaR_used-run_map$sigmaR)<=tol_CV))
    input_check <- readRDS(input_file)
    stopifnot(identical(input_check$input_checks, "PASS"))
    stopifnot(identical(input_check$execution_mode, "full"))
    stopifnot(identical(input_check$robustness_id, robustness_id))
    stopifnot(identical(as.character(input_check$reference_scenario_id), as.character(sc$reference_scenario_id)))
    stopifnot(isTRUE(all.equal(input_check$run_map, run_map, check.attributes=FALSE, tolerance=tol_CV)))
    required_meta <- c("robustness_id", "reference_scenario_id", "execution_mode", "scenario_id",
                       "Fcap", "Besc", "block", "n_trajectories", "sigmaR_scenario",
                       "seasonal_exploitation", "assessment_error", "forecast_uncertainty")
    stopifnot(nrow(metadata)==1, all(required_meta %in% names(metadata)))
    stopifnot(metadata$robustness_id==robustness_id, metadata$reference_scenario_id==sc$reference_scenario_id)
    stopifnot(metadata$execution_mode=="full", metadata$n_trajectories==expected_block_size)
    stopifnot(metadata$sigmaR_scenario==sc$recruitment_PE_mode)
    stopifnot(metadata$seasonal_exploitation==sc$OM_seasonal_share_mode)
    stopifnot(metadata$assessment_error==sc$assessment_error_SSB_F)
    stopifnot(metadata$forecast_uncertainty==sc$forecast_uncertainty)

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
    stopifnot(all(run_map$forecast_CV>=0))
    
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
    biol <- OM$biols$ANE
    stopifnot(all(as.character(proj.yr:last.yr) %in% dimnames(biol@n)$year))
    stopifnot(identical(dimnames(biol@n)$season, as.character(1:4)))
    Mat <- predict(biol@mat)
    SSB_q2 <- quantSums(biol@n[,,,"2",,]*biol@wt[,,,"2",,]*Mat[,,,"2",,])
    SSB_values <- as.numeric(SSB_q2[,as.character(proj.yr:last.yr),,,,])
    R_values <- as.numeric(biol@n["0",as.character(proj.yr:last.yr),,"3",,])
    stopifnot(length(SSB_values)==30*nb, length(R_values)==30*nb)
    stopifnot(all(is.finite(SSB_values)), all(SSB_values>=0))
    stopifnot(all(is.finite(R_values)), all(R_values>=0))

    
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
      REF_map_pairing_PASS=TRUE,
      saved_input_checks_PASS=TRUE,
      min_SSB_Q2=min(SSB_values),
      max_SSB_Q2=max(SSB_values),
      min_R_Q3=min(R_values),
      max_R_Q3=max(R_values),
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
write.csv(block_validation, file.path(output_dir, "21c_refined_new_validation_by_block.csv"), row.names = FALSE)
write.csv(scenario_validation, file.path(output_dir, "21c_refined_new_validation_by_scenario.csv"), row.names = FALSE)
write.csv(branch_summary, file.path(output_dir, "21c_refined_new_HCR_branch_summary.csv"), row.names = FALSE)
write.csv(TAC_summary, file.path(output_dir, "21c_refined_new_TAC_implementation_summary.csv"), row.names = FALSE)
write.csv(TAC_flags, file.path(output_dir, "21c_refined_new_TAC_implementation_flags.csv"), row.names = FALSE)
write.csv(largest_TAC_differences, file.path(output_dir, "21c_refined_new_largest_TAC_Catch_differences.csv"), row.names = FALSE)


validation_object <- list(
  settings=list(
    robustness_id=robustness_id,
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


saveRDS(
  validation_object,
  file.path(output_dir, "21c_refined_new_validation_summary.rds")
)

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
cat("Robustness configuration:",robustness_id,"\n")
cat("REF map pairing: PASS\n")
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


