# ============================================================
# 22c_DIAGNOSE_RECRUITMENT_PROCESS_ERROR.R
# Diagnostic analysis of recruitment process variability
#
# Stock:
#   European anchovy ane.27.9aS (Gulf of Cadiz)
#   Shortcut MSE in FLBEIA + FLasher conditioned from SS3
#
# Role in workflow:
#   22_EVALUATE_CANDIDATE_MPS.R
#     -> 22b_DIAGNOSE_BLIM_METHOD_SENSITIVITY.R
#     -> 22c_DIAGNOSE_RECRUITMENT_PROCESS_ERROR.R
#
# Purpose:
#   Diagnose whether the low probability of SSB < Blim under
#   the Reference OM is associated with limited recruitment
#   process variability and, in particular, with a low frequency
#   of sustained sequences of poor recruitment.
#
# IMPORTANT:
#   - This is a post-processing diagnostic.
#   - No new MSE simulations are performed.
#   - sigmaR is NOT modified.
#   - No alternative recruitment scenario is introduced.
#   - Recruitment is extracted in Q3.
#   - SSB is extracted in Q2.
#
# Reference OM:
#   Common recruitment process-error sigmaR approximately 0.216.
#
# Main questions:
#   1. What recruitment variability is actually realised?
#   2. How frequent are low-recruitment years?
#   3. How frequent are consecutive low-recruitment sequences?
#   4. Are low-SSB events associated with previous poor
#      recruitment?
#   5. Do trajectories reaching Blim differ systematically
#      from trajectories remaining above Blim?
#
# Recruitment diagnostic scale:
#
#   Rrel(i,t) = R(i,t) / median_t[R(i,t)]
#
# where the median is calculated within each source OM.
#
# Diagnostic low-recruitment thresholds:
#   Rrel < 0.75
#   Rrel < 0.50
#   Rrel < 0.25
#
# These thresholds are diagnostic descriptors only.
# They are NOT alternative OM assumptions.
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

scenario_file <- here("data","mse","scenarios","candidate_MP_grid.csv")
scenario_root <- here("data","mse","MP_runs")
diagnostic_dir <- here("data","mse","diagnostics","recruitment_process_error")

dir.create(diagnostic_dir,recursive=TRUE,showWarnings=FALSE)

block_ids <- 1:10
Blim <- 4721
Bpa <- 6561
risk_threshold <- 0.05

low_R_thresholds <- c(
  R075=0.75,
  R050=0.50,
  R025=0.25)

horizons <- list(
  Short=2025:2034,
  Medium=2035:2044,
  Long=2045:2054,
  Full=2025:2054)


# ============================================================
# LOAD CANDIDATE MP GRID
# ============================================================

stopifnot(file.exists(scenario_file))
scenario_grid <- read.csv(scenario_file,stringsAsFactors=FALSE)
required_columns <- c(
  "scenario_number",
  "scenario_id",
  "Fcap",
  "Besc",
  "Blim",
  "Bpa",
  "first_projection_year",
  "last_projection_year")

stopifnot(all(required_columns %in% names(scenario_grid)))

scenario_grid <- scenario_grid %>%arrange(scenario_number)
n_scenarios <- nrow(scenario_grid)

stopifnot(n_distinct(scenario_grid$scenario_number)==n_scenarios)
stopifnot(n_distinct(scenario_grid$scenario_id)==n_scenarios)
stopifnot(all(abs(scenario_grid$Blim-Blim)<1))
stopifnot(all(abs(scenario_grid$Bpa-Bpa)<1))

cat("\n============================================================\n")
cat("RECRUITMENT PROCESS-ERROR DIAGNOSTIC\n")
cat("============================================================\n")

cat("Candidate MPs:",n_scenarios,"\n")
cat("Blim:",Blim,"t | Bpa:",Bpa,"t\n")

# ============================================================
# HELPERS
# ============================================================

add_global_map <- function(dat,map_b) {
  
  out <- dat %>%
    mutate(
      iter_local=as.integer(as.character(iter))
    ) %>%
    left_join(
      map_b %>%
        dplyr::select(
          iter_local,
          iter1000,
          om,
          replicate
        ),
      by="iter_local"
    )
  
  stopifnot(
    !anyNA(out$iter1000)
  )
  
  stopifnot(
    !anyNA(out$om)
  )
  
  stopifnot(
    !anyNA(out$replicate)
  )
  
  out
}


standardise_dynamic <- function(x,map_b,proj_years) {
  
  as.data.frame(x) %>%
    rename(
      value=data
    ) %>%
    mutate(
      year=as.integer(as.character(year))
    ) %>%
    add_global_map(
      map_b
    ) %>%
    filter(
      year %in% proj_years
    ) %>%
    transmute(
      year,
      iter1000,
      om,
      replicate,
      value
    )
}


# ============================================================
# FUNCTION TO CALCULATE RUN LENGTH
# ============================================================

max_run_length <- function(x) {
  
  r <- rle(x)
  
  if(!any(r$values)) {
    return(0L)
  }
  
  max(
    r$lengths[r$values]
  )
}


# ============================================================
# EXTRACT RECRUITMENT AND SSB FOR ONE CANDIDATE MP
# ============================================================

extract_candidate_RS <- function(scenario_row) {
  
  scenario_number <- scenario_row$scenario_number[[1]]
  
  scenario_id <- scenario_row$scenario_id[[1]]
  
  Fcap <- scenario_row$Fcap[[1]]
  
  Besc <- scenario_row$Besc[[1]]
  
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
  
  mp_dir <- file.path(
    scenario_root,
    scenario_id
  )
  
  mp_files <- file.path(
    mp_dir,
    sprintf(
      "block_%02d.rds",
      block_ids
    )
  )
  
  map_files <- file.path(
    mp_dir,
    sprintf(
      "block_%02d_map.csv",
      block_ids
    )
  )
  
  if(!dir.exists(mp_dir)) {
    
    stop(
      "Scenario directory not found: ",
      mp_dir
    )
  }
  
  if(!all(file.exists(mp_files))) {
    
    stop(
      "Missing RDS blocks for scenario ",
      scenario_number
    )
  }
  
  if(!all(file.exists(map_files))) {
    
    stop(
      "Missing map files for scenario ",
      scenario_number
    )
  }
  
  
  # ==========================================================
  # CHECK COMPLETE ITERATION MAP
  # ==========================================================
  
  iteration_maps <- map_dfr(
    block_ids,
    ~read.csv(
      map_files[.x]
    ) %>%
      mutate(
        block=.x
      )
  )
  
  stopifnot(
    nrow(iteration_maps)==1000
  )
  
  stopifnot(
    n_distinct(iteration_maps$iter1000)==1000
  )
  
  stopifnot(
    n_distinct(iteration_maps$om)==100
  )
  
  stopifnot(
    all(table(iteration_maps$om)==10)
  )
  
  stopifnot(
    all(sort(unique(iteration_maps$replicate))==1:10)
  )
  
  
  # ==========================================================
  # EXTRACT ONE BLOCK
  # ==========================================================
  
  extract_block_RS <- function(block_id) {
    
    OM_b <- readRDS(
      mp_files[block_id]
    )
    
    map_b <- read.csv(
      map_files[block_id]
    )
    
    biol <- OM_b$biols$ANE
    
    
    # --------------------------------------------------------
    # Recruitment: age 0, Q3
    # --------------------------------------------------------
    
    R_q3 <- biol@n[
      "0",
      ,
      ,
      "3",
      ,
      
    ]
    
    
    R_b <- standardise_dynamic(
      R_q3,
      map_b,
      proj_years
    ) %>%
      rename(
        R=value
      )
    
    
    # --------------------------------------------------------
    # SSB: Q2
    # --------------------------------------------------------
    
    Mat <- predict(
      biol@mat
    )
    
    SSB_q2 <- quantSums(
      biol@n[,,,"2",,] *
        biol@wt[,,,"2",,] *
        Mat[,,,"2",,]
    )
    
    
    SSB_b <- standardise_dynamic(
      SSB_q2,
      map_b,
      proj_years
    ) %>%
      rename(
        SSB=value
      )
    
    
    # --------------------------------------------------------
    # Join R and SSB
    # --------------------------------------------------------
    
    RS_b <- R_b %>%
      inner_join(
        SSB_b,
        by=c(
          "year",
          "iter1000",
          "om",
          "replicate"
        ),
        relationship="one-to-one"
      )
    
    RS_b
  }
  
  
  # ==========================================================
  # EXTRACT COMPLETE ENSEMBLE
  # ==========================================================
  
  RS_MP <- map_dfr(
    block_ids,
    extract_block_RS
  )
  
  
  # ==========================================================
  # VALIDATION
  # ==========================================================
  
  expected_rows <- length(proj_years)*1000
  
  stopifnot(
    nrow(RS_MP)==expected_rows
  )
  
  stopifnot(
    n_distinct(RS_MP$iter1000)==1000
  )
  
  stopifnot(
    n_distinct(RS_MP$om)==100
  )
  
  stopifnot(
    all(sort(unique(RS_MP$replicate))==1:10)
  )
  
  stopifnot(
    !anyNA(RS_MP$R)
  )
  
  stopifnot(
    !anyNA(RS_MP$SSB)
  )
  
  stopifnot(
    all(is.finite(RS_MP$R))
  )
  
  stopifnot(
    all(is.finite(RS_MP$SSB))
  )
  
  stopifnot(
    all(RS_MP$R>=0)
  )
  
  stopifnot(
    all(RS_MP$SSB>=0)
  )
  
  
  cat(
    "R + SSB extraction passed:",
    nrow(RS_MP),
    "rows |",
    n_distinct(RS_MP$iter1000),
    "trajectories\n"
  )
  
  
  # ==========================================================
  # ADD SCENARIO INFORMATION
  # ==========================================================
  
  RS_MP <- RS_MP %>%
    mutate(
      scenario_number=scenario_number,
      scenario_id=scenario_id,
      Fcap=Fcap,
      Besc=Besc,
      .before=1
    )
  
  RS_MP
}


# ============================================================
# EXTRACT ALL CANDIDATE MPs
# ============================================================

RS_all <- map_dfr(
  seq_len(nrow(scenario_grid)),
  ~extract_candidate_RS(
    scenario_grid[.x,,drop=FALSE]
  )
)


cat("\n============================================================\n")
cat("ALL RECRUITMENT AND SSB TRAJECTORIES EXTRACTED\n")
cat("============================================================\n")

cat(
  "Scenarios:",
  n_distinct(RS_all$scenario_number),
  "\n"
)

cat(
  "Rows:",
  nrow(RS_all),
  "\n"
)


# ============================================================
# STRUCTURAL CHECKS
# ============================================================

stopifnot(
  n_distinct(RS_all$scenario_number)==n_scenarios
)

stopifnot(
  all(
    table(RS_all$scenario_number)==30000
  )
)

stopifnot(
  !anyNA(RS_all$R)
)

stopifnot(
  !anyNA(RS_all$SSB)
)

stopifnot(
  all(RS_all$R>=0)
)

stopifnot(
  all(RS_all$SSB>=0)
)

cat(
  "Structural checks passed\n"
)


# ============================================================
# CHECK WHETHER RECRUITMENT IS IDENTICAL ACROSS MPs
#
# Recruitment process deviations should be paired across
# candidate MPs. This check determines whether the realised
# recruitment is in fact identical after management feedback.
#
# Absolute R may differ because SSB can affect expected
# recruitment through the stock-recruit relationship.
# ============================================================

R_pair_check <- RS_all %>%
  group_by(
    year,
    iter1000
  ) %>%
  summarise(
    n_scenarios=n_distinct(scenario_number),
    R_min=min(R),
    R_max=max(R),
    R_range=R_max-R_min,
    .groups="drop"
  )

cat("\n============================================================\n")
cat("RECRUITMENT ACROSS MP SCENARIOS\n")
cat("============================================================\n")

cat(
  "Maximum absolute R range across MPs:",
  max(R_pair_check$R_range),
  "\n"
)


# ============================================================
# OM-SPECIFIC RECRUITMENT SCALE
#
# Calculate the characteristic recruitment level of each OM.
#
# The scale is calculated separately by scenario because
# realised recruitment can be affected by management through
# the stock-recruit relationship.
# ============================================================

R_OM_scale <- RS_all %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    om
  ) %>%
  summarise(
    R_OM_median=median(R),
    R_OM_mean=mean(R),
    .groups="drop"
  )


stopifnot(
  all(R_OM_scale$R_OM_median>0)
)


# ============================================================
# RELATIVE RECRUITMENT
# ============================================================

RS_relative <- RS_all %>%
  left_join(
    R_OM_scale,
    by=c(
      "scenario_number",
      "scenario_id",
      "Fcap",
      "Besc",
      "om"
    ),
    relationship="many-to-one"
  ) %>%
  mutate(
    Rrel=R/R_OM_median,
    below_Blim=SSB<Blim,
    below_Bpa=SSB<Bpa
  )


stopifnot(
  !anyNA(RS_relative$Rrel)
)

stopifnot(
  all(is.finite(RS_relative$Rrel))
)

stopifnot(
  all(RS_relative$Rrel>=0)
)


# ============================================================
# ADD LOW-RECRUITMENT INDICATORS
# ============================================================

RS_relative <- RS_relative %>%
  mutate(
    Rrel_lt_075=Rrel<low_R_thresholds["R075"],
    Rrel_lt_050=Rrel<low_R_thresholds["R050"],
    Rrel_lt_025=Rrel<low_R_thresholds["R025"]
  )


# ============================================================
# ANNUAL RECRUITMENT DISTRIBUTION
# ============================================================

R_annual <- RS_relative %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    year
  ) %>%
  summarise(
    n=dplyr::n(),
    R_min=min(R),
    R_q01=quantile(R,0.01),
    R_q05=quantile(R,0.05),
    R_q10=quantile(R,0.10),
    R_median=median(R),
    R_mean=mean(R),
    R_sd=sd(R),
    R_CV=sd(R)/mean(R),
    Rrel_min=min(Rrel),
    Rrel_q01=quantile(Rrel,0.01),
    Rrel_q05=quantile(Rrel,0.05),
    Rrel_q10=quantile(Rrel,0.10),
    Rrel_median=median(Rrel),
    P_Rrel_lt_075=mean(Rrel_lt_075),
    P_Rrel_lt_050=mean(Rrel_lt_050),
    P_Rrel_lt_025=mean(Rrel_lt_025),
    P_Blim=mean(below_Blim),
    .groups="drop"
  )


stopifnot(
  all(R_annual$n==1000)
)


# ============================================================
# TRAJECTORY-LEVEL RECRUITMENT DIAGNOSTICS
# ============================================================

R_trajectory <- RS_relative %>%
  arrange(
    scenario_number,
    iter1000,
    year
  ) %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    iter1000,
    om,
    replicate
  ) %>%
  summarise(
    R_min=min(R),
    R_mean=mean(R),
    R_median=median(R),
    R_CV=sd(R)/mean(R),
    Rrel_min=min(Rrel),
    Rrel_mean=mean(Rrel),
    Rrel_q05=quantile(Rrel,0.05),
    n_Rrel_lt_075=sum(Rrel_lt_075),
    n_Rrel_lt_050=sum(Rrel_lt_050),
    n_Rrel_lt_025=sum(Rrel_lt_025),
    max_run_lt_075=max_run_length(Rrel_lt_075),
    max_run_lt_050=max_run_length(Rrel_lt_050),
    max_run_lt_025=max_run_length(Rrel_lt_025),
    min_SSB=min(SSB),
    ever_Blim=any(below_Blim),
    ever_Bpa=any(below_Bpa),
    .groups="drop"
  )


# ============================================================
# FREQUENCY OF LOW-RECRUITMENT SEQUENCES
# ============================================================

R_sequence_summary <- R_trajectory %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc
  ) %>%
  summarise(
    P_any_R075=mean(n_Rrel_lt_075>0),
    P_any_R050=mean(n_Rrel_lt_050>0),
    P_any_R025=mean(n_Rrel_lt_025>0),
    P_run2_R075=mean(max_run_lt_075>=2),
    P_run3_R075=mean(max_run_lt_075>=3),
    P_run2_R050=mean(max_run_lt_050>=2),
    P_run3_R050=mean(max_run_lt_050>=3),
    P_run2_R025=mean(max_run_lt_025>=2),
    P_run3_R025=mean(max_run_lt_025>=3),
    median_max_run_R075=median(max_run_lt_075),
    max_run_R075=max(max_run_lt_075),
    median_max_run_R050=median(max_run_lt_050),
    max_run_R050=max(max_run_lt_050),
    median_max_run_R025=median(max_run_lt_025),
    max_run_R025=max(max_run_lt_025),
    .groups="drop"
  )


# ============================================================
# DEFINE EVALUATION HORIZONS
# ============================================================

horizon_table <- imap_dfr(
  horizons,
  ~tibble(
    horizon=.y,
    year=.x
  )
)


RS_horizon <- RS_relative %>%
  inner_join(
    horizon_table,
    by="year",
    relationship="many-to-many"
  )


# ============================================================
# RECRUITMENT DISTRIBUTION BY HORIZON
# ============================================================

R_horizon <- RS_horizon %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc,
    horizon
  ) %>%
  summarise(
    R_min=min(R),
    R_q01=quantile(R,0.01),
    R_q05=quantile(R,0.05),
    R_q10=quantile(R,0.10),
    R_median=median(R),
    R_mean=mean(R),
    R_CV=sd(R)/mean(R),
    Rrel_min=min(Rrel),
    Rrel_q01=quantile(Rrel,0.01),
    Rrel_q05=quantile(Rrel,0.05),
    Rrel_q10=quantile(Rrel,0.10),
    Rrel_median=median(Rrel),
    P_Rrel_lt_075=mean(Rrel_lt_075),
    P_Rrel_lt_050=mean(Rrel_lt_050),
    P_Rrel_lt_025=mean(Rrel_lt_025),
    P_Blim=mean(below_Blim),
    .groups="drop"
  )


# ============================================================
# LAGGED RECRUITMENT
#
# SSB in year t occurs in Q2 and recruitment age 0 enters in Q3.
# Therefore R(t) cannot explain SSB(t).
#
# Recruitment is lagged so that SSB(t) is compared with
# recruitment realised in previous years.
# ============================================================

RS_lagged <- RS_relative %>%
  arrange(
    scenario_number,
    iter1000,
    year
  ) %>%
  group_by(
    scenario_number,
    iter1000
  ) %>%
  mutate(
    Rrel_lag1=lag(Rrel,1),
    Rrel_lag2=lag(Rrel,2),
    Rrel_lag3=lag(Rrel,3),
    R_lag1=lag(R,1),
    R_lag2=lag(R,2),
    R_lag3=lag(R,3)
  ) %>%
  ungroup()


# ============================================================
# RECRUITMENT PRECEDING SSB < BLIM
# ============================================================

R_before_Blim <- RS_lagged %>%
  filter(
    below_Blim
  ) %>%
  summarise(
    n_Blim_events=dplyr::n(),
    n_trajectories=n_distinct(iter1000),
    median_Rrel_lag1=median(Rrel_lag1,na.rm=TRUE),
    q05_Rrel_lag1=quantile(Rrel_lag1,0.05,na.rm=TRUE),
    median_Rrel_lag2=median(Rrel_lag2,na.rm=TRUE),
    q05_Rrel_lag2=quantile(Rrel_lag2,0.05,na.rm=TRUE),
    median_Rrel_lag3=median(Rrel_lag3,na.rm=TRUE),
    q05_Rrel_lag3=quantile(Rrel_lag3,0.05,na.rm=TRUE)
  )


# ============================================================
# EVENT-LEVEL COMPARISON:
# SSB BELOW VS ABOVE BLIM
# ============================================================

R_Blim_event_comparison <- RS_lagged %>%
  filter(
    !is.na(Rrel_lag1)
  ) %>%
  mutate(
    SSB_state=if_else(
      below_Blim,
      "Below_Blim",
      "Above_Blim"
    )
  ) %>%
  group_by(
    SSB_state
  ) %>%
  summarise(
    n=dplyr::n(),
    Rrel_lag1_median=median(Rrel_lag1,na.rm=TRUE),
    Rrel_lag1_q05=quantile(Rrel_lag1,0.05,na.rm=TRUE),
    Rrel_lag2_median=median(Rrel_lag2,na.rm=TRUE),
    Rrel_lag2_q05=quantile(Rrel_lag2,0.05,na.rm=TRUE),
    Rrel_lag3_median=median(Rrel_lag3,na.rm=TRUE),
    Rrel_lag3_q05=quantile(Rrel_lag3,0.05,na.rm=TRUE),
    .groups="drop"
  )


# ============================================================
# TRAJECTORY COMPARISON:
# EVER BELOW BLIM VS NEVER BELOW BLIM
# ============================================================

R_trajectory_comparison <- R_trajectory %>%
  mutate(
    trajectory_state=if_else(
      ever_Blim,
      "Ever_Blim",
      "Never_Blim"
    )
  ) %>%
  group_by(
    trajectory_state
  ) %>%
  summarise(
    n_trajectories=dplyr::n(),
    median_min_Rrel=median(Rrel_min),
    q05_min_Rrel=quantile(Rrel_min,0.05),
    median_n_R075=median(n_Rrel_lt_075),
    median_n_R050=median(n_Rrel_lt_050),
    median_n_R025=median(n_Rrel_lt_025),
    median_max_run_R075=median(max_run_lt_075),
    median_max_run_R050=median(max_run_lt_050),
    median_max_run_R025=median(max_run_lt_025),
    median_min_SSB=median(min_SSB),
    .groups="drop"
  )


# ============================================================
# TRAJECTORY-LEVEL ASSOCIATION BETWEEN RECRUITMENT AND MIN SSB
# ============================================================

R_SSB_association <- R_trajectory %>%
  group_by(
    scenario_number,
    scenario_id,
    Fcap,
    Besc
  ) %>%
  summarise(
    cor_minRrel_minSSB=cor(
      Rrel_min,
      min_SSB,
      method="spearman"
    ),
    cor_nR075_minSSB=cor(
      n_Rrel_lt_075,
      min_SSB,
      method="spearman"
    ),
    cor_nR050_minSSB=cor(
      n_Rrel_lt_050,
      min_SSB,
      method="spearman"
    ),
    cor_runR075_minSSB=cor(
      max_run_lt_075,
      min_SSB,
      method="spearman"
    ),
    cor_runR050_minSSB=cor(
      max_run_lt_050,
      min_SSB,
      method="spearman"
    ),
    .groups="drop"
  )


# ============================================================
# GLOBAL RECRUITMENT DIAGNOSTIC
# ============================================================

global_R_summary <- RS_relative %>%
  summarise(
    n=dplyr::n(),
    R_min=min(R),
    R_q01=quantile(R,0.01),
    R_q05=quantile(R,0.05),
    R_q10=quantile(R,0.10),
    R_median=median(R),
    R_mean=mean(R),
    R_CV=sd(R)/mean(R),
    Rrel_min=min(Rrel),
    Rrel_q01=quantile(Rrel,0.01),
    Rrel_q05=quantile(Rrel,0.05),
    Rrel_q10=quantile(Rrel,0.10),
    Rrel_median=median(Rrel),
    P_Rrel_lt_075=mean(Rrel_lt_075),
    P_Rrel_lt_050=mean(Rrel_lt_050),
    P_Rrel_lt_025=mean(Rrel_lt_025),
    P_Blim=mean(below_Blim)
  )


# ============================================================
# FULL-HORIZON SUMMARY BY MP
# ============================================================

R_full <- R_horizon %>%
  filter(
    horizon=="Full"
  ) %>%
  arrange(
    Fcap,
    Besc
  )


# ============================================================
# PRINT MAIN DIAGNOSTICS
# ============================================================

cat("\n============================================================\n")
cat("GLOBAL REALISED RECRUITMENT DISTRIBUTION\n")
cat("============================================================\n")

print(global_R_summary)


cat("\n============================================================\n")
cat("FULL-HORIZON RECRUITMENT DIAGNOSTIC BY MP\n")
cat("============================================================\n")

print(
  R_full %>%
    dplyr::select(
      scenario_number,
      Fcap,
      Besc,
      Rrel_q01,
      Rrel_q05,
      Rrel_q10,
      P_Rrel_lt_075,
      P_Rrel_lt_050,
      P_Rrel_lt_025,
      P_Blim
    ),
  n=Inf
)


cat("\n============================================================\n")
cat("LOW-RECRUITMENT SEQUENCES\n")
cat("============================================================\n")

print(
  R_sequence_summary %>%
    dplyr::select(
      scenario_number,
      Fcap,
      Besc,
      P_run2_R075,
      P_run3_R075,
      P_run2_R050,
      P_run3_R050,
      P_run2_R025,
      P_run3_R025
    ),
  n=Inf
)


cat("\n============================================================\n")
cat("RECRUITMENT BEFORE SSB < BLIM\n")
cat("============================================================\n")

print(
  R_before_Blim,
  width=Inf
)


cat("\n============================================================\n")
cat("TRAJECTORIES WITH VS WITHOUT BLIM EVENTS\n")
cat("============================================================\n")

print(
  R_trajectory_comparison,
  width=Inf
)


cat("\n============================================================\n")
cat("RECRUITMENT - MINIMUM SSB ASSOCIATION\n")
cat("============================================================\n")

print(
  R_SSB_association,
  n=Inf
)


# ============================================================
# FINAL VALIDATION
# ============================================================

stopifnot(
  n_distinct(R_full$scenario_number)==n_scenarios
)

stopifnot(
  all(
    R_full$P_Rrel_lt_075>=0 &
      R_full$P_Rrel_lt_075<=1
  )
)

stopifnot(
  all(
    R_full$P_Rrel_lt_050>=0 &
      R_full$P_Rrel_lt_050<=1
  )
)

stopifnot(
  all(
    R_full$P_Rrel_lt_025>=0 &
      R_full$P_Rrel_lt_025<=1
  )
)

stopifnot(
  all(
    R_sequence_summary$P_run2_R075>=0 &
      R_sequence_summary$P_run2_R075<=1
  )
)

stopifnot(
  all(
    R_sequence_summary$P_run3_R075>=0 &
      R_sequence_summary$P_run3_R075<=1
  )
)

cat("\nFINAL VALIDATION PASSED\n")


# ============================================================
# SAVE OUTPUTS
# ============================================================

diagnostic_outputs <- list(
  recruitment_OM_scale=R_OM_scale,
  annual_recruitment=R_annual,
  recruitment_by_horizon=R_horizon,
  recruitment_trajectory=R_trajectory,
  sequence_summary=R_sequence_summary,
  recruitment_before_Blim=R_before_Blim,
  Blim_event_comparison=R_Blim_event_comparison,
  trajectory_comparison=R_trajectory_comparison,
  R_SSB_association=R_SSB_association,
  global_R_summary=global_R_summary,
  recruitment_pair_check=R_pair_check
)


saveRDS(
  diagnostic_outputs,
  file.path(
    diagnostic_dir,
    "recruitment_process_error_diagnostic.rds"
  )
)


write.csv(
  R_OM_scale,
  file.path(
    diagnostic_dir,
    "recruitment_OM_scale.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_annual,
  file.path(
    diagnostic_dir,
    "recruitment_annual_distribution.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_horizon,
  file.path(
    diagnostic_dir,
    "recruitment_by_horizon.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_trajectory,
  file.path(
    diagnostic_dir,
    "recruitment_trajectory_diagnostics.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_sequence_summary,
  file.path(
    diagnostic_dir,
    "recruitment_low_sequence_summary.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_before_Blim,
  file.path(
    diagnostic_dir,
    "recruitment_before_Blim.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_Blim_event_comparison,
  file.path(
    diagnostic_dir,
    "recruitment_Blim_event_comparison.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_trajectory_comparison,
  file.path(
    diagnostic_dir,
    "recruitment_trajectory_Blim_comparison.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_SSB_association,
  file.path(
    diagnostic_dir,
    "recruitment_SSB_association.csv"
  ),
  row.names=FALSE
)


write.csv(
  global_R_summary,
  file.path(
    diagnostic_dir,
    "recruitment_global_summary.csv"
  ),
  row.names=FALSE
)


write.csv(
  R_pair_check,
  file.path(
    diagnostic_dir,
    "recruitment_pair_check.csv"
  ),
  row.names=FALSE
)


# ============================================================
# END
# ============================================================

cat("\n============================================================\n")
cat("RECRUITMENT PROCESS-ERROR DIAGNOSTIC COMPLETED\n")
cat("============================================================\n")

cat("Candidate MPs evaluated:",n_scenarios,"\n")
cat("Recruitment season: Q3\n")
cat("SSB season: Q2\n")
cat("No alternative sigmaR was simulated\n")
cat("Outputs saved in:",diagnostic_dir,"\n")