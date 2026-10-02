# ============================================================
# 12_RUN_OM_REFERENCE_VALIDATION.R
# Open-loop validation of the reference operating model
#
# Purpose:
#   Run open-loop projections of the conditioned operating model
#   to validate biological and fishery dynamics before applying
#   a management procedure.
#
# Validation scenarios:
#   C0    - No fishing
#   C7000 - Constant annual catch advice of 7,000 t
#
# INPUTS:
#   data/Rdata/
#     - iteration_map_1000.csv
#     - biols_conditioned_1000iter.rds
#     - SRs_reference_1000iter_process_error.rds
#     - fleets_conditioned_1000iter.rds
#     - fleets_ctrl_reference_1000iter.rds
#
# OUTPUTS:
#   data/mse/OM_validation_C0/
#     - OM_validation_C0_block_XX.rds
#     - OM_validation_C0_block_XX_map.csv
#     - OM_validation_C0_block_XX_metadata.csv
#
#   data/mse/OM_validation_C7000/
#     - OM_validation_C7000_block_XX.rds
#     - OM_validation_C7000_block_XX_map.csv
#     - OM_validation_C7000_block_XX_metadata.csv
#
# Simulation period:
#   Historical: 1989-2024
#   Projection: 2025-2054
#
# Simulation structure:
#   Conditioned OMs:      100
#   Replicates per OM:     10
#   Total trajectories:  1000
#   Block size:            100
#   Blocks per scenario:    10
#
# Configuration:
#   - Perfect observation
#   - No assessment
#   - Fixed catch advice
#
# Usage:
#   Rscript RUN_OM_REFERENCE_VALIDATION.R C0 <block 1-10>
#   Rscript RUN_OM_REFERENCE_VALIDATION.R C7000 <block 1-10>
# ============================================================
rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)

# ============================================================
#'*DIRECTORIES*
# ============================================================

rds_boot_dir <- file.path(here(),"data","Rdata")

input_files <- file.path(rds_boot_dir,c("iteration_map_1000.csv",
                                        "biols_conditioned_1000iter.rds",
                                        "SRs_reference_1000iter_process_error.rds",
                                        "fleets_conditioned_1000iter.rds",
                                        "fleets_ctrl_reference_1000iter.rds"))

if(!all(file.exists(input_files))) stop("Missing input files:\n",paste(input_files[!file.exists(input_files)],collapse="\n"))

# ============================================================
#'*PARAMETERS*
# ============================================================

first.yr <- 1989
last.obs.yr <- 2024
proj.yr <- 2025
last.yr <- 2054

ni <- 1000
ns <- 4
block_size <- 100

stks <- "ANE"

Blim <- 4721
Bpa <- 6561

# ============================================================
#'*COMMAND-LINE ARGUMENTS*
# ============================================================

args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=2) stop("Usage: Rscript RUN_OM_REFERENCE.R <C0|C7000> <block 1-10>")

scenario_id <- args[1]
block <- as.integer(args[2])

if(!scenario_id %in% c("C0","C7000")) stop("Scenario must be C0 or C7000")

validation_TAC <- switch(scenario_id,C0=0,C7000=7000)
scenario <- switch(scenario_id,C0="no_fishing",C7000="constant_catch")

n_blocks <- ceiling(ni/block_size)

if(is.na(block) || block<1 || block>n_blocks) stop("Invalid block number. Must be between 1 and ",n_blocks)

out_dir <- file.path(here(),"data","mse",paste0("OM_validation_",scenario_id))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

cat("\n============================================\n")
cat("OPEN-LOOP OM VALIDATION\n")
cat("Scenario:",scenario_id,"-",scenario,"\n")
cat("Catch advice:",validation_TAC,"t\n")
cat("Block:",block,"/",n_blocks,"\n")
cat("============================================\n")

# ============================================================
#'*LOAD ITERATION MAP*
# ============================================================
iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map_1000.csv"))

stopifnot(all(c("iter1000","om","replicate") %in% names(iteration_map)),
          nrow(iteration_map)==ni,
          n_distinct(iteration_map$iter1000)==ni,
          n_distinct(iteration_map$om)==100)

cat("\nTotal trajectories:",nrow(iteration_map),"\n")
cat("Number of conditioned OMs:",n_distinct(iteration_map$om),"\n")

# ============================================================
#'*SELECT BLOCK*
# ============================================================

first_iter <- (block-1)*block_size+1
last_iter <- min(block*block_size,ni)

ii <- first_iter:last_iter
nb <- length(ii)

cat("\n============================================\n",
    "RUNNING BLOCK",block,"OF",n_blocks,"\n",
    "Global iterations:",first_iter,"-",last_iter,"\n",
    "Number of trajectories:",nb,"\n",
    "============================================\n")


# ============================================================
#'*RUN MAP*
#
# Maps local FLBEIA iteration 1:nb to:
#   - original iter1000
#   - conditioned OM
#   - stochastic replicate
# ============================================================

run_map <- iteration_map %>%
  filter(iter1000 %in% ii) %>%
  arrange(match(iter1000,ii)) %>%
  mutate(iter_local=seq_len(n()),scenario_id=scenario_id,scenario=scenario)

stopifnot(nrow(run_map)==nb,
          all(run_map$iter1000==ii),
          n_distinct(run_map$iter_local)==nb)

cat("\nRun map:\n")
print(run_map)

# ============================================================
#'*LOAD CONDITIONED OM*
# ============================================================
biols <- readRDS(file.path(rds_boot_dir,"biols_conditioned_1000iter.rds"))
SRs <- readRDS(file.path(rds_boot_dir,"SRs_reference_1000iter_process_error.rds"))
fleets <- readRDS(file.path(rds_boot_dir,"fleets_conditioned_1000iter.rds"))
fleets.ctrl <- readRDS(file.path(rds_boot_dir,"fleets_ctrl_reference_1000iter.rds"))

stopifnot(
  "ANE" %in% names(biols),
  "ANE" %in% names(SRs),
  "SEINE" %in% names(fleets),
  dim(biols$ANE@n)[6]==ni,
  dim(SRs$ANE@params)[4]==ni,
  dim(SRs$ANE@uncertainty)[6]==ni,
  dim(fleets$SEINE@effort)[6]==ni,
  dim(catch.q(fleets$SEINE@metiers$ALL@catches$ANE))[6]==ni
)

cat("Conditioned OM objects loaded: OK\n")

# ============================================================
#'*SUBSETTING FUNCTIONS*
# ============================================================

sub_flq <- function(x,ii) {
  jj <- if(dim(x)[6]==1) rep(1,length(ii)) else ii
  out <- x[,,,,,jj,drop=FALSE]
  dimnames(out)$iter <- as.character(seq_along(ii))
  out
}

sub_FLBiol <- function(x,ii) {
  out <- x
  
  for(s in slotNames(x)) {
    z <- slot(x,s)
    
    if(inherits(z,"FLQuant")) slot(out,s,check=FALSE) <- sub_flq(z,ii)
    
    if(inherits(z,"predictModel")) {
      z_new <- z
      z_new@.Data[[1]] <- sub_flq(z@.Data[[1]],ii)
      slot(out,s,check=FALSE) <- z_new
    }
  }
  
  validObject(out)
  out
}

sub_FLSR <- function(x,ii) {
  out <- x
  
  for(s in c("rec","ssb","uncertainty","proportion")) {
    slot(out,s,check=FALSE) <- sub_flq(slot(x,s),ii)
  }
  
  pars <- x@params[,,,ii,drop=FALSE]
  dimnames(pars)$iter <- as.character(seq_along(ii))
  slot(out,"params",check=FALSE) <- pars
  
  validObject(out)
  out
}

sub_fleet <- function(x,ii) {
  out <- x
  
  # Fleet-level FLQuant slots
  for(s in slotNames(x)) {
    z <- slot(x,s)
    if(inherits(z,"FLQuant")) slot(out,s,check=FALSE) <- sub_flq(z,ii)
  }
  
  # Metier-level FLQuant slots
  met <- x@metiers$ALL
  
  for(s in slotNames(met)) {
    z <- slot(met,s)
    if(inherits(z,"FLQuant")) slot(met,s,check=FALSE) <- sub_flq(z,ii)
  }
  
  # Catch-level FLQuant slots
  ca <- met@catches$ANE
  
  for(s in slotNames(ca)) {
    z <- slot(ca,s)
    if(inherits(z,"FLQuant")) slot(ca,s,check=FALSE) <- sub_flq(z,ii)
  }
  
  catches_b <- met@catches
  catches_b$ANE <- ca
  slot(met,"catches",check=FALSE) <- catches_b
  
  metiers_b <- out@metiers
  metiers_b$ALL <- met
  slot(out,"metiers",check=FALSE) <- metiers_b
  
  validObject(out)
  out
}
# ============================================================
#'*CONDITIONED OBJECTS FOR CURRENT RUN*
# ============================================================

biols_b <- biols
biols_b$ANE <- sub_FLBiol(biols$ANE,ii)

SRs_b <- SRs
SRs_b$ANE <- sub_FLSR(SRs$ANE,ii)

fleets_b <- fleets
fleets_b$SEINE <- sub_fleet(fleets$SEINE,ii)

fleets.ctrl_b <- fleets.ctrl
fleets.ctrl_b$seasonal.share[[1]] <- sub_flq(fleets.ctrl$seasonal.share[[1]],ii)

stopifnot(
  dim(biols_b$ANE@n)[6]==nb,
  dim(SRs_b$ANE@params)[4]==nb,
  dim(SRs_b$ANE@uncertainty)[6]==nb,
  dim(fleets_b$SEINE@effort)[6]==nb,
  dim(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE))[6]==nb,
  dim(fleets.ctrl_b$seasonal.share[[1]])[6]==nb
)

cat("Conditioned objects subset to current block: OK\n")

#'*============================================================*
#'*CONTROLS*
#'*============================================================*

biols.ctrl <- create.biols.ctrl(stksnames=stks,growth.model="ASPG_Baranov")
main.ctrl <- list(sim.years=c(initial=last.obs.yr,final=last.yr))

covars <- NULL
covars.ctrl <- NULL
indices <- NULL

#'*============================================================*
#'*PERFECT OBSERVATION*
#
# Open-loop validation: no observation error.
#'*============================================================*

flq.ANE <- FLQuant(dimnames=list(age="all",year=first.yr:last.yr,unit=1,season=1:ns,iter=1:nb))

obs.ctrl <- create.obs.ctrl(stksnames=stks,flq.ANE=flq.ANE,stkObs.models="perfectObs")
obs.ctrl$ANE$obs.curryr <- TRUE

#'*============================================================*
#'*NO ASSESSMENT*
#
# Open-loop validation: no assessment model or assessment error.
#'*============================================================*

assess.ctrl <- create.assess.ctrl(stksnames=stks,assess.models="NoAssessment")
assess.ctrl$ANE$ass.curryr <- TRUE

#'*============================================================*
#'*ADVICE DATA*
#'*============================================================*

ANE_advice.TAC.flq <- FLQuant(c(seasonSums(quantSums(catchWStock(fleets_b,stock="ANE")))),dimnames=list(quant=stks,year=first.yr:last.yr),iter=nb)

# ------------------------------------------------------------
# Historical TACs
# ------------------------------------------------------------

ANE_advice.TAC.flq[,ac(1989:2018),] <- NA
ANE_advice.TAC.flq[,ac(2019:2024),] <- c(5278,8856,9459,4383,1892,1733)

# ------------------------------------------------------------
# Fixed catch advice during projection
#
# Open-loop validation scenario, not a management procedure.
# ------------------------------------------------------------

ANE_advice.TAC.flq[,ac(proj.yr:last.yr),] <- validation_TAC

# ------------------------------------------------------------
# Quota share
# ------------------------------------------------------------

ANE_advice.quota.share.flq <- FLQuant(1,dimnames=list(quant=stks,year=first.yr:last.yr),iter=nb)
ANE_advice.avg.yrs <- c(2023,2024)

stksTac.data <- list(ANE=ls(pattern="^ANE"))
yrs <- c(first.yr=first.yr,proj.yr=proj.yr,last.yr=last.yr)

advice_fixTAC <- create.advice.data(yrs,ns,nb,stksTac.data,fleets_b)
units(advice_fixTAC$TAC) <- "t"

stopifnot(all(as.numeric(advice_fixTAC$TAC[,ac(proj.yr:last.yr),,,,])==validation_TAC))

cat("Projection catch advice:",validation_TAC,"t: OK\n")

#'*============================================================*
#'*ADVICE CONTROL*
#'*============================================================*

advice_fixTAC.ctrl <- list(ANE=list(HCR.model="fixedAdvice"))

#'*============================================================*
#'*PRE-RUN CHECKS*
#'*============================================================*

cat("\n============================================\n",
    "PRE-RUN CHECKS\n",
    "============================================\n")

stopifnot(
  dim(biols_b$ANE@n)[6]==nb,
  dim(SRs_b$ANE@params)[4]==nb,
  dim(SRs_b$ANE@uncertainty)[6]==nb,
  dim(fleets_b$SEINE@effort)[6]==nb,
  dim(fleets_b$SEINE@metiers$ALL@vcost)[6]==nb,
  dim(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE))[6]==nb,
  dim(fleets.ctrl_b$seasonal.share[[1]])[6]==nb,
  all(is.finite(as.numeric(fleets_b$SEINE@metiers$ALL@vcost))),
  all(is.finite(as.numeric(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE))))
)

cat("Block inputs: OK\n")
# ============================================================
# CHECK FUTURE RECRUITMENT PROCESS ERROR
# ============================================================

future_unc <- as.numeric(SRs_b$ANE@uncertainty[,ac(proj.yr:last.yr),,,,])

stopifnot(length(future_unc)>0,
          all(is.finite(future_unc)),
          all(future_unc>0))

cat("Recruitment multiplier mean =",mean(future_unc),"\n")
cat("Recruitment multiplier SD =",sd(future_unc),"\n")
cat("Recruitment multiplier range =",range(future_unc),"\n")

# ============================================================
# RUN FLBEIA
# ============================================================

cat("\n============================================\n",
    "RUNNING OPEN-LOOP OM VALIDATION\n",
    "Scenario:",scenario_id,"-",scenario,"\n",
    "Catch advice:",validation_TAC,"t\n",
    "Block:",block,"/",n_blocks,"\n",
    "============================================\n")

start_time <- Sys.time()

OM_b <- FLBEIA(
  biols=biols_b,
  fleets=fleets_b,
  SRs=SRs_b,
  BDs=NULL,
  covars=covars,
  indices=indices,
  advice=advice_fixTAC,
  main.ctrl=main.ctrl,
  biols.ctrl=biols.ctrl,
  fleets.ctrl=fleets.ctrl_b,
  covars.ctrl=covars.ctrl,
  obs.ctrl=obs.ctrl,
  assess.ctrl=assess.ctrl,
  advice.ctrl=advice_fixTAC.ctrl
)

end_time <- Sys.time()

cat("\nRuntime:",round(as.numeric(difftime(end_time,start_time,units="mins")),2),"minutes\n")
# ============================================================
# POST-RUN CHECKS
# ============================================================

stopifnot(dim(OM_b$biols$ANE@n)[6]==nb)

future_n <- as.numeric(OM_b$biols$ANE@n[,ac(proj.yr:last.yr),,,,])
stopifnot(length(future_n)>0,all(is.finite(future_n)),all(future_n>=0))

cat("Projected abundance: OK\n")

# ============================================================
# CHECK REALISED RECRUITMENT
#
# Recruitment = age 0, Q3
# ============================================================

future_R <- as.numeric(OM_b$biols$ANE@n["0",ac(proj.yr:last.yr),,"3",,])
stopifnot(length(future_R)>0,all(is.finite(future_R)),all(future_R>=0))

cat("Projected recruitment: OK\n")
cat("Median projected recruitment =",median(future_R),"\n")

# ============================================================
# CHECK VALIDATION SCENARIO
# ============================================================

cobj_out <- OM_b$fleets$SEINE@metiers$ALL@catches$ANE
future_landings <- as.numeric(cobj_out@landings[,ac(proj.yr:last.yr),,,,])

stopifnot(length(future_landings)>0,all(is.finite(future_landings)),all(future_landings>=0))

if(scenario_id=="C0") {
  stopifnot(max(abs(future_landings))<1e-8)
  cat("No-fishing scenario: realised landings = 0: OK\n")
}

if(scenario_id=="C7000") {
  TAC_proj <- as.numeric(OM_b$advice$TAC[,ac(proj.yr:last.yr),,,,])
  stopifnot(length(TAC_proj)>0,all(is.finite(TAC_proj)),all(TAC_proj==validation_TAC))
  cat("Constant-catch scenario: catch advice =",validation_TAC,"t: OK\n")
}

# ============================================================
# SAVE OUTPUT
# ============================================================

OM_file <- file.path(out_dir,sprintf("OM_validation_%s_block_%02d.rds",scenario_id,block))
map_out_file <- file.path(out_dir,sprintf("OM_validation_%s_block_%02d_map.csv",scenario_id,block))
metadata_file <- file.path(out_dir,sprintf("OM_validation_%s_block_%02d_metadata.csv",scenario_id,block))

saveRDS(OM_b,OM_file)
write.csv(run_map,map_out_file,row.names=FALSE)

# ============================================================
# SAVE RUN METADATA
# ============================================================

run_metadata <- data.frame(
  scenario=scenario,
  scenario_id=scenario_id,
  first_year=first.yr,
  last_historical_year=last.obs.yr,
  first_projection_year=proj.yr,
  last_projection_year=last.yr,
  n_trajectories=nb,
  n_conditioned_OMs=n_distinct(run_map$om),
  validation_TAC=validation_TAC,
  Blim=Blim,
  Bpa=Bpa,
  runtime_minutes=as.numeric(difftime(end_time,start_time,units="mins"))
)

write.csv(run_metadata,metadata_file,row.names=FALSE)

# ============================================================
# COMPLETION MESSAGE
# ============================================================

cat("\n============================================\n")
cat("OPEN-LOOP OM VALIDATION BLOCK COMPLETED\n")
cat("============================================\n")
cat("Scenario:",scenario_id,"-",scenario,"\n")
cat("Block:",block,"/",n_blocks,"\n")
cat("Global iterations:",first_iter,"-",last_iter,"\n")
cat("Historical:",first.yr,"-",last.obs.yr,"\n")
cat("Projection:",proj.yr,"-",last.yr,"\n")
cat("Conditioned OMs:",n_distinct(run_map$om),"\n")
cat("Trajectories:",nb,"\n")
cat("Catch advice:",validation_TAC,"t\n")
cat("Runtime:",round(as.numeric(difftime(end_time,start_time,units="mins")),2),"minutes\n")

cat("\nSaved OM:\n",OM_file,"\n")
cat("\nSaved iteration map:\n",map_out_file,"\n")
cat("\nSaved metadata:\n",metadata_file,"\n")

