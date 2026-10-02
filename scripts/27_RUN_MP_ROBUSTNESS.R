# ============================================================
# 27_RUN_MP_ROBUSTNESS.R
#
# Generic closed-loop runner for Fcap + Bescapement MSE
#
# Historical period : 1989-2024
# Projection period : 2025-2054
#
# Configuration:
#   - Perfect observation
#   - Assessment error in SSB and F controlled by scenario table
#   - No assessment error in recruitment
#   - Probabilistic Fcap + Bescapement HCR
#   - Empirical forecast uncertainty (S2) paired by global trajectory
#
# Usage:
#   Rscript scripts/27_RUN_MP_ROBUSTNESS.R <robustness_id> <reference_scenario_number> <block> [prepare|test|full]
# ============================================================

rm(list=ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)

# ============================================================
# DIRECTORIES
# ============================================================

rds_boot_dir <- file.path(here(),"data","Rdata")
scenario_dir <- file.path(here(),"data","mse","scenarios")
default_scenario_file <- file.path(scenario_dir,"candidate_MP_grid.csv")

input_files <- file.path(rds_boot_dir,c("iteration_map_1000.csv",
                                        "biols_conditioned_1000iter.rds",
                                        "fleets_conditioned_1000iter.rds",
                                        "fleets_ctrl_reference_1000iter.rds",
                                        "ane_stock_conditioned_1000iter.rds"))
if(!all(file.exists(input_files))) stop("Missing input files:\n",
                                        paste(input_files[!file.exists(input_files)],collapse="\n"))
# ============================================================
# PARAMETERS
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
nFyears <- 5

# ============================================================
# LOAD CUSTOM FUNCTIONS
# ============================================================

source(file.path(here(),"functions","perfectObs4seas.R"))
source(file.path(here(),"functions","FcapBpaHCR_ane9aS_PROBABILITY_OPTION.R"))

environment(FcapBpaHCR_ane) <- asNamespace("FLBEIA")
assignInNamespace("annualTAC",FcapBpaHCR_ane,ns="FLBEIA")

txt <- strsplit(paste(deparse(body(getFromNamespace("annualTAC","FLBEIA"))),collapse="\n"),"\n")[[1]]
F_lines <- grep("Fage.*Fseason",txt,value=TRUE)

cat("\nLoaded F equations:\n")
print(F_lines)
stopifnot(length(F_lines)==4,!any(grepl("dt",F_lines)))

# ============================================================
# ROBUSTNESS DESIGN AND EXECUTION MODE
# CLI: <robustness_id> <reference_scenario_number> <block> [prepare|test|full]
# Interactive defaults: S-C1, reference MP 10, block 1, prepare.
# Override with MSE_ROBUSTNESS_ID, MSE_REFERENCE_MP, MSE_BLOCK, MSE_MODE.
# prepare validates/saves block inputs without running FLBEIA.
# test runs 10 trajectories; full runs the 100-trajectory block.
# ============================================================

if(interactive()) {
  robustness_id <- Sys.getenv("MSE_ROBUSTNESS_ID", "S-C1")
  scenario_number <- as.integer(Sys.getenv("MSE_REFERENCE_MP", "10"))
  block <- as.integer(Sys.getenv("MSE_BLOCK", "1"))
  execution_mode <- Sys.getenv("MSE_MODE", "prepare")
} else {
  args <- commandArgs(trailingOnly=TRUE)
  if(!length(args) %in% c(3,4)) stop("Usage: Rscript scripts/27_RUN_MP_ROBUSTNESS.R <robustness_id> <reference_scenario_number> <block> [prepare|test|full]")
  robustness_id <- args[1]
  scenario_number <- as.integer(args[2])
  block <- as.integer(args[3])
  execution_mode <- if(length(args)==4) args[4] else "prepare"
}
stopifnot(execution_mode %in% c("prepare", "test", "full"))
stopifnot(!is.na(scenario_number), scenario_number>0, !is.na(block), block>0)
test_mode <- execution_mode!="full"
scenario_file <- here("outputs", "mse", "robustness", "design", "Table26_02_MP_robustness_design.csv")
stopifnot(file.exists(scenario_file))
scenario_table <- read.csv(scenario_file, stringsAsFactors=FALSE)
required_cols <- c("robustness_id", "reference_scenario_number", "reference_scenario_id", "robustness_run_id", "Fcap", "Besc", "recruitment_PE_mode", "OM_seasonal_share_mode", paste0("OM_share_Q",1:4), "assessment_error_SSB_F", "forecast_uncertainty", "MP_nFyears", "annual_risk_threshold")
stopifnot(all(required_cols %in% names(scenario_table)))
sc <- scenario_table[scenario_table$robustness_id==robustness_id & scenario_table$reference_scenario_number==scenario_number,,drop=FALSE]
stopifnot(nrow(sc)==1)
stopifnot(robustness_id %in% c("REF", "R-OM", "S-C1", "S-C2", "I0", "I1"))
scenario_id <- sc$robustness_run_id
Fcap <- sc$Fcap
Besc <- sc$Besc
sigmaR_scenario <- sc$recruitment_PE_mode
seasonal_exploitation <- sc$OM_seasonal_share_mode
OM_share <- as.numeric(unlist(sc[1,paste0("OM_share_Q",1:4)],use.names=FALSE))
# Legacy Fprop fields retained only for output compatibility; not MP propf.
Fprop <- OM_share
assessment_error_ON <- as.logical(sc$assessment_error_SSB_F)
assessment_uncertainty <- if(assessment_error_ON) "historical" else "perfect"
forecast_uncertainty <- sc$forecast_uncertainty
# Keep the same HCR branch; CV=0 yields deterministic escapement.
forecast_probability <- TRUE
forecast_alpha <- sc$annual_risk_threshold
stopifnot(is.finite(Fcap), Fcap>0, is.finite(Besc), Besc>0)
stopifnot(sc$MP_nFyears==nFyears, forecast_alpha==0.05)
stopifnot(identical(assessment_error_ON, robustness_id!="I0"))
stopifnot(forecast_uncertainty==if(robustness_id %in% c("I0","I1")) "S0" else "S2")
stopifnot(sigmaR_scenario==if(robustness_id=="R-OM") "source_OM_sigmaR" else "reference_common")
expected_season <- switch(robustness_id, `S-C1`="fixed_C1", `S-C2`="fixed_C2", "conditioned_OM_recent10")
stopifnot(seasonal_exploitation==expected_season)
if(robustness_id %in% c("S-C1","S-C2")) {
  expected_share <- if(robustness_id=="S-C1") c(0.271,0.534,0.164,0.031) else c(0.235,0.309,0.299,0.157)
  stopifnot(isTRUE(all.equal(OM_share,expected_share,tolerance=1e-12)))
} else stopifnot(all(is.na(OM_share)))

n_blocks <- ceiling(ni/block_size)
if(block>n_blocks) stop("Invalid block number. Must be between 1 and ",n_blocks)
out_dir <- if(execution_mode=="full") here("data","mse","robustness_runs",scenario_id) else here("outputs","mse","robustness","validation",execution_mode,paste0("n",Sys.getenv("MSE_N_TEST","10")),scenario_id)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

SR_reference_file <- file.path(rds_boot_dir,"SRs_reference_1000iter_process_error.rds")
SR_file <- if(robustness_id=="R-OM") file.path(rds_boot_dir,"SRs_om_specific_1000iter_process_error.rds") else SR_reference_file
if(!file.exists(SR_file)) stop("Reference recruitment process-error file not found: ",SR_file)

AE_file <- file.path(rds_boot_dir,"assessment_error_reference_1000iter.rds")
if(!file.exists(AE_file)) stop("Assessment-error file not found: ",AE_file)

forecast_CV_file <- file.path(scenario_dir,"forecast_CV_trajectory_map.csv")
if(!file.exists(forecast_CV_file)) stop("Forecast CV trajectory map not found: ",forecast_CV_file)

cat("\n============================================\n")
cat("MSE SCENARIO RUN\n")
cat("Scenario:",scenario_number,"\n")
cat("Scenario ID:",scenario_id,"\n")
cat("Fcap:",Fcap,"\n")
cat("Besc:",Besc,"t\n")
cat("Recruitment PE mode:",sigmaR_scenario,"\n")
cat("Seasonal exploitation:",seasonal_exploitation,"\n")
cat("Seasonal proportions:",paste(round(Fprop,4),collapse=", "),"\n")
cat("Assessment:",assessment_uncertainty,"\n")
cat("Forecast uncertainty:",forecast_uncertainty,"\n")
cat("Forecast probability:",forecast_probability,"\n")
cat("Forecast alpha:",forecast_alpha,"\n")
cat("Block:",block,"/",n_blocks,"\n")
cat("============================================\n")

# ============================================================
# LOAD ITERATION MAP
# ============================================================

iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map_1000.csv"))
stopifnot(all(c("iter1000","om","replicate") %in% names(iteration_map)),
          nrow(iteration_map)==ni,
          n_distinct(iteration_map$iter1000)==ni,
          n_distinct(iteration_map$om)==100)

# Require the original order used by script 09 and source-OM sigmaR.
stopifnot(identical(as.integer(iteration_map$iter1000),seq_len(ni)))
stopifnot("sigmaR" %in% names(iteration_map))
sigma_by_OM <- iteration_map %>% distinct(om,sigmaR)
stopifnot(nrow(sigma_by_OM)==100, all(is.finite(sigma_by_OM$sigmaR)), all(sigma_by_OM$sigmaR>0))
sigmaR_ref <- median(sigma_by_OM$sigmaR)
sigmaR <- if(robustness_id=="R-OM") NA_real_ else sigmaR_ref

forecast_CV_map <- read.csv(forecast_CV_file,stringsAsFactors=FALSE)
stopifnot(all(c("global_iter","forecast_CV") %in% names(forecast_CV_map)))
stopifnot(nrow(forecast_CV_map)==ni)
stopifnot(identical(forecast_CV_map$global_iter,seq_len(ni)))
stopifnot(!anyDuplicated(forecast_CV_map$global_iter))
stopifnot(all(is.finite(forecast_CV_map$forecast_CV)))
stopifnot(all(forecast_CV_map$forecast_CV>0))

# ============================================================
# SELECT BLOCK
# ============================================================

first_iter <- (block-1)*block_size+1
last_iter <- min(block*block_size,ni)
ii <- first_iter:last_iter
n_test <- as.integer(Sys.getenv("MSE_N_TEST", "10"))
stopifnot(length(n_test)==1, !is.na(n_test), n_test>=1, n_test<=block_size)
if(test_mode) ii <- head(ii,n_test)
nb <- length(ii)

run_map <- iteration_map %>% 
            dplyr::filter(iter1000 %in% ii) %>% 
            dplyr::arrange(match(iter1000,ii)) %>% 
            dplyr::mutate(iter_local=seq_len(dplyr::n()))

forecast_CV_b <- forecast_CV_map %>% 
                  dplyr::filter(global_iter %in% ii) %>% 
                  dplyr::arrange(match(global_iter,ii))

stopifnot(nrow(forecast_CV_b)==nb)
stopifnot(all(forecast_CV_b$global_iter==ii))

run_map <- run_map %>% dplyr::left_join(forecast_CV_b,by=c("iter1000"="global_iter"))

stopifnot(nrow(run_map)==nb,all(run_map$iter1000==ii),
          dplyr::n_distinct(run_map$iter_local)==nb)
stopifnot(all(is.finite(run_map$forecast_CV)))
stopifnot(all(run_map$forecast_CV>0))
run_map$forecast_CV_reference <- run_map$forecast_CV
if(forecast_uncertainty=="S0") run_map$forecast_CV <- 0
stopifnot(all(run_map$forecast_CV>=0))
if(forecast_uncertainty=="S0") stopifnot(all(run_map$forecast_CV==0))
if(forecast_uncertainty=="S2") stopifnot(identical(run_map$forecast_CV,run_map$forecast_CV_reference))
run_map$sigmaR_used <- if(robustness_id=="R-OM") run_map$sigmaR else sigmaR_ref
if(test_mode) stopifnot(nb==n_test)
if(!test_mode) stopifnot(nb==block_size)

cat("\n============================================\n")
cat("RUNNING BLOCK",block,"OF",n_blocks,"\n")
cat("Global iterations:",min(ii),"-",max(ii),"\n")
cat("Number of trajectories:",nb,"\n")
cat("Test mode:",test_mode,"\n")
cat("Forecast CV median:",median(run_map$forecast_CV),"\n")
cat("Forecast CV range:",range(run_map$forecast_CV),"\n")
cat("============================================\n")

# ============================================================
# LOAD CONDITIONED OM
# ============================================================

biols <- readRDS(file.path(rds_boot_dir,"biols_conditioned_1000iter.rds"))
SRs <- readRDS(SR_file)
fleets <- readRDS(file.path(rds_boot_dir,"fleets_conditioned_1000iter.rds"))
fleets.ctrl <- readRDS(file.path(rds_boot_dir,"fleets_ctrl_reference_1000iter.rds"))
assessment_error <- readRDS(AE_file)
if(!assessment_error_ON) assessment_error$multiplier[] <- 1

stopifnot(dim(assessment_error$multiplier)[1]==30)
stopifnot(dim(assessment_error$multiplier)[2]==ni)
stopifnot(dim(assessment_error$multiplier)[3]==3)
stopifnot(identical(dimnames(assessment_error$multiplier)$year,as.character(2024:2053)))
stopifnot(identical(dimnames(assessment_error$multiplier)$variable,c("SSB","F","R")))


stopifnot("ANE" %in% names(biols),"ANE" %in% names(SRs),"SEINE" %in% 
            names(fleets),dim(biols$ANE@n)[6]==ni,
          dim(SRs$ANE@params)[4]==ni,
          dim(SRs$ANE@uncertainty)[6]==ni,dim(fleets$SEINE@effort)[6]==ni,
          dim(catch.q(fleets$SEINE@metiers$ALL@catches$ANE))[6]==ni)

# ============================================================
# SUBSETTING FUNCTIONS
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
  for(s in c("rec","ssb","uncertainty","proportion")) slot(out,s,check=FALSE) <- sub_flq(slot(x,s),ii)
  pars <- x@params[,,,ii,drop=FALSE]
  dimnames(pars)$iter <- as.character(seq_along(ii))
  slot(out,"params",check=FALSE) <- pars
  validObject(out)
  out
}

sub_fleet <- function(x,ii) {
  out <- x
  for(s in slotNames(x)) {
    z <- slot(x,s)
    if(inherits(z,"FLQuant")) slot(out,s,check=FALSE) <- sub_flq(z,ii)
  }
  met <- x@metiers$ALL
  for(s in slotNames(met)) {
    z <- slot(met,s)
    if(inherits(z,"FLQuant")) slot(met,s,check=FALSE) <- sub_flq(z,ii)
  }
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

sub_assessment_error <- function(x,ii) {
  out <- x
  out$multiplier <- x$multiplier[,ii,,drop=FALSE]
  dimnames(out$multiplier)$iter <- as.character(seq_along(ii))
  out
}
# ============================================================
# CONDITIONED OBJECTS FOR CURRENT BLOCK
# ============================================================

biols_b <- biols
biols_b$ANE <- sub_FLBiol(biols$ANE,ii)

SRs_b <- SRs
SRs_b$ANE <- sub_FLSR(SRs$ANE,ii)

fleets_b <- fleets
fleets_b$SEINE <- sub_fleet(fleets$SEINE,ii)

fleets.ctrl_b <- fleets.ctrl
fleets.ctrl_b$seasonal.share[[1]] <- sub_flq(fleets.ctrl$seasonal.share[[1]],ii)

assessment_error_b <- sub_assessment_error(assessment_error,ii)

stopifnot(dim(assessment_error_b$multiplier)[2]==nb)
stopifnot(all(is.finite(assessment_error_b$multiplier)))
stopifnot(all(assessment_error_b$multiplier>0))

stopifnot(dim(biols_b$ANE@n)[6]==nb,
          dim(SRs_b$ANE@params)[4]==nb,
          dim(SRs_b$ANE@uncertainty)[6]==nb,dim(fleets_b$SEINE@effort)[6]==nb,
          dim(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE))[6]==nb,
          dim(fleets.ctrl_b$seasonal.share[[1]])[6]==nb)

cat("Conditioned objects subset to current block: OK\n")

# ============================================================
# ROBUSTNESS INPUT MODIFICATIONS AND PAIRED VERIFICATION
# ============================================================

shares_reference <- fleets.ctrl_b$seasonal.share[[1]]
if(robustness_id %in% c("S-C1","S-C2")) {
  for(q in seq_len(ns)) {
    fleets.ctrl_b$seasonal.share[[1]][,ac(proj.yr:last.yr),,ac(q),,] <- OM_share[q]
    actual <- as.numeric(fleets.ctrl_b$seasonal.share[[1]][,ac(proj.yr:last.yr),,ac(q),,])
    stopifnot(all(abs(actual-OM_share[q])<1e-12))
  }
}
share_years <- as.numeric(dimnames(shares_reference)$year)
hist_share_years <- ac(share_years[share_years<=last.obs.yr])
stopifnot(identical(fleets.ctrl_b$seasonal.share[[1]][,hist_share_years,,,,],shares_reference[,hist_share_years,,,,]))
seasonal_check <- seasonSums(fleets.ctrl_b$seasonal.share[[1]][,ac(proj.yr:last.yr),,,,])
stopifnot(all(is.finite(as.numeric(seasonal_check))))
stopifnot(all(abs(as.numeric(seasonal_check)-1)<1e-8))
if(!robustness_id %in% c("S-C1","S-C2")) stopifnot(identical(fleets.ctrl_b$seasonal.share[[1]],shares_reference))

# Compare the complete R-OM FLSRsim against REF. Only uncertainty
# in projection Q3 may differ; BH parameters and other slots stay REF.
SRs_ref_b <- readRDS(SR_reference_file)
SRs_ref_b$ANE <- sub_FLSR(SRs_ref_b$ANE,ii)
for(slot_name in setdiff(slotNames(SRs_ref_b$ANE),"uncertainty")) {
  stopifnot(identical(slot(SRs_b$ANE,slot_name),slot(SRs_ref_b$ANE,slot_name)))
}
unc_expected <- SRs_ref_b$ANE@uncertainty
z_file <- file.path(rds_boot_dir,"recruitment_standard_normal_deviates_1000iter.rds")
stopifnot(file.exists(z_file))
z_R <- readRDS(z_file)
stopifnot(identical(rownames(z_R),ac(proj.yr:last.yr)))
stopifnot(identical(colnames(z_R),ac(seq_len(ni))))
z_b <- z_R[,ii,drop=FALSE]
expected_ref <- exp(sigmaR_ref*z_b-0.5*sigmaR_ref^2)
actual_ref <- as.numeric(unc_expected[,ac(proj.yr:last.yr),,ac(3),,])
stopifnot(isTRUE(all.equal(actual_ref,as.numeric(expected_ref),tolerance=1e-12)))
expected_R <- exp(sweep(sweep(z_b,2,run_map$sigmaR_used,"*"),2,0.5*run_map$sigmaR_used^2,"-"))
if(robustness_id=="R-OM") {
  for(i in seq_len(nb)) unc_expected[,ac(proj.yr:last.yr),,ac(3),,i] <- expected_R[,i]
}
stopifnot(isTRUE(all.equal(SRs_b$ANE@uncertainty,unc_expected,tolerance=1e-12)))
stopifnot(isTRUE(all.equal(as.numeric(SRs_b$ANE@uncertainty[,ac(proj.yr:last.yr),,ac(3),,]),as.numeric(expected_R),tolerance=1e-12)))

AE_reference <- sub_assessment_error(readRDS(AE_file),ii)
if(assessment_error_ON) stopifnot(identical(assessment_error_b,AE_reference))
if(!assessment_error_ON) stopifnot(all(assessment_error_b$multiplier==1))
cat("Seasonal assignment, historical preservation and recruitment/AE/CV pairing: PASS\n")

# ============================================================
# CONDITIONED FLSTOCK FOR perfectObs4seas
# ============================================================

ane.stock.all <- readRDS(file.path(rds_boot_dir,"ane_stock_conditioned_1000iter.rds"))
ane.stock <- ane.stock.all

for(s in slotNames(ane.stock.all)) {
  z <- slot(ane.stock.all,s)
  if(inherits(z,"FLQuant")) slot(ane.stock,s,check=FALSE) <- sub_flq(z,ii)
}

range(ane.stock,"minfbar") <- 3
range(ane.stock,"maxfbar") <- 3
validObject(ane.stock)

# ============================================================
# CONTROLS
# ============================================================

biols.ctrl <- create.biols.ctrl(stksnames=stks,growth.model="ASPG_Baranov")
main.ctrl <- list(sim.years=c(initial=last.obs.yr,final=last.yr))
covars <- NULL
covars.ctrl <- NULL
indices <- NULL

# ============================================================
# PERFECT OBSERVATION
# ============================================================

flq.ANE <- FLQuant(dimnames=list(age="all",year=first.yr:last.yr,unit=1,season=1:ns,iter=1:nb))
obs.ctrl <- create.obs.ctrl(stksnames=stks,flq.ANE=flq.ANE,stkObs.models="perfectObs")
obs.ctrl$ANE$obs.curryr <- TRUE

# ============================================================
# NO INTERNAL ASSESSMENT MODEL
# Assessment error is applied within the HCR
# ============================================================

assess.ctrl <- create.assess.ctrl(stksnames=stks,assess.models="NoAssessment")
assess.ctrl$ANE$ass.curryr <- TRUE

# ============================================================
# ADVICE DATA
# ============================================================

ANE_advice.TAC.flq <- FLQuant(c(seasonSums(quantSums(catchWStock(fleets_b,stock="ANE")[,ac(first.yr:last.yr),,,,]))),
                              dimnames=list(quant=stks,year=first.yr:last.yr),iter=nb)
ANE_advice.TAC.flq[,ac(1989:2018),] <- NA
ANE_advice.TAC.flq[,ac(2019:2024),] <- c(5278,8856,9459,4383,1892,1733)
ANE_advice.TAC.flq[,ac(proj.yr:last.yr),] <- NA

ANE_advice.quota.share.flq <- FLQuant(1,dimnames=list(quant=stks,year=first.yr:last.yr),iter=nb)
ANE_advice.avg.yrs <- c(2023,2024)

stksTac.data <- list(ANE=ls(pattern="^ANE"))
yrs <- c(first.yr=first.yr,proj.yr=proj.yr,last.yr=last.yr)

advice_HCR <- create.advice.data(yrs,ns,nb,stksTac.data,fleets_b)
units(advice_HCR$TAC) <- "t"

cat("Initial future TAC from create.advice.data():\n")
print(advice_HCR$TAC[,ac(proj.yr:last.yr),,,,])

advice_HCR$TAC[,ac(proj.yr:last.yr),,,,] <- NA
stopifnot(all(is.na(as.numeric(advice_HCR$TAC[,ac(proj.yr:last.yr),,,,]))))
# ============================================================
# REFERENCE POINTS
# ============================================================

ref.pts <- matrix(c(Blim,Bpa),ncol=1,dimnames=list(c("Blim","Bpa"),"value"))

# ============================================================
# FCAP + BESCAPEMENT ADVICE CONTROL
# ============================================================

advice_HCR.ctrl <- list(ANE=list(HCR.model="annualTAC",
                                 Fcap=Fcap,
                                 Besc=Besc,
                                 nFyears=nFyears,
                                 spawn.season=2,
                                 rec.season=3,
                                 nyears=1,
                                 wts.nyears=3,
                                 fbar.nyears=3,
                                 f.rescale=TRUE,
                                 ref.pts=ref.pts,
                                 sr=SRs_b$ANE,
                                 catch.q=catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE),
                                 assessment.error=assessment_error_b,
                                 assess_error_R=FALSE,
                                 use.probability=forecast_probability,
                                 forecast.CV=run_map$forecast_CV,
                                 alpha=forecast_alpha,
                                 diagnostics=FALSE))

stopifnot(length(advice_HCR.ctrl$ANE$forecast.CV)==nb)
stopifnot(identical(as.numeric(advice_HCR.ctrl$ANE$forecast.CV),as.numeric(run_map$forecast_CV)))

cat("Forecast CV pairing with HCR control: OK\n")
# ============================================================
# PRE-RUN CHECKS
# ============================================================

cat("\n============================================\n")
cat("PRE-RUN CHECKS\n")
cat("============================================\n")

stopifnot(dim(biols_b$ANE@n)[6]==nb,
          dim(SRs_b$ANE@params)[4]==nb,
          dim(SRs_b$ANE@uncertainty)[6]==nb,
          dim(fleets_b$SEINE@effort)[6]==nb,
          dim(fleets_b$SEINE@metiers$ALL@vcost)[6]==nb,
          dim(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE))[6]==nb,
          dim(fleets.ctrl_b$seasonal.share[[1]])[6]==nb,
          all(is.finite(as.numeric(fleets_b$SEINE@metiers$ALL@vcost))),
          all(is.finite(as.numeric(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE)))))

cat("Block inputs: OK\n")
cat("Fcap =",Fcap,"\n")
cat("Besc =",Besc,"t\n")
cat("Assessment error =",ifelse(assessment_error_ON,"ON (SSB + F)","OFF / perfect"),"\n")
cat("Assessment error R: OFF\n")
cat("Forecast probability =",forecast_probability,"\n")
cat("Forecast alpha =",forecast_alpha,"\n")
cat("Forecast CV median/range =",median(run_map$forecast_CV),range(run_map$forecast_CV),"\n")
stopifnot(length(run_map$forecast_CV)==nb)
if(assessment_error_ON) cat("SSB multiplier range =",range(assessment_error_b$multiplier[,,"SSB"]),"\n")
if(assessment_error_ON) cat("F multiplier range =",range(assessment_error_b$multiplier[,,"F"]),"\n")
if(!assessment_error_ON) stopifnot(all(assessment_error_b$multiplier==1))

# ============================================================
# CHECK FUTURE RECRUITMENT PROCESS ERROR
# ============================================================

future_unc <- as.numeric(SRs_b$ANE@uncertainty[,ac(proj.yr:last.yr),,ac(3),,])

stopifnot(length(future_unc)>0,all(is.finite(future_unc)),all(future_unc>0))

log_future_unc <- log(future_unc)

cat("Recruitment multiplier mean =",mean(future_unc),"\n")
cat("Recruitment multiplier median =",median(future_unc),"\n")
cat("Recruitment multiplier range =",range(future_unc),"\n")
cat("Recruitment log-error mean =",mean(log_future_unc),"\n")
cat("Recruitment log-error SD =",sd(log_future_unc),"\n")
cat("Recruitment PE mode:",sigmaR_scenario,"\n")
# ============================================================
# RUN FLBEIA
# ============================================================

cat("\n============================================\n")
cat("RUNNING CANDIDATE MP CLOSED-LOOP MSE\n")
cat("Fcap:",Fcap,"\n")
cat("Besc:",Besc,"t\n")
cat("Block:",block,"/",n_blocks,"\n")
cat("============================================\n")

validation_file <- file.path(out_dir,sprintf("block_%02d_input_validation.rds",block))
if(file.exists(validation_file)) stop("Validation output already exists: ",validation_file)
saveRDS(list(robustness_id=robustness_id, reference_scenario_id=sc$reference_scenario_id,
             execution_mode=execution_mode, scenario=sc, run_map=run_map,
             recruitment_Q3=SRs_b$ANE@uncertainty[,ac(proj.yr:last.yr),,ac(3),,],
             seasonal_share=fleets.ctrl_b$seasonal.share[[1]],
             assessment_error=assessment_error_b, advice_control=advice_HCR.ctrl,
             input_checks="PASS", closed_loop_validation="PENDING"),validation_file)
if(execution_mode=="prepare") {
  cat("\nINPUT PREPARATION COMPLETE. No FLBEIA simulation executed.\n")
  cat("Validation file:",validation_file,"\n")
} else {
OM_target <- file.path(out_dir,sprintf("block_%02d.rds",block))
if(file.exists(OM_target)) stop("Run output already exists: ",OM_target)

start_time <- Sys.time()

OM_b <- FLBEIA(biols=biols_b,
               fleets=fleets_b,
               SRs=SRs_b,
               BDs=NULL,
               covars=covars,
               indices=indices,
               advice=advice_HCR,
               main.ctrl=main.ctrl,
               biols.ctrl=biols.ctrl,
               fleets.ctrl=fleets.ctrl_b,
               covars.ctrl=covars.ctrl,
               obs.ctrl=obs.ctrl,
               assess.ctrl=assess.ctrl,
               advice.ctrl=advice_HCR.ctrl)

end_time <- Sys.time()

cat("\nRuntime:",round(as.numeric(difftime(end_time,start_time,units="mins")),2),"minutes\n")

# ============================================================
# POST-RUN STRUCTURAL CHECKS
# ============================================================

stopifnot(dim(OM_b$biols$ANE@n)[6]==nb)

future_n <- as.numeric(OM_b$biols$ANE@n[,ac(proj.yr:last.yr),,,,])
stopifnot(length(future_n)>0,all(is.finite(future_n)),all(future_n>=0))

future_R <- as.numeric(OM_b$biols$ANE@n["0",ac(proj.yr:last.yr),,"3",,])
stopifnot(length(future_R)>0,all(is.finite(future_R)),all(future_R>=0))

cat("Projected abundance: OK\n")
cat("Projected recruitment: OK\n")

# ============================================================
# CHECK MP ADVICE
# ============================================================

TAC_proj <- as.numeric(OM_b$advice$TAC["ANE",ac(proj.yr:last.yr),,,,])
stopifnot(length(TAC_proj)>0,all(is.finite(TAC_proj)),all(TAC_proj>=0))

cat("Future TAC advice: finite and non-negative: OK\n")
cat("Median TAC =",median(TAC_proj),"\n")
cat("TAC range =",range(TAC_proj),"\n")

# ============================================================
# CHECK REALISED CATCH
# ============================================================

cobj_out <- OM_b$fleets$SEINE@metiers$ALL@catches$ANE
C_real <- seasonSums(quantSums(cobj_out@landings))
C_real_proj <- as.numeric(C_real[,ac(proj.yr:last.yr),,,,])

stopifnot(length(C_real_proj)>0,all(is.finite(C_real_proj)),all(C_real_proj>=0))

cat("Realised catch: finite and non-negative: OK\n")
cat("Median realised catch =",median(C_real_proj),"\n")
cat("Realised catch range =",range(C_real_proj),"\n")

# ============================================================
# CHECK TAC IMPLEMENTATION
# ============================================================

TAC_realisation <- as.numeric(C_real[,ac(proj.yr:last.yr),,,,])-TAC_proj
TAC_ratio <- ifelse(TAC_proj>0,as.numeric(C_real[,ac(proj.yr:last.yr),,,,])/TAC_proj,NA_real_)

stopifnot(all(is.finite(TAC_realisation)))

cat("TAC implementation: OK\n")
cat("Catch - TAC median/range =",median(TAC_realisation),range(TAC_realisation),"\n")
cat("Catch/TAC median/range =",median(TAC_ratio,na.rm=TRUE),range(TAC_ratio,na.rm=TRUE),"\n")

# ============================================================
# CHECK CLOSED-LOOP OUTPUTS
# ============================================================

proj_years <- ac(proj.yr:last.yr)
Fadv_proj <- as.numeric(OM_b$advice$Fadv["ANE",proj_years,,,,])

stopifnot(length(Fadv_proj)>0,all(is.finite(Fadv_proj)),all(Fadv_proj>=0),all(Fadv_proj<=Fcap+1e-8))

cat("Fadv: finite and non-negative: OK\n")
cat("Fadv <= Fcap: OK\n")
cat("Median Fadv =",median(Fadv_proj),"\n")
cat("Fadv range =",range(Fadv_proj),"\n")


# ============================================================
# SAVE OUTPUT
# ============================================================

OM_file <- file.path(out_dir,sprintf("block_%02d.rds",block))
map_file <- file.path(out_dir,sprintf("block_%02d_map.csv",block))
metadata_file <- file.path(out_dir,sprintf("block_%02d_metadata.csv",block))

saveRDS(OM_b,OM_file)
write.csv(run_map,map_file,row.names=FALSE)

run_metadata <- data.frame(robustness_id=robustness_id,
                           reference_scenario_id=sc$reference_scenario_id,
                           execution_mode=execution_mode,
                           scenario_number=scenario_number,
                           scenario_id=scenario_id,
                           Fcap=Fcap,
                           Besc=Besc,
                           Blim=Blim,
                           Bpa=Bpa,
                           sigmaR_scenario=sigmaR_scenario,
                           sigmaR=sigmaR,
                           seasonal_exploitation=seasonal_exploitation,
                           OM_share_Q1=OM_share[1],
                           OM_share_Q2=OM_share[2],
                           OM_share_Q3=OM_share[3],
                           OM_share_Q4=OM_share[4],
                           assessment_uncertainty=assessment_uncertainty,
                           assessment_error=assessment_error_ON,
                           assessment_error_R=FALSE,forecast_uncertainty=forecast_uncertainty,
                           forecast_probability=forecast_probability,
                           forecast_alpha=forecast_alpha,
                           forecast_CV_median=median(run_map$forecast_CV),
                           forecast_CV_min=min(run_map$forecast_CV),
                           forecast_CV_max=max(run_map$forecast_CV),
                           first_year=first.yr,
                           last_historical_year=last.obs.yr,
                           first_projection_year=proj.yr,
                           last_projection_year=last.yr,
                           n_trajectories=nb,
                           n_conditioned_OMs=n_distinct(run_map$om),
                           block=block,
                           test_mode=test_mode,
                           runtime_minutes=as.numeric(difftime(end_time,start_time,units="mins")))
write.csv(run_metadata,metadata_file,row.names=FALSE)

# ============================================================
# COMPLETION MESSAGE
# ============================================================

cat("\n============================================\n")
cat("MSE SCENARIO BLOCK COMPLETED\n")
cat("============================================\n")
cat("Block:",block,"/",n_blocks,"\n")
cat("Global iterations:",min(ii),"-",max(ii),"\n")
cat("Historical:",first.yr,"-",last.obs.yr,"\n")
cat("Projection:",proj.yr,"-",last.yr,"\n")
cat("Conditioned OMs:",n_distinct(run_map$om),"\n")
cat("Trajectories:",nb,"\n")
cat("Fcap:",Fcap,"\n")
cat("Besc:",Besc,"t\n")
cat("Assessment error =",ifelse(assessment_error_ON,"ON (SSB + F)","OFF / perfect"),"\n")
cat("Assessment error R: OFF\n")
cat("Forecast uncertainty:",forecast_uncertainty,"\n")
cat("Forecast probability:",forecast_probability,"\n")
cat("Forecast alpha:",forecast_alpha,"\n")
cat("Forecast CV median/range:",median(run_map$forecast_CV),range(run_map$forecast_CV),"\n")
if(assessment_error_ON) cat("SSB multiplier range =",range(assessment_error_b$multiplier[,,"SSB"]),"\n")
if(assessment_error_ON) cat("F multiplier range =",range(assessment_error_b$multiplier[,,"F"]),"\n")
if(!assessment_error_ON) stopifnot(all(assessment_error_b$multiplier==1))
cat("Runtime:",round(as.numeric(difftime(end_time,start_time,units="mins")),2),"minutes\n")




} # execution_mode != prepare
