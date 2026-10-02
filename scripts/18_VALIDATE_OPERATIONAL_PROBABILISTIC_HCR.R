# ============================================================
# 18_VALIDATE_OPERATIONAL_PROBABILISTIC_HCR.R
#
# Operational closed-loop validation of the probabilistic
# Fcap-Bescapement HCR for anchovy 9a South.
#
# Reference implementation:
#   FcapBpaHCR_ane9aS_CORRECTED.R
#
# Probabilistic extension:
#   FcapBpaHCR_ane9aS_PROBABILITY_OPTION.R
#
# Validation sequence:
#
#   1. Reference deterministic HCR
#   2. Probabilistic implementation with probability OFF
#   3. Probabilistic implementation with probability ON and CV = 0
#   4. S1: current-assessment forecast uncertainty
#   5. S2: contemporary median forecast uncertainty
#
# Tests 1-3 verify backward compatibility:
#
#   Probability OFF == Reference
#   Probability ON, CV=0 == Reference
#
# Tests 4-5 verify that forecast uncertainty modifies management
# advice only through the probabilistic escapement constraint.
#
# All scenarios use:
#   - the same conditioned OM trajectory
#   - the same recruitment process error
#   - perfect observation
#   - no assessment error
#   - a single management step: 2024 -> 2025
#
# Therefore, differences among scenarios are attributable only
# to the forecast-uncertainty implementation.
# ============================================================

rm(list=ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)

# ============================================================
# DIRECTORIES
# ============================================================

rds_dir <- here("data","Rdata")
function_dir <- here("functions")

# ============================================================
# SETTINGS
# ============================================================

first.yr <- 1989
last.obs.yr <- 2024
proj.yr <- 2025
last.yr <- 2025

ns <- 4
stks <- "ANE"

Fcap <- 1.5
Blim <- 4721
Bpa <- 6561
Besc <- 6561

nFyears <- 5
alpha <- 0.05


forecast_uncertainty <- readRDS(file.path(rds_dir,"forecast_uncertainty_characterisation.rds"))

CV_S0 <- 0
CV_S1 <- forecast_uncertainty$base$SSB_CV
CV_S2 <- median(forecast_uncertainty$contemporary$SSB_2025_CV)

stopifnot(length(CV_S1)==1)
stopifnot(is.finite(CV_S1))
stopifnot(CV_S1>0)

stopifnot(length(CV_S2)==1)
stopifnot(is.finite(CV_S2))
stopifnot(CV_S2>0)

test_iter <- 1

# ============================================================
# FILES
# ============================================================

reference_file <- file.path(function_dir,"FcapBpaHCR_ane9aS_CORRECTED.R")
probability_file <- file.path(function_dir,"FcapBpaHCR_ane9aS_PROBABILITY_OPTION.R")

input_files <- file.path(rds_dir,c(
  "biols_conditioned_1000iter.rds",
  "SRs_reference_1000iter_process_error.rds",
  "fleets_conditioned_1000iter.rds",
  "fleets_ctrl_reference_1000iter.rds",
  "ane_stock_conditioned_1000iter.rds",
  "assessment_error_reference_1000iter.rds"
))

stopifnot(file.exists(reference_file))
stopifnot(file.exists(probability_file))
stopifnot(all(file.exists(input_files)))

cat("\nInput files found: OK\n")

# ============================================================
# LOAD CUSTOM OBSERVATION FUNCTION
# ============================================================

source(file.path(function_dir,"perfectObs4seas.R"))

# ============================================================
# LOAD CONDITIONED OBJECTS
# ============================================================

biols_all <- readRDS(file.path(rds_dir,"biols_conditioned_1000iter.rds"))
SRs_all <- readRDS(file.path(rds_dir,"SRs_reference_1000iter_process_error.rds"))
fleets_all <- readRDS(file.path(rds_dir,"fleets_conditioned_1000iter.rds"))
fleets.ctrl_all <- readRDS(file.path(rds_dir,"fleets_ctrl_reference_1000iter.rds"))
assessment_error_all <- readRDS(file.path(rds_dir,"assessment_error_reference_1000iter.rds"))

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

# ============================================================
# SELECT ONE CONDITIONED TRAJECTORY
# ============================================================

ii <- test_iter
nb <- 1

biols_test <- biols_all
biols_test$ANE <- sub_FLBiol(biols_all$ANE,ii)

SRs_test <- SRs_all
SRs_test$ANE <- sub_FLSR(SRs_all$ANE,ii)

fleets_test <- fleets_all
fleets_test$SEINE <- sub_fleet(fleets_all$SEINE,ii)

fleets.ctrl_test <- fleets.ctrl_all
fleets.ctrl_test$seasonal.share[[1]] <- sub_flq(fleets.ctrl_all$seasonal.share[[1]],ii)


cat("\nConditioned trajectory:",test_iter,"\n")

# ============================================================
# ASSESSMENT ERROR OFF
# ============================================================

assessment_error_test <- assessment_error_all
assessment_error_test$multiplier <- assessment_error_all$multiplier[,ii,,drop=FALSE]
dimnames(assessment_error_test$multiplier)$iter <- "1"
assessment_error_test$multiplier[] <- 1

stopifnot(all(assessment_error_test$multiplier==1))

cat("Assessment error: OFF\n")

# ============================================================
# CHECK CONDITIONED OBJECTS
# ============================================================

stopifnot(dim(biols_test$ANE@n)[6]==1)
stopifnot(dim(SRs_test$ANE@params)[4]==1)
stopifnot(dim(SRs_test$ANE@uncertainty)[6]==1)
stopifnot(dim(fleets_test$SEINE@effort)[6]==1)
stopifnot(dim(catch.q(fleets_test$SEINE@metiers$ALL@catches$ANE))[6]==1)

seasonal_check <- seasonSums(fleets.ctrl_test$seasonal.share[[1]][,ac(proj.yr),,,,])

stopifnot(all(abs(as.numeric(seasonal_check)-1)<1e-8))

cat("Conditioned objects: OK\n")

# ============================================================
# FLBEIA CONTROLS
# ============================================================

biols.ctrl <- create.biols.ctrl(stksnames=stks,growth.model="ASPG_Baranov")

main.ctrl <- list(sim.years=c(initial=last.obs.yr,final=last.yr))

covars <- NULL
covars.ctrl <- NULL
indices <- NULL

# ============================================================
# PERFECT OBSERVATION
# ============================================================

flq.ANE <- FLQuant(dimnames=list(age="all",year=first.yr:last.yr,unit=1,season=1:ns,iter=1))

obs.ctrl <- create.obs.ctrl(stksnames=stks,flq.ANE=flq.ANE,stkObs.models="perfectObs")

obs.ctrl$ANE$obs.curryr <- TRUE

# ============================================================
# NO INTERNAL ASSESSMENT MODEL
# ============================================================

assess.ctrl <- create.assess.ctrl(stksnames=stks,assess.models="NoAssessment")

assess.ctrl$ANE$ass.curryr <- TRUE

# ============================================================
# ADVICE DATA
# ============================================================

ANE_advice.TAC.flq <- FLQuant(c(seasonSums(quantSums(catchWStock(fleets_test,stock="ANE")[,ac(first.yr:last.yr),,,,]))),dimnames=list(quant=stks,year=first.yr:last.yr),iter=1)

ANE_advice.TAC.flq[,ac(1989:2018),] <- NA

ANE_advice.TAC.flq[,ac(2019:2024),] <- c(5278,8856,9459,4383,1892,1733)

ANE_advice.TAC.flq[,ac(proj.yr:last.yr),] <- NA

ANE_advice.quota.share.flq <- FLQuant(1,dimnames=list(quant=stks,year=first.yr:last.yr),iter=1)

ANE_advice.avg.yrs <- c(2023,2024)

stksTac.data <- list(ANE=ls(pattern="^ANE"))

yrs <- c(first.yr=first.yr,proj.yr=proj.yr,last.yr=last.yr)

advice_template <- create.advice.data(yrs,ns,1,stksTac.data,fleets_test)

units(advice_template$TAC) <- "t"

advice_template$TAC[,ac(proj.yr:last.yr),,,,] <- NA

stopifnot(all(is.na(as.numeric(advice_template$TAC[,ac(proj.yr:last.yr),,,,]))))

# ============================================================
# REFERENCE POINTS
# ============================================================

ref.pts <- matrix(c(Blim,Bpa),ncol=1,dimnames=list(c("Blim","Bpa"),"value"))

# ============================================================
# COMMON ADVICE CONTROL
# ============================================================

make_advice_ctrl <- function() {
  list(ANE=list(
    HCR.model="annualTAC",
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
    sr=SRs_test$ANE,
    catch.q=catch.q(fleets_test$SEINE@metiers$ALL@catches$ANE),
    assessment.error=assessment_error_test,
    assess_error_R=FALSE,
    diagnostics=TRUE
  ))
}

# ============================================================
# RUN ONE OPERATIONAL SCENARIO
#
# Each scenario starts from the same conditioned trajectory.
# The requested HCR implementation is loaded immediately before
# running FLBEIA to ensure that annualTAC() uses the correct
# management procedure.
#
# Returned quantities:
#
#   Fadv  : fishing mortality selected by the MP
#   TAC   : catch advice produced by the MP
#   SSB   : realised OM spawning biomass in Q2 of 2025
#   Catch : realised OM catch during 2025
#
# Note:
#   SSB is the realised operating-model SSB, not the internal
#   one-year forecast SSB used by the HCR.
# ============================================================

run_scenario <- function(label,hcr_file,use_probability=NULL,forecast_CV=NULL,diagnostics=FALSE) {
  
  cat("\n============================================================\n")
  cat(label,"\n")
  cat("============================================================\n")
  
  source(hcr_file)
  
  environment(FcapBpaHCR_ane) <- asNamespace("FLBEIA")
  assignInNamespace("annualTAC",FcapBpaHCR_ane,ns="FLBEIA")
  
  advice.ctrl <- make_advice_ctrl()
  
  advice.ctrl$ANE$diagnostics <- diagnostics
  
  if(!is.null(use_probability)) advice.ctrl$ANE$use.probability <- use_probability
  if(!is.null(forecast_CV)) advice.ctrl$ANE$forecast.CV <- forecast_CV
  if(!is.null(use_probability)) advice.ctrl$ANE$alpha <- alpha
  
  OM <- FLBEIA(
    biols=biols_test,
    fleets=fleets_test,
    SRs=SRs_test,
    BDs=NULL,
    covars=covars,
    indices=indices,
    advice=advice_template,
    main.ctrl=main.ctrl,
    biols.ctrl=biols.ctrl,
    fleets.ctrl=fleets.ctrl_test,
    covars.ctrl=covars.ctrl,
    obs.ctrl=obs.ctrl,
    assess.ctrl=assess.ctrl,
    advice.ctrl=advice.ctrl
  )
  
  Fadv <- as.numeric(OM$advice$Fadv["ANE",ac(proj.yr),,,,])
  TAC <- as.numeric(OM$advice$TAC["ANE",ac(proj.yr),,,,])
  
  biol <- OM$biols$ANE
  Mat <- predict(biol@mat)
  
  SSB <- as.numeric(quantSums(biol@n[,ac(proj.yr),,"2",,]*biol@wt[,ac(proj.yr),,"2",,]*Mat[,ac(proj.yr),,"2",,]))
  
  cobj <- OM$fleets$SEINE@metiers$ALL@catches$ANE
  
  Catch <- as.numeric(seasonSums(quantSums(cobj@landings[,ac(proj.yr),,,,])))
  
  stopifnot(length(Fadv)==1)
  stopifnot(length(TAC)==1)
  stopifnot(length(SSB)==1)
  stopifnot(length(Catch)==1)
  
  stopifnot(is.finite(Fadv))
  stopifnot(is.finite(TAC))
  stopifnot(is.finite(SSB))
  stopifnot(is.finite(Catch))
  
  stopifnot(Fadv>=0)
  stopifnot(Fadv<=Fcap+1e-8)
  stopifnot(TAC>=0)
  stopifnot(SSB>=0)
  stopifnot(Catch>=0)
  
  cat("Fadv  =",Fadv,"\n")
  cat("TAC   =",TAC,"\n")
  cat("SSB   =",SSB,"\n")
  cat("Catch =",Catch,"\n")
  
  list(
    label=label,
    OM=OM,
    Fadv=Fadv,
    TAC=TAC,
    SSB=SSB,
    Catch=Catch
  )
}

# ============================================================
# TEST 1 - REFERENCE
# ============================================================

res_reference <- run_scenario(
  label="REFERENCE - CORRECTED HCR",
  hcr_file=reference_file
)

# ============================================================
# TEST 2 - PROBABILISTIC IMPLEMENTATION, PROBABILITY OFF
#
# Backward-compatibility test:
# the probabilistic implementation must reproduce the reference
# deterministic HCR exactly when probability is disabled.
# ============================================================
res_off <- run_scenario(
  label="PROBABILITY OPTION - OFF",
  hcr_file=probability_file,
  use_probability=FALSE,
  forecast_CV=0
)

# ============================================================
# TEST 3 - PROBABILISTIC IMPLEMENTATION, CV = 0
#
# Zero-uncertainty limit:
# when CV = 0, the probabilistic escapement criterion reduces
# exactly to the deterministic criterion.
# ============================================================

res_CV0 <- run_scenario(
  label="PROBABILITY OPTION - CV = 0",
  hcr_file=probability_file,
  use_probability=TRUE,
  forecast_CV=CV_S0
)

# ============================================================
# BASELINE EQUIVALENCE
# ============================================================

stopifnot(isTRUE(all.equal(res_reference$Fadv,res_off$Fadv,tolerance=1e-10)))
stopifnot(isTRUE(all.equal(res_reference$TAC,res_off$TAC,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(res_reference$SSB,res_off$SSB,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(res_reference$Catch,res_off$Catch,tolerance=1e-8)))

cat("\nPASS: probability OFF reproduces reference closed loop\n")

stopifnot(isTRUE(all.equal(res_reference$Fadv,res_CV0$Fadv,tolerance=1e-10)))
stopifnot(isTRUE(all.equal(res_reference$TAC,res_CV0$TAC,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(res_reference$SSB,res_CV0$SSB,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(res_reference$Catch,res_CV0$Catch,tolerance=1e-8)))

cat("PASS: probability CV=0 reproduces reference closed loop\n")

# ============================================================
# TEST 4 - S1: CURRENT-ASSESSMENT FORECAST UNCERTAINTY
#
# Uses the first-year forecast SSB CV estimated from the current
# base assessment.
#
# CV = 0.2395548
#
# This test verifies the operational response of the HCR when
# the escapement constraint is evaluated probabilistically.
# ============================================================
res_S1 <- run_scenario(
  label="S1 - CURRENT-ASSESSMENT FORECAST UNCERTAINTY",
  hcr_file=probability_file,
  use_probability=TRUE,
  forecast_CV=CV_S1
)

# ============================================================
# TEST 5 - S2: CONTEMPORARY FORECAST UNCERTAINTY
#
# Uses the median first-year forecast SSB CV obtained from the
# contemporary bootstrap assessment ensemble.
#
# CV = 0.29615
#
# For the conditioned trajectory used in this validation,
# Fcap does not satisfy the probabilistic escapement criterion,
# while F = 0 does.
#
# Therefore, the HCR is expected to select an intermediate F
# such that:
#
#   P(SSBforecast < Besc) = alpha
#
# equivalently:
#
#   Qalpha(SSBforecast) = Besc
#
# This expected branch is specific to the validation trajectory.
# ============================================================
res_S2 <- run_scenario(
  label="S2 - CONTEMPORARY MEDIAN FORECAST UNCERTAINTY",
  hcr_file=probability_file,
  use_probability=TRUE,
  forecast_CV=CV_S2,
  diagnostics=FALSE
)
# ============================================================
# FORECAST-UNCERTAINTY RESPONSE CHECKS
#
# These checks apply only to the conditioned trajectory used
# in this operational validation.
#
# S1:
#   the probabilistic escapement criterion remains compatible
#   with Fcap.
#
# S2:
#   Fcap does not satisfy the probabilistic escapement criterion,
#   whereas F = 0 does. Therefore, the HCR selects an intermediate
#   fishing mortality within (0,Fcap).
# ============================================================
stopifnot(isTRUE(all.equal(res_S1$Fadv,Fcap,tolerance=1e-10)))
stopifnot(res_S2$Fadv>0)
stopifnot(res_S2$Fadv<Fcap)
stopifnot(res_S2$Fadv<res_S1$Fadv)

cat("\nPASS: S1 retains Fcap for the validation trajectory\n")
cat("PASS: S2 reduces F below Fcap for the validation trajectory\n")

# ============================================================
# SUMMARY
# ============================================================

summary_results <- data.frame(
  Scenario=c(
    "Reference",
    "Probability OFF",
    "Probability CV=0",
    "S1 current assessment",
    "S2 contemporary median"),
  CV=c(NA,NA,CV_S0,CV_S1,CV_S2),
  Fadv=c(res_reference$Fadv,res_off$Fadv,res_CV0$Fadv,res_S1$Fadv,res_S2$Fadv),
  TAC=c(res_reference$TAC,res_off$TAC,res_CV0$TAC,res_S1$TAC,res_S2$TAC),
  SSB=c(res_reference$SSB,res_off$SSB,res_CV0$SSB,res_S1$SSB,res_S2$SSB),
  Catch=c(res_reference$Catch,res_off$Catch,res_CV0$Catch,res_S1$Catch,res_S2$Catch)
)

print(summary_results,row.names=FALSE)

# ============================================================
# FINAL STRUCTURAL CHECKS
#
# Verify that all operational outputs are finite and remain
# within their expected physical and management bounds.
# ============================================================

stopifnot(all(is.finite(summary_results$Fadv)))
stopifnot(all(summary_results$Fadv>=0))
stopifnot(all(summary_results$Fadv<=Fcap+1e-8))

stopifnot(all(is.finite(summary_results$TAC)))
stopifnot(all(summary_results$TAC>=0))

stopifnot(all(is.finite(summary_results$SSB)))
stopifnot(all(summary_results$SSB>=0))

stopifnot(all(is.finite(summary_results$Catch)))
stopifnot(all(summary_results$Catch>=0))

cat("\nPASS: all operational outputs are finite and within bounds\n")

# ============================================================
# COMPLETION
# ============================================================

cat("\n============================================================\n")
cat("OPERATIONAL PROBABILISTIC HCR VALIDATION COMPLETED\n")
cat("============================================================\n")
cat("Trajectory:",test_iter,"\n")
cat("Projection:",last.obs.yr,"->",proj.yr,"\n")
cat("Assessment error: OFF\n")
cat("Observation error: OFF\n")
cat("Fcap:",Fcap,"\n")
cat("Besc:",Besc,"\n")
cat("Risk threshold alpha:",alpha,"\n")
cat("S1 forecast CV:",CV_S1,"\n")
cat("S2 forecast CV:",CV_S2,"\n")
cat("------------------------------------------------------------\n")
cat("Backward compatibility: PASS\n")
cat("Zero-CV equivalence: PASS\n")
cat("Operational uncertainty response: PASS\n")
cat("============================================================\n")

