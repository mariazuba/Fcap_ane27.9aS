# ============================================================
# 15_VALIDATE_OBSERVATION_PROJECTION_INFRASTRUCTURE.R
# Observation and projection infrastructure validation
# Shortcut MSE - Anchovy ane.27.9aS
#
# Purpose:
#   Validate the technical infrastructure required before the
#   Management Procedure is introduced into the closed loop.
#
# Scope:
#   1. Check that perfectObs4seas() returns the expected seasonal
#      FLStock structure and exposes the realised fishing mortality
#      used by the Management Procedure through the harvest slot.
#
#   2. Check recruitment timing through FLasher::stf() and
#      FLasher::fwd(), including the Q3 recruitment implementation.
#
# Important:
#   - This script does NOT validate the complete closed loop.
#   - It does NOT apply assessment error to the perceived stock.
#   - It does NOT evaluate the Fcap + Besc Management Procedure.
#   - True/perceived-stock separation and assessment-error
#     propagation into management advice are deliberately deferred
#     until Fcap + Besc has been implemented.
#
# Workflow:
#
#   Validated reference OM
#           |
#           v
#   perfectObs4seas()
#           |
#           +----> realised-F observation check
#           |
#           v
#   stf() / fwd()
#           |
#           +----> recruitment-timing check
#           |
#           v
#   Observation/projection infrastructure validated
#           |
#           v
#   Assessment-error model (script 14)
#           |
#           v
#   Fcap + Besc MP implementation
#           |
#           v
#   Complete closed-loop validation
#
# Required existing project objects/functions:
#   functions/perfectObs4seas.R
#   OM_validation_C7000_block_01.rds
#   fit_bh_boot_unconstrained_100.rds
#

# NOTE:
#   The original tests are preserved below. The three project
#   dependencies above were implicit in the original code and are
#   therefore NOT guessed or silently reconstructed here.
#   Define them using the same validated objects/helpers already
#   used in the OM-validation workflow before running Section 4.
# ============================================================

rm(list=ls())

library(FLCore)
library(FLFishery)
library(FLasher)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(MASS)
library(callr)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

test_year <- 2025L

rds_dir <- here("data","Rdata")
dir_C7000 <- here("data","mse","OM_validation_C7000")

# ============================================================
# 2. LOAD VALIDATED INPUTS AND FUNCTIONS
# ============================================================

source(here("functions","perfectObs4seas.R"))

OM_test <- readRDS(
  file.path(
    dir_C7000,
    "OM_validation_C7000_block_01.rds"
  )
)

stk_base <- readRDS(
  file.path(
    rds_dir,
    "ane_stock_conditioned_1000iter.rds"
  )
)

bh_unconstrained <- readRDS(
  file.path(
    rds_dir,
    "fit_bh_boot_unconstrained_100.rds"
  )
)

sr_test <- bh_unconstrained[[1]]

# ============================================================
# 3. BASIC INPUT CHECKS
# ============================================================

stopifnot(inherits(OM_test$biols$ANE,"FLBiol"))
stopifnot(length(OM_test$fleets)>0)
stopifnot(inherits(stk_base,"FLStock"))
stopifnot(dims(stk_base)$season==4)

range(stk_base,"minfbar") <- 3
range(stk_base,"maxfbar") <- 3

validObject(stk_base)

n_iter_OM <- dim(OM_test$biols$ANE@n)[6]
n_iter_stock <- dims(stk_base)$iter

cat("OM test iterations =",n_iter_OM,"\n")
cat("Historical stock iterations =",n_iter_stock,"\n")
# ============================================================
# 3.1 MATCH CONDITIONED STOCK TO C7000 BLOCK
# ============================================================

map_C7000 <- read.csv(
  file.path(
    dir_C7000,
    "OM_validation_C7000_block_01_map.csv"
  )
)

stopifnot("iter1000" %in% names(map_C7000))
stopifnot(nrow(map_C7000)==n_iter_OM)

ii <- map_C7000$iter1000

stopifnot(length(ii)==n_iter_OM)
stopifnot(all(ii>=1))
stopifnot(all(ii<=n_iter_stock))
stopifnot(!anyDuplicated(ii))

sub_flq <- function(x,ii) {
  
  jj <- if(dim(x)[6]==1) {
    rep(1,length(ii))
  } else {
    ii
  }
  
  out <- x[,,,,,jj,drop=FALSE]
  
  dimnames(out)$iter <- as.character(seq_along(ii))
  
  out
}

stk_base_block <- stk_base

for(s in slotNames(stk_base)) {
  
  z <- slot(stk_base,s)
  
  if(inherits(z,"FLQuant")) {
    slot(stk_base_block,s,check=FALSE) <- sub_flq(z,ii)
  }
}

range(stk_base_block,"minfbar") <- 3
range(stk_base_block,"maxfbar") <- 3

validObject(stk_base_block)

cat(
  "C7000 block iterations =",
  n_iter_OM,
  "\n"
)

cat(
  "Matched stock iterations =",
  dims(stk_base_block)$iter,
  "\n"
)

cat(
  "Global iterations =",
  min(ii),
  "-",
  max(ii),
  "\n"
)

stopifnot(dims(stk_base_block)$iter==n_iter_OM)
# ============================================================
# 4. VALIDATE PERFECT OBSERVATION: REALISED F
# ============================================================
stk_obs_test <- perfectObs4seas(
  biol=OM_test$biols$ANE,
  fleets=OM_test$fleets,
  covars=NULL,
  obs.ctrl=NULL,
  year=test_year,
  stknm="ANE",
  historical.stock=stk_base_block
)
# ============================================================
# 4.1 CHECK FLSTOCK STRUCTURE
# ============================================================

perfectObs_structure <- tibble(
  class=class(stk_obs_test)[1],
  n_age=dim(stk_obs_test)[1],
  n_year=dim(stk_obs_test)[2],
  n_season=dim(stk_obs_test)[4],
  n_iter=dim(stk_obs_test)[6],
  first_year=min(as.numeric(dimnames(stk_obs_test)$year)),
  last_year=max(as.numeric(dimnames(stk_obs_test)$year))
)

print(perfectObs_structure)

stopifnot(inherits(stk_obs_test,"FLStock"))
stopifnot(dims(stk_obs_test)$season==4)

# ============================================================
# 4.2 CHECK REALISED F EXPOSED BY perfectObs4seas()
# ============================================================

F_obs_test <- harvest(stk_obs_test)[,ac(test_year),,,,]

stopifnot(all(is.finite(as.numeric(F_obs_test))))
stopifnot(all(as.numeric(F_obs_test)>=0))

F_real_test <- sapply(
  seq_len(dims(stk_obs_test)$iter),
  function(i) {
    mean(
      as.numeric(
        harvest(stk_obs_test)["3",ac(test_year),,,,i]
      )
    )
  }
)

stopifnot(all(is.finite(F_real_test)))
stopifnot(all(F_real_test>=0))

perfectObs_F_summary <- tibble(
  year=test_year,
  n_iter=length(F_real_test),
  median_F=median(F_real_test),
  min_F=min(F_real_test),
  max_F=max(F_real_test)
)

print(perfectObs_F_summary)

# ============================================================

# 5. VALIDATE RECRUITMENT TIMING: stf() + FLasher::fwd()
# ============================================================

bh_unconstrained <- readRDS(here("data","Rdata","fit_bh_boot_unconstrained_100.rds"))
sr_test <- bh_unconstrained[[1]]

stk_stf_test <- FLasher::stf(
  stk_obs_test,
  nyears=1,
  wts.nyears=3,
  fbar.nyears=3,
  f.rescale=TRUE
)

rbind(
  `2024`=as.numeric(stock.n(stk_stf_test)["0","2024",,,,1]),
  `2025`=as.numeric(stock.n(stk_stf_test)["0","2025",,,,1])
)

stk_stf_1 <- iter(stk_stf_test,1)

df_rec_test <- data.frame(
  year=rep(2025,4),
  quant=rep("f",4),
  value=rep(0,4),
  season=1:4
)

ctrl_rec_test <- FLasher::fwdControl(df_rec_test,iters=1)


run_flasher_fwd <- function(stk,ctrl,sr,maxF=1000) {
  callr::r(
    func=function(stk,ctrl,sr,maxF) {
      library(methods)
      library(FLCore)
      library(FLFishery)
      library(FLasher)
      FLasher::fwd(stk,control=ctrl,sr=sr,maxF=maxF)
    },
    args=list(stk=stk,ctrl=ctrl,sr=sr,maxF=maxF)
  )
}


stk_fwd_test <- run_flasher_fwd(
  stk=stk_stf_1,
  ctrl=ctrl_rec_test,
  sr=sr_test,
  maxF=1000
)

rbind(
  before_fwd=as.numeric(stock.n(stk_stf_1)["0","2025",,,,1]),
  after_fwd=as.numeric(stock.n(stk_fwd_test)["0","2025",,,,1])
)

m_fwd <- getMethod(
  "fwd",
  signature=signature(
    object="FLStock",
    fishery="missing",
    control="fwdControl"
  ),
  where=getNamespace("FLasher")
)

m_fwd

B_test <- as(stk_stf_1,"FLBiol")

sr_pm <- predictModel(
  model=model(sr_test),
  params=params(sr_test)
)

class(sr_pm)
sr_pm

rec(B_test) <- sr_pm

B_test@rec

sr_test_q3 <- rec(stk_stf_1)[,,,3]

stk_fwd_test_q3 <- run_flasher_fwd(
  stk=stk_stf_1,
  ctrl=ctrl_rec_test,
  sr=sr_test_q3,
  maxF=1000)

rbind(
  before_fwd=as.numeric(stock.n(stk_stf_1)["0","2025",,,,1]),
  after_fwd_q3=as.numeric(stock.n(stk_fwd_test_q3)["0","2025",,,,1])
)

# ============================================================
# 6. VALIDATION SUMMARY
# ============================================================

observation_projection_validation <- list(
  perfectObs_structure=perfectObs_structure,
  perfectObs_F_summary=perfectObs_F_summary
)

print(perfectObs_structure)
print(perfectObs_F_summary)

message(
  "Observation/projection infrastructure checks completed. ",
  "Assessment-error propagation and complete closed-loop validation ",
  "remain deferred until implementation of the Fcap + Besc MP."
)

