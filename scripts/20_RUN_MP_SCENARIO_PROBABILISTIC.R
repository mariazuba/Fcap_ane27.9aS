# ============================================================
# 20_RUN_MP_SCENARIO_PROBABILISTIC.R
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
#   Rscript scripts/20_RUN_MP_SCENARIO_PROBABILISTIC.R <scenario_number> <block> [scenario_table.csv]
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
# COMMAND-LINE ARGUMENTS AND LOCAL TEST MODE
# ============================================================

test_mode <- interactive()

if(test_mode) {
  scenario_number <- 10
  block <- 1
  scenario_file <- default_scenario_file
}

if(!test_mode) {
  args <- commandArgs(trailingOnly=TRUE)
  if(!length(args) %in% c(2,3)) stop("Usage: Rscript 20_RUN_MP_SCENARIO_PROBABILISTIC.R <scenario_number> <block> [scenario_table.csv]")
  scenario_number <- as.integer(args[1])
  block <- as.integer(args[2])
  scenario_file <- if(length(args)==3) args[3] else default_scenario_file
}

if(is.na(scenario_number) || scenario_number<1) stop("Invalid scenario number")
if(is.na(block) || block<1) stop("Invalid block number")
if(!file.exists(scenario_file)) stop("Scenario table not found: ",scenario_file)


scenario_table <- read.csv(scenario_file,stringsAsFactors=FALSE)

required_cols <- c("scenario_number",
                   "scenario_id",
                   "Fcap",
                   "Besc",
                   "sigmaR_scenario",
                   "sigmaR",
                   "seasonal_exploitation",
                   "Fprop_Q1",
                   "Fprop_Q2",
                   "Fprop_Q3",
                   "Fprop_Q4",
                   "assessment_uncertainty",
                   "assessment_error",
                   "forecast_uncertainty",
                   "forecast_probability",
                   "forecast_alpha")

missing_cols <- setdiff(required_cols,names(scenario_table))
if(length(missing_cols)>0) stop("Missing columns in scenario table: ",paste(missing_cols,collapse=", "))
sc <- scenario_table %>% dplyr::filter(scenario_number==!!scenario_number)
if(nrow(sc)!=1) stop("Scenario number must identify exactly one row")

scenario_id <- sc$scenario_id
Fcap <- sc$Fcap
Besc <- sc$Besc
sigmaR_scenario <- sc$sigmaR_scenario
sigmaR <- sc$sigmaR
seasonal_exploitation <- sc$seasonal_exploitation
Fprop <- c(sc$Fprop_Q1,sc$Fprop_Q2,sc$Fprop_Q3,sc$Fprop_Q4)
assessment_uncertainty <- sc$assessment_uncertainty
assessment_error_ON <- as.logical(sc$assessment_error)
forecast_uncertainty <- sc$forecast_uncertainty
forecast_probability <- as.logical(sc$forecast_probability)
forecast_alpha <- sc$forecast_alpha

stopifnot(Fcap>0,Besc>0,sigmaR>0,all(Fprop>=0),
          abs(sum(Fprop)-1)<1e-8,!is.na(assessment_error_ON),
          !is.na(forecast_probability),forecast_alpha>0,forecast_alpha<1)

n_blocks <- ceiling(ni/block_size)
if(block>n_blocks) stop("Invalid block number. Must be between 1 and ",n_blocks)
out_dir <- file.path(here(),"data","mse","MP_runs",scenario_id)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

SR_file <- file.path(rds_boot_dir,"SRs_reference_1000iter_process_error.rds")
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
cat("sigmaR:",sigmaR,"\n")
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
if(test_mode) ii <- head(ii,10)
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
if(test_mode) stopifnot(nb==10)
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
# REFERENCE SEASONAL EXPLOITATION
#
# Stage 1 uses the OM-specific conditioned seasonal shares
# stored in fleets.ctrl.
#
# Fprop_Q1-Q4 from the scenario table describe the reference
# seasonal exploitation scenario and are retained as metadata;
# they do not overwrite the conditioned OM-specific shares.
# ============================================================

seasonal_check <- seasonSums(fleets.ctrl_b$seasonal.share[[1]][,ac(proj.yr:last.yr),,,,])

stopifnot(all(abs(as.numeric(seasonal_check)-1)<1e-8))

cat("Reference seasonal exploitation retained from fleets.ctrl\n")
cat("Seasonal exploitation scenario:",seasonal_exploitation,"\n")

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
cat("Reference recruitment process error: bootstrap-conditioned OM\n")
# ============================================================
# RUN FLBEIA
# ============================================================

cat("\n============================================\n")
cat("RUNNING CANDIDATE MP CLOSED-LOOP MSE\n")
cat("Fcap:",Fcap,"\n")
cat("Besc:",Besc,"t\n")
cat("Block:",block,"/",n_blocks,"\n")
cat("============================================\n")

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

run_metadata <- data.frame(scenario_number=scenario_number,
                           scenario_id=scenario_id,
                           Fcap=Fcap,
                           Besc=Besc,
                           Blim=Blim,
                           Bpa=Bpa,
                           sigmaR_scenario=sigmaR_scenario,
                           sigmaR=sigmaR,
                           seasonal_exploitation=seasonal_exploitation,
                           Fprop_Q1=Fprop[1],
                           Fprop_Q2=Fprop[2],
                           Fprop_Q3=Fprop[3],
                           Fprop_Q4=Fprop[4],
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



