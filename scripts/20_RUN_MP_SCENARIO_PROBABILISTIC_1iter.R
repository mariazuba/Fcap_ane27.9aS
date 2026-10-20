# ============================================================
# 20_TEST_FCAP_BESC_DETERMINISTIC.R
#
# Deterministic diagnostic of the Fcap-Besc HCR
#
# Historical period: 1989-2024
# Projection period: 2025-2054
#
# Objective:
# Check whether the HCR maintains the escapement biomass
# and how the operating model responds to high fishing mortality.
#
# Configuration:
# - One OM and one iteration
# - Perfect observation
# - No assessment error
# - No recruitment process error during projections
# - No probabilistic forecast
#
# Only Fcap and Besc are changed between tests.
# ============================================================

# Clear the R environment and set a seed for reproducibility.
rm(list = ls())
set.seed(12345)

# Load packages required for FLBEIA simulations.
library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)

# ============================================================
# 1. DIRECTORIES
# ============================================================

# Define the directory containing the conditioned OM inputs.
rds_boot_dir <- file.path(here(), "data", "Rdata")

# ============================================================
# 2. GENERAL PARAMETERS
# ============================================================

# Define the historical and projection periods.
first.yr <- 1989
last.obs.yr <- 2024
proj.yr <- 2025
last.yr <- 2036 #2054

# Define the number of iterations in the original input files.
# Only one iteration will be selected for the simulation.
ni <- 1000

# Define the number of seasons and stock name.
ns <- 4
stks <- "ANE"

# Define the biological reference points (tonnes).
Blim <- 4721
Bpa <- 6561

# Number of years used to calculate reference fishing mortality.
nFyears <- 5

# ============================================================
# 3. SELECT OM AND HCR PARAMETERS
# ============================================================

# Select global iteration 1, corresponding to the reference OM 1.
# Only one iteration will be simulated.
ii <- 1L
nb <- 1L

# Define the maximum fishing mortality allowed by the HCR.
# Increase Fcap in later tests to examine the HCR response.
Fcap <- 2

# Define the escapement biomass.
# Here Besc is equal to Blim.
Besc <- Blim

# Disable assessment error in SSB, F and recruitment.
assessment_error_ON <- FALSE

# Disable the probabilistic forecast option.
# The HCR will use its non-probabilistic configuration.
forecast_probability <- FALSE

# Retain alpha for compatibility with the HCR interface.
# It should not affect the result when probability is disabled.
forecast_alpha <- 0.05

# ============================================================
# 4. LOAD CUSTOM FUNCTIONS
# ============================================================

# Load the observation function used in the simulations.
source(file.path(here(), "functions", "perfectObs4seas.R"))

# Load the Fcap-Besc HCR function.
source(file.path(here(), "functions", "FcapBpaHCR_ane9aS_PROBABILITY_OPTION_Fquarter.R"))

# Set the function environment to the FLBEIA namespace.
environment(FcapBpaHCR_ane) <- asNamespace("FLBEIA")

# Register the custom HCR as the annualTAC function in FLBEIA.
assignInNamespace("annualTAC", FcapBpaHCR_ane, ns = "FLBEIA")

# ============================================================
# 5. SELECT SINGLE ITERATION
# ============================================================

# Load the map linking global iterations to OMs and replicates.
iteration_map <- read.csv(file.path(rds_boot_dir, "iteration_map_1000.csv"))

# Select the global iteration defined by ii.
run_map <- iteration_map[iteration_map$iter1000 == ii, , drop = FALSE]

# Assign local iteration number 1.
run_map$iter_local <- 1L

# Check that exactly one iteration was selected.
stopifnot(nrow(run_map) == 1)

# Retain a positive CV only for compatibility with the HCR.
# Forecast uncertainty is not intended to be used in this test.
run_map$forecast_CV <- 0

# ============================================================
# 6. OUTPUT DIRECTORY
# ============================================================

# Create a unique name identifying the OM and HCR parameters.
# scenario_id <- paste0("TEST_OM", run_map$om, "_Fcap", Fcap, "_Besc", Besc)
scenario_id <- paste0("TEST_OM", run_map$om, "_Mconstant_Fcap", Fcap)
# Define the directory where simulation results will be saved.
out_dir <- file.path(here(), "data", "mse", "MP_runs", scenario_id)

# Create the output directory if it does not exist.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 7. LOAD CONDITIONED OPERATING MODEL
# ============================================================

# Load biological inputs from the SS3 bootstrap OMs.
biols <- readRDS(file.path(rds_boot_dir, "biols_Mconstant_1000iter.rds"))

# Load stock-recruitment relationships and process-error multipliers.
SRs <- readRDS(file.path(rds_boot_dir, "SRs_reference_1000iter_process_error.rds"))

# Load the conditioned fishing fleets.
fleets <- readRDS(file.path(rds_boot_dir, "fleets_conditioned_1000iter.rds"))

# Load fleet controls, including seasonal fishing patterns.
fleets.ctrl <- readRDS(file.path(rds_boot_dir, "fleets_ctrl_reference_1000iter.rds"))

# Load assessment-error multipliers.
assessment_error <- readRDS(file.path(rds_boot_dir, "assessment_error_reference_1000iter.rds"))

# ============================================================
# 8. DISABLE ASSESSMENT ERROR
# ============================================================

# Set all assessment-error multipliers to 1.
# Estimated SSB, F and recruitment are not modified by assessment error.
assessment_error$multiplier[] <- 1

# Recruitment process error will be removed after subsetting SRs.
# Future recruitment multipliers will be set to 1 for 2025-2054.

# ============================================================
# 9. BASIC INPUT CHECKS
# ============================================================

# Check that the biological input contains 1000 iterations.
stopifnot(dim(biols$ANE@n)[6] == ni)

# Check that stock-recruitment parameters contain 1000 iterations.
stopifnot(dim(SRs$ANE@params)[4] == ni)

# Check that recruitment uncertainty contains 1000 iterations.
stopifnot(dim(SRs$ANE@uncertainty)[6] == ni)

# Check that fleet effort contains 1000 iterations.
stopifnot(dim(fleets$SEINE@effort)[6] == ni)

# Confirm that assessment error has been disabled.
stopifnot(all(assessment_error$multiplier == 1))

# ============================================================
# 10. TEST CONFIGURATION SUMMARY
# ============================================================

# Print the configuration before running FLBEIA.
cat("\n============================================\n")
cat("DETERMINISTIC HCR DIAGNOSTIC\n")
cat("============================================\n")
cat("OM:", run_map$om, "\n")
cat("Global iteration:", ii, "\n")
cat("Number of iterations:", nb, "\n")
cat("Fcap:", Fcap, "\n")
cat("Besc:", Besc, "t\n")
cat("Assessment error: OFF\n")
cat("Recruitment process error: TO BE DISABLED AFTER SUBSETTING\n")
cat("Probabilistic forecast: OFF\n")
cat("============================================\n")
# ============================================================
# 11. SUBSETTING FUNCTIONS
# ============================================================

# Select one iteration from an FLQuant object.
# If the object has only one iteration, keep that iteration.
# Rename the selected iteration as 1.

sub_flq <- function(x, ii) {
  jj <- if (dim(x)[6] == 1) 1L else ii
  out <- x[, , , , , jj, drop = FALSE]
  dimnames(out)$iter <- "1"
  return(out)
}

# Select one iteration from the biological model.
# Keep numbers-at-age, weights, maturity, natural mortality
# and other biological quantities from the selected OM.

sub_FLBiol <- function(x, ii) {
  out <- x
  for (s in slotNames(x)) {
    z <- slot(x, s)
    if (inherits(z, "FLQuant")) slot(out, s, check = FALSE) <- sub_flq(z, ii)
    if (inherits(z, "predictModel")) {
      z@.Data[[1]] <- sub_flq(z@.Data[[1]], ii)
      slot(out, s, check = FALSE) <- z
    }
  }
  validObject(out)
  return(out)
}

# Select one iteration from the stock-recruitment model.
# Keep the Beverton-Holt parameters of the selected OM.

sub_FLSR <- function(x, ii) {
  out <- x
  for (s in c("rec", "ssb", "uncertainty", "proportion")) slot(out, s, check = FALSE) <- sub_flq(slot(x, s), ii)
  pars <- x@params[, , , ii, drop = FALSE]
  dimnames(pars)$iter <- "1"
  slot(out, "params", check = FALSE) <- pars
  validObject(out)
  return(out)
}

# Select one iteration from the fishing fleet.
# Keep effort, fishing characteristics and catch information
# corresponding to the selected OM.

sub_fleet <- function(x, ii) {
  out <- x
  for (s in slotNames(x)) {
    z <- slot(x, s)
    if (inherits(z, "FLQuant")) slot(out, s, check = FALSE) <- sub_flq(z, ii)
  }
  met <- x@metiers$ALL
  for (s in slotNames(met)) {
    z <- slot(met, s)
    if (inherits(z, "FLQuant")) slot(met, s, check = FALSE) <- sub_flq(z, ii)
  }
  ca <- met@catches$ANE
  for (s in slotNames(ca)) {
    z <- slot(ca, s)
    if (inherits(z, "FLQuant")) slot(ca, s, check = FALSE) <- sub_flq(z, ii)
  }
  catches_b <- met@catches
  catches_b$ANE <- ca
  slot(met, "catches", check = FALSE) <- catches_b
  metiers_b <- out@metiers
  metiers_b$ALL <- met
  slot(out, "metiers", check = FALSE) <- metiers_b
  validObject(out)
  return(out)
}

# Select one iteration of assessment-error multipliers.
# These multipliers were previously set to 1.

sub_assessment_error <- function(x, ii) {
  out <- x
  out$multiplier <- x$multiplier[, ii, , drop = FALSE]
  dimnames(out$multiplier)$iter <- "1"
  return(out)
}

# ============================================================
# 12. SELECT ONE ITERATION OF THE OPERATING MODEL
# ============================================================

# Select biological inputs for the chosen OM.
biols_b <- biols
biols_b$ANE <- sub_FLBiol(biols$ANE, ii)
M.forecast <- m(biols_b$ANE)[, ac(proj.yr), , , , , drop=FALSE]

# Select the stock-recruitment relationship for the same OM.
SRs_b <- SRs
SRs_b$ANE <- sub_FLSR(SRs$ANE, ii)

# Remove recruitment process error during projections.
# A multiplier of 1 means deterministic recruitment.
# Beverton-Holt parameters remain unchanged.
SRs_b$ANE@uncertainty[, ac(proj.yr:last.yr), , , , ] <- 1

# Select fishing fleet information for the same OM.
fleets_b <- fleets
fleets_b$SEINE <- sub_fleet(fleets$SEINE, ii)

# Keep the original seasonal fishing pattern of the OM.
fleets.ctrl_b <- fleets.ctrl
fleets.ctrl_b$seasonal.share[[1]] <- sub_flq(fleets.ctrl$seasonal.share[[1]], ii)

# Select assessment-error multipliers for the same iteration.
# All multipliers are equal to 1.
assessment_error_b <- sub_assessment_error(assessment_error, ii)

# ============================================================
# 13. CHECK DETERMINISTIC OM
# ============================================================

# Confirm that all model components contain one iteration.
stopifnot(dim(biols_b$ANE@n)[6] == 1)
stopifnot(dim(SRs_b$ANE@params)[4] == 1)
stopifnot(dim(SRs_b$ANE@uncertainty)[6] == 1)
stopifnot(dim(fleets_b$SEINE@effort)[6] == 1)
stopifnot(dim(catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE))[6] == 1)
stopifnot(dim(fleets.ctrl_b$seasonal.share[[1]])[6] == 1)

# Confirm that recruitment process error is disabled.
stopifnot(all(SRs_b$ANE@uncertainty[, ac(proj.yr:last.yr), , , , ] == 1))

# Confirm that assessment error is disabled.
stopifnot(all(assessment_error_b$multiplier == 1))

cat("OM:", run_map$om, "\n")
cat("Number of iterations:", nb, "\n")
cat("Recruitment process error: OFF\n")
cat("Assessment error: OFF\n")
cat("Single OM configuration: OK\n")

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

# Keep the original seasonal fishing pattern of the selected OM.
# No changes are made to the seasonal exploitation proportions.

# Calculate the sum of seasonal proportions for each projection year.
seasonal_check <- seasonSums(fleets.ctrl_b$seasonal.share[[1]][, ac(proj.yr:last.yr), , , , ])

# Check that the four seasonal proportions sum to 1.
stopifnot(all(abs(as.numeric(seasonal_check) - 1) < 1e-8))

cat("Original seasonal exploitation pattern retained: OK\n")

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
                                 M.forecast=M.forecast,
                                 catch.q=catch.q(fleets_b$SEINE@metiers$ALL@catches$ANE),
                                 seasonal.share=fleets.ctrl_b$seasonal.share$ANE,
                                 assessment.error=assessment_error_b,
                                 assess_error_R=FALSE,
                                 use.probability=forecast_probability,
                                 forecast.CV=run_map$forecast_CV,
                                 alpha=forecast_alpha,
                                 diagnostics=TRUE))

stopifnot(length(advice_HCR.ctrl$ANE$forecast.CV)==nb)
stopifnot(identical(as.numeric(advice_HCR.ctrl$ANE$forecast.CV),as.numeric(run_map$forecast_CV)))

cat("Forecast CV pairing with HCR control: OK\n")

stopifnot(
  isTRUE(all.equal(as.numeric(advice_HCR.ctrl$ANE$Fcap), as.numeric(Fcap))),
  isTRUE(all.equal(as.numeric(advice_HCR.ctrl$ANE$Besc), as.numeric(Besc)))
)

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

# stopifnot(Fcap == 20)
stopifnot(Besc == Blim)
stopifnot(assessment_error_ON == FALSE)
stopifnot(forecast_probability == FALSE)
stopifnot(all(assessment_error_b$multiplier == 1))

# ============================================================
# RUN FLBEIA
# ============================================================

# Run one deterministic trajectory of the selected OM.
cat("\n============================================\n")
cat("RUNNING DETERMINISTIC HCR TEST\n")
cat("OM:", run_map$om, "\n")
cat("Global iteration:", ii, "\n")
cat("Fcap:", Fcap, "\n")
cat("Besc:", Besc, "t\n")
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


# Calculate simulation time.
end_time <- Sys.time()
cat("Runtime:", round(as.numeric(difftime(end_time, start_time, units="mins")), 2), "minutes\n")

# ============================================================
# SAVE SIMULATION RESULTS
# ============================================================

saveRDS(OM_b, file.path(out_dir, paste0("OM_Mconstant_Fcap_",Fcap,".rds")))

# ============================================================
# DIAGNOSTICS - CONSTANT F TEST
# ============================================================

N <- OM_b$biols[[1]]@n
W <- OM_b$biols[[1]]@wt
M <- OM_b$biols[[1]]@m
F <- OM_b$stocks[[1]]@harvest

stk <- OM_b$stocks[[1]]

years <- as.character(proj.yr:last.yr)
ages <- dimnames(N)$age





# SSB
plot(ssb(OM_b$stocks[[1]])[, ac(1989:2030), , "2"]) +
  geom_flpar(data = FLPar(Blim = 4721, Bpa = 6561), x = 2025)+
  facet_wrap(~season) +
  theme_bw()

# rec
plot(rec(OM_b$stocks[[1]])[, ac(1989:2030), , "3"]) +
  facet_wrap(~season) +
  theme_bw()
# catch
plot(catch(OM_b$stocks[[1]])[, ac(1989:2030), , ]) +
  facet_wrap(~season) +
  theme_bw()
#fbar
plot(fbar(OM_b$stocks[[1]])[, ac(1989:2030), , ]) +
  facet_wrap(~season) +
  theme_bw()

# N-at-age
ggplot(OM_b$biols[[1]]@n[, ac(2020:2030)],aes(x = year, y = data,group=age, color = age)) +
  geom_line(linewidth = 1) +
  facet_wrap(~season, nrow = 2) +
  labs(title = "Abundance at age by quarter",x = "Year",y = "Abundance (thousands)", color = "Age") +
  theme_bw()

# harvest
ggplot(OM_b$stocks[[1]]@harvest[, ac(1989:2030)],aes(x = year, y = data, group=age,colour = age)) +
  geom_line() +facet_wrap(~season) +theme_bw()

# natural mortality
ggplot(OM_b$biols[[1]]@m[, ac(1989:2030)],aes(x = year, y = data, group=age,colour = age)) +
  geom_line() +facet_wrap(~season) +theme_bw()



#variables guardadas de la HCR
HCR_results <- lapply(list.files("HCR_diagnostics", pattern="\\.rds$", full.names=TRUE),readRDS)
names(HCR_results) <- sapply(HCR_results, function(x) as.character(x$year))

N_HCR <- simplify2array(lapply(HCR_results, `[[`, "N"))
W_HCR <- simplify2array(lapply(HCR_results, `[[`, "W"))
Mat_HCR <- simplify2array(lapply(HCR_results, `[[`, "Mat"))
F_HCR <- simplify2array(lapply(HCR_results, `[[`, "F"))
M_HCR <- simplify2array(lapply(HCR_results, `[[`, "M"))


# HCR
N_HCR_df <- as.data.frame.table(N_HCR, responseName="N")
names(N_HCR_df) <- c("age", "season", "year", "N")
N_HCR_df$Model <- "HCR"

# OM
N_OM_df <- as.data.frame(N[, ac(2025:2030)])
N_OM_df <- N_OM_df |>
  select(age, season, year, N=data) |>
  mutate(Model="OM")

# Combinar
N_comparison <- bind_rows(
  N_HCR_df |> mutate(across(c(age, season, year), as.character)),
  N_OM_df |> mutate(across(c(age, season, year), as.character))) |>
  mutate(
    year = as.numeric(year),
    age = factor(age, levels = c("0", "1", "2", "3")),
    season = factor(season, levels = as.character(1:4),
                    labels = paste0("Q", 1:4)))

# Figura
ggplot(N_comparison, aes(x=year, y=N, color=Model)) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  facet_grid(age ~ season, scales="free_y") +
  scale_x_continuous(breaks=2025:2030) +
  labs(title="Numbers-at-age: HCR vs OM",x="Year",y="Abundance",color="Model") +
  theme_bw() +
  theme(strip.text=element_text(size=12),axis.text=element_text(size=10),axis.text.x=element_text(angle=45, hjust=1),legend.position="bottom")


# Natural mortality
M_OM <- OM_b$biols[[1]]@m

M_HCR_df <- as.data.frame.table(M_HCR, responseName="M")
names(M_HCR_df) <- c("age", "season", "year", "M")
M_HCR_df$Model <- "HCR"

M_OM_df <- as.data.frame(M_OM[, ac(2025:2030)]) |>
  dplyr::select(age, season, year, M=data) |>
  dplyr::mutate(Model="OM")

M_comparison <- dplyr::bind_rows(
  M_HCR_df |> dplyr::mutate(dplyr::across(c(age, season, year), as.character)),
  M_OM_df |> dplyr::mutate(dplyr::across(c(age, season, year), as.character))) |>
  dplyr::mutate(
    year=as.numeric(year),
    age=factor(age, levels=c("0","1","2","3")),
    season=factor(season, levels=as.character(1:4), labels=paste0("Q",1:4)))


ggplot(M_comparison, aes(x=year, y=M, color=Model)) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  facet_grid(age ~ season, scales="free_y") +
  scale_x_continuous(breaks=2025:2030) +
  labs(
    title="Natural mortality-at-age: HCR vs OM",
    x="Year",
    y="Natural mortality (M)",
    color="Model"
  ) +
  theme_bw() +
  theme(
    strip.text=element_text(size=12),
    axis.text.x=element_text(angle=45, hjust=1),
    legend.position="bottom"
  )





# Reclutamiento

R_HCR <- sapply(HCR_results, `[[`, "R")

R_OM <- as.numeric(
  OM_b$biols[[1]]@n["0", ac(2025:2030), , "3"]
)

R_comparison <- data.frame(
  Year=2025:2030,
  HCR=as.numeric(R_HCR),
  OM=R_OM
)

R_comparison$Diff_pct <- 100 *
  (R_comparison$HCR - R_comparison$OM) /
  ifelse(R_comparison$OM == 0, NA, R_comparison$OM)

print(R_comparison)



ggplot(R_comparison, aes(x=Year)) +
  geom_line(aes(y=HCR, color="HCR"), linewidth=1) +
  geom_line(aes(y=OM, color="OM"), linewidth=1) +
  geom_point(aes(y=HCR, color="HCR"), size=2) +
  geom_point(aes(y=OM, color="OM"), size=2) +
  scale_x_continuous(breaks=2025:2030) +
  labs(
    title="Recruitment Q3: HCR vs OM",
    x="Year",
    y="Recruitment",
    color="Model"
  ) +
  theme_bw()


# SSB 
SSB_HCR <- sapply(HCR_results, `[[`, "SSB")

SSB_OM <- as.numeric(
  ssb(OM_b$stocks[[1]])[, ac(2025:2030), , "2"]
)

SSB_comparison <- data.frame(
  Year=2025:2030,
  HCR=as.numeric(SSB_HCR),
  OM=SSB_OM
)

SSB_comparison$Diff_pct <- 100 *
  (SSB_comparison$HCR - SSB_comparison$OM) /
  SSB_comparison$OM

print(SSB_comparison)

ggplot(SSB_comparison, aes(x=Year)) +
  geom_line(aes(y=HCR, color="HCR"), linewidth=1) +
  geom_line(aes(y=OM, color="OM"), linewidth=1) +
  geom_point(aes(y=HCR, color="HCR"), size=2) +
  geom_point(aes(y=OM, color="OM"), size=2) +
  scale_x_continuous(breaks=2025:2030) +
  labs(
    title="Spawning stock biomass Q2: HCR vs OM",
    x="Year",
    y="SSB (tonnes)",
    color="Model"
  ) +
  theme_bw()

# ============================================================
# SSB: HCR vs OM (FLBiol)
# ============================================================

biol <- OM_b$biols[[1]]

SSB_HCR <- sapply(HCR_results, `[[`, "SSB")

SSB_OM <- sapply(2025:2030, function(y) {
  sum(
    biol@n[, ac(y), , "2"] *
      wt(biol)[, ac(y), , "2"] *
      mat(biol)[, ac(y), , "2"] *
      fec(biol)[, ac(y), , "2"]
  )
})

SSB_comparison <- data.frame(
  Year=2025:2030,
  HCR=as.numeric(SSB_HCR[as.character(2025:2030)]),
  OM=as.numeric(SSB_OM)
)

SSB_comparison$Diff_pct <- 100 *
  (SSB_comparison$HCR - SSB_comparison$OM) /
  SSB_comparison$OM

print(SSB_comparison)

HCR_results <- lapply(
  file.path("HCR_diagnostics",paste0("HCR_projection_",2025:2030,".rds")),
  readRDS
)

names(HCR_results) <- as.character(2025:2030)
# ============================================================
# QUARTERLY F: HCR vs OM
# ============================================================

F_HCR <- simplify2array(lapply(HCR_results, `[[`, "F"))
F_OM <- harvest(OM_b$stocks[[1]])

F_comparison <- expand.grid(Quarter=1:4, Year=2025:2030)

F_comparison$F_HCR <- mapply(function(q,y) F_HCR["2",q,as.character(y)], F_comparison$Quarter, F_comparison$Year)

F_comparison$F_OM <- mapply(function(q,y) as.numeric(F_OM["2",as.character(y),,as.character(q)])/4, F_comparison$Quarter, F_comparison$Year)

F_comparison$Diff_pct <- 100*(F_comparison$F_HCR-F_comparison$F_OM)/ifelse(F_comparison$F_OM==0,NA,F_comparison$F_OM)

print(F_comparison, row.names=FALSE)

# ============================================================
# COMPARE BEVERTON-HOLT PARAMETERS: HCR vs OM
# ============================================================

BH_HCR <- advice_HCR.ctrl$ANE$sr@params
BH_OM <- OM_b$SRs$ANE@params

BH_comparison <- data.frame(
  Year = proj.yr:last.yr,
  a_HCR = as.numeric(BH_HCR["a", ac(proj.yr:last.yr), "3", 1]),
  a_OM = as.numeric(BH_OM["a", ac(proj.yr:last.yr), "3", 1]),
  b_HCR = as.numeric(BH_HCR["b", ac(proj.yr:last.yr), "3", 1]),
  b_OM = as.numeric(BH_OM["b", ac(proj.yr:last.yr), "3", 1])
)

BH_comparison$Diff_a <- BH_comparison$a_HCR - BH_comparison$a_OM
BH_comparison$Diff_b <- BH_comparison$b_HCR - BH_comparison$b_OM

print(BH_comparison)

stopifnot(all.equal(BH_comparison$a_HCR, BH_comparison$a_OM))
stopifnot(all.equal(BH_comparison$b_HCR, BH_comparison$b_OM))
