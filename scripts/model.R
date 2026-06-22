################################################################################
# Run one MSE scenario
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(FLBEIA)
  library(FLCore)
  library(FLFishery)
  library(dplyr)
  library(here)
})

wd <- here()
setwd(wd)

#===============================================================================
# Source functions
#===============================================================================

source("functions/perfectObs4seas.R")
source("functions/FcapBpaHCR_ane9aS.R")

#===============================================================================
# Experiment configuration
#===============================================================================

experiment_name <- "test_1000iter_30year"

base_file      <- file.path("data/mse", paste0(experiment_name, ".RData"))
obs_file       <- "data/mse/obs_perfect.RData"
assess_file    <- "data/mse/assess_none.RData"
advice_file    <- "data/mse/adv_fcap.RData"
scenarios_file <- "config/scenarios.csv"

required_files <- c(
  base_file,
  obs_file,
  assess_file,
  advice_file,
  scenarios_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

load(base_file)
load(obs_file)
load(assess_file)
load(advice_file)

scenarios <- read.csv(scenarios_file, stringsAsFactors = FALSE)

#===============================================================================
# Read scenario id
#===============================================================================
args        <- commandArgs(trailingOnly = TRUE);  if (length(args) == 0) {stop("scenario_id required. Example: Rscript scripts/model.R 1")}
scenario_id <- suppressWarnings(as.integer(args[1])); if (is.na(scenario_id)) {stop("scenario_id must be an integer. Received: ", args[1])}
sc          <- scenarios[scenarios$scenario_id == scenario_id, ]; if (nrow(sc) != 1) {stop("scenario_id not found or duplicated: ", scenario_id)}
#===============================================================================
# Load scenario-specific objects
#===============================================================================

sr_file <- file.path("data/mse", paste0("SR_", sc$SR, ".RData")); if (!file.exists(sr_file)) {stop("SR file not found: ", sr_file)}

load(sr_file)

# Seasonal catch allocation object
catch_prop_file <- file.path("data/mse/catch_prop", paste0("fleets_ctrl_", sc$catch_prop, ".RData")); if (!file.exists(catch_prop_file)) {stop("catch_prop file not found: ", catch_prop_file)}

load(catch_prop_file); if (!exists("fleets.ctrl.cp")) {stop("Object fleets.ctrl.cp not found in: ", catch_prop_file)}

fleets.ctrl <- fleets.ctrl.cp

#===============================================================================
# Output directories
#===============================================================================

dir.create("outputs/mse/res", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/mse/summary", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/mse/logs", recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# Scenario-specific advice
#===============================================================================

Fcap_value <- as.numeric(sc$Fcap)
Besc_value <- as.numeric(sc$Besc)

advice.ctrl.i <- advice.ctrl
advice.ctrl.i$ANE$Fcap <- Fcap_value
advice.ctrl.i$ANE$Besc <- Besc_value

# Update propf from selected fleets.ctrl
seasonal_prop <- as.numeric(fleets.ctrl$seasonal.share[[1]][, ac(proj.yr), , , , 1])
seasonal_prop <- seasonal_prop / sum(seasonal_prop)

advice.ctrl.i$ANE$propf <- seasonal_prop

#===============================================================================
# Output name
#===============================================================================

fmt_num <- function(x, digits = 2) {gsub("\\.", "p", format(round(as.numeric(x), digits), nsmall = digits))}

out_name <- sprintf("sc%03d_SR_%s_Fcap_%s_Besc_%s_CP_%s",
  scenario_id,
  sc$SR,
  fmt_num(Fcap_value, 2),
  round(Besc_value),
  sc$catch_prop
)

res_file  <- file.path("outputs/mse/res", paste0(out_name, ".rds"))
meta_file <- file.path("outputs/mse/summary", paste0(out_name, "_metadata.rds"))

# Avoid overwriting completed runs
if (file.exists(res_file) && file.exists(meta_file)) {
  message("Scenario already completed. Skipping: ", scenario_id)
  quit(save = "no", status = 0)
}

#===============================================================================
# Run MSE
#===============================================================================

message(
  "Running scenario ", scenario_id,
  " | SR = ", sc$SR,
  " | Fcap = ", Fcap_value,
  " | Besc = ", Besc_value,
  " | catch_prop = ", sc$catch_prop,
  " | ni = ", ni,
  " | proj_nyr = ", proj.nyr
)

message(
  "Seasonal propf = ",
  paste(round(seasonal_prop, 3), collapse = ", ")
)

t0 <- Sys.time()

res_i <- FLBEIA(
  biols       = biols,
  SRs         = SRs,
  BDs         = NULL,
  fleets      = fleets,
  covars      = NULL,
  indices     = NULL,
  advice      = advice_TACNA,
  main.ctrl   = main.ctrl,
  biols.ctrl  = biols.ctrl,
  fleets.ctrl = fleets.ctrl,
  covars.ctrl = NULL,
  obs.ctrl    = obs.ctrl,
  assess.ctrl = assess.ctrl,
  advice.ctrl = advice.ctrl.i
)

t1 <- Sys.time()

#===============================================================================
# Save outputs
#===============================================================================

saveRDS(
  res_i,
  file = res_file
)

saveRDS(
  list(
    scenario_id   = scenario_id,
    scenario      = sc,
    SR            = sc$SR,
    Fcap          = Fcap_value,
    Bescapement   = Besc_value,
    catch_prop    = sc$catch_prop,
    seasonal_prop = seasonal_prop,
    runtime       = t1 - t0,
    start_time    = t0,
    end_time      = t1,
    ni            = ni,
    proj.yrs      = proj.yrs,
    Blim          = Blim,
    Bpa           = Bpa,
    res_file      = res_file
  ),
  file = meta_file
)

message("Finished scenario ", scenario_id)
message("Runtime: ", round(as.numeric(difftime(t1, t0, units = "mins")), 2), " min")