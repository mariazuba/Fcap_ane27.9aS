rm(list = ls())

library(FLBEIA)
library(FLCore)
library(FLFishery)
library(dplyr)
library(here)

wd <- here()
setwd(wd)

source("functions/perfectObs4seas.R")
source("functions/FcapBpaHCR_ane9aS.R")

experiment_name <- "test_100iter_10year"
load(file.path("data/mse", paste0(experiment_name, ".RData")))


# load("data/mse/SR_bh.RData")
load("data/mse/obs_perfect.RData")
load("data/mse/assess_none.RData")
load("data/mse/adv_fcap.RData")

scenarios <- read.csv("config/scenarios.csv")

args <- commandArgs(trailingOnly = TRUE)
# args <- 1

if (length(args) == 0) {
  stop("Debe indicar scenario_id. Ejemplo: Rscript scripts/model.R 1")
}

scenario_id <- as.integer(args[1])
sc <- scenarios[scenarios$scenario_id == scenario_id, ]



if (nrow(sc) != 1) {
  stop("scenario_id no encontrado o duplicado: ", scenario_id)
}

sr_file <- file.path("data/mse", paste0("SR_", sc$SR, ".RData"))

if (!file.exists(sr_file)) {
  stop("No existe el archivo SR: ", sr_file)
}

load(sr_file)

dir.create("outputs/mse/res", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/mse/summary", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/mse/logs", recursive = TRUE, showWarnings = FALSE)

Fcap_value <- as.numeric(sc$Fcap)
Besc_value <- as.numeric(sc$Besc)

message(
  "Running scenario ", scenario_id,
  " | SR = ", sc$SR,
  " | Fcap = ", Fcap_value,
  " | Besc = ", Besc_value,
  " | catch_prop = ", sc$catch_prop
)


advice.ctrl.i <- advice.ctrl
advice.ctrl.i$ANE$Fcap <- Fcap_value
advice.ctrl.i$ANE$Besc <- Besc_value

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

out_name <- sprintf(
  "sc%03d_SR_%s_Fcap_%s_Besc_%s",
  scenario_id,
  sc$SR,
  gsub("\\.", "p", format(Fcap_value, nsmall = 2)),
  round(Besc_value)
)

saveRDS(res_i,file = file.path("outputs/mse/res", paste0(out_name, ".rds")))

saveRDS(
  list(
    scenario = sc,
    SR = sc$SR,
    Fcap = Fcap_value,
    Bescapement = Besc_value,
    catch_prop = sc$catch_prop,
    runtime = t1 - t0,
    start_time = t0,
    end_time = t1,
    ni = ni,
    proj.yrs = proj.yrs,
    Blim = Blim,
    Bpa = Bpa
  ),
  file = file.path("outputs/mse/summary", paste0(out_name, "_metadata.rds"))
)

message("Finished scenario ", scenario_id)
