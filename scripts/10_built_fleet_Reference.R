# ============================================================
# 10_built_fleet_Reference.R
# Fleet conditioning for 1000 bootstrap-conditioned OMs
#
# Purpose:
#   Build the fishery component of the operating model for the
#   1000 bootstrap-conditioned biological trajectories.
#
#   The script:
#     - Expands the historical FLStock from 100 to 1000 iterations.
#     - Builds the SEINE fleet and ALL metier.
#     - Transfers historical catches and weights from the
#       conditioned SS3 stock.
#     - Incorporates external fleet variables.
#     - Extends fleet quantities through 2025-2054.
#     - Configures the SMFB effort and Baranov catch models.
#     - Estimates catchability from the final historical years.
#     - Defines seasonal catch shares from historical catches.
#     - Checks consistency between fleet and biological weights.
#
# INPUTS:
#   data/Rdata/
#     - biols_conditioned_1000iter.rds
#     - ane_stock_conditioned_100iter.rds
#     - iteration_map_1000.csv
#
#   boot/data/
#     - seine_vcost.csv
#     - seine.met1_effshare.csv
#     - seine_effort.csv
#     - seine_capacity.csv
#     - seine_crewshare.csv
#     - seine_fcost.csv
#
# OUTPUTS:
#   data/Rdata/
#     - ane_stock_conditioned_1000iter.rds
#     - fleets_conditioned_1000iter.rds
#     - fleets_ctrl_reference_1000iter.rds
#
# Simulation period:
#   Historical: 1989-2024
#   Projection: 2025-2054
#
# Fishery configuration:
#   Fleet: SEINE
#   Metier: ALL
#   Stock: ANE
#   Effort model: SMFB
#   Catch model: Baranov
#   Effort restriction: catch
#
# Notes:
#   - Historical fleet catches are inherited from the
#     bootstrap-conditioned SS3 trajectories.
#   - Future catch weights and fleet variables are based on
#     recent historical means.
#   - Seasonal catch shares for the projection are based on
#     the mean of the final 10 historical years.
# ============================================================

rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)

# ============================================================
# DIRECTORIES
# ============================================================

rds_boot_dir <- file.path(here(), "data", "Rdata")

# External fleet data
seine_vcost_csv         <- file.path(here(), "boot", "data", "seine_vcost.csv")
seine.met1_effshare_csv <- file.path(here(), "boot", "data", "seine.met1_effshare.csv")
seine_effort_csv        <- file.path(here(), "boot", "data", "seine_effort.csv")
seine_capacity_csv      <- file.path(here(), "boot", "data", "seine_capacity.csv")
seine_crewshare_csv     <- file.path(here(), "boot", "data", "seine_crewshare.csv")
seine_fcost_csv         <- file.path(here(), "boot", "data", "seine_fcost.csv")

# ============================================================
# PARAMETERS
# ============================================================

first.yr    <- 1989
last.obs.yr <- 2024

proj.yr  <- 2025
proj.nyr <- 30
last.yr  <- proj.yr + proj.nyr - 1

hist.yrs <- first.yr:last.obs.yr
proj.yrs <- proj.yr:last.yr

ass.yr <- last.obs.yr

ni <- 1000
ns <- 4

stks <- "ANE"

ANE.age.min <- 0
ANE.age.max <- 3
ANE.unit    <- 1

# ============================================================
# LOAD CONDITIONED OBJECTS
# ============================================================

biols <- readRDS(file.path(rds_boot_dir,"biols_conditioned_1000iter.rds"))
ane.stock.100 <- readRDS(file.path(rds_boot_dir,"ane_stock_conditioned_100iter.rds"))
iteration_map_1000 <- read.csv(file.path(rds_boot_dir,"iteration_map_1000.csv"))
source_iter <- iteration_map_1000$om

# ============================================================
# EXPAND FLStock 100 -> 1000
# ============================================================

expand_FLQuant <- function(x, source_iter) {
  if(dim(x)[6] == 1) {
    out <- x[, , , , ,rep(1, length(source_iter)),drop = FALSE]
  } else {
    out <- x[, , , , ,source_iter,drop = FALSE]
  }
  dimnames(out)$iter <-as.character(seq_along(source_iter))
  out
}


expand_FLStock <- function(x, source_iter) {
  out <- x
  for(s in slotNames(x)) {
    z <- slot(x, s)
    if(inherits(z, "FLQuant")) {
      slot(out, s) <-expand_FLQuant(z,source_iter)
    }
  }
  out
}


ane.stock <- expand_FLStock(ane.stock.100,source_iter)
wt_expand_diff <- max(abs(as.numeric(stock.wt(ane.stock)-landings.wt(ane.stock))),na.rm=TRUE)
stopifnot(wt_expand_diff < 1e-10)
# ============================================================
# FLEET STRUCTURE
# ============================================================

fls <- "SEINE"
SEINE.mets <- "ALL"
SEINE.ALL.stks <- "ANE"

fle <- "SEINE"
met <- "ALL"

maxyear <- last.obs.yr

# 1. Capturas históricas
catch <- FLCatchExt(name = names(biols),
                    landings.n =landings.n(ane.stock[, ac(hist.yrs)]),
                    landings =landings(ane.stock[, ac(hist.yrs)]),
                    landings.wt =landings.wt(ane.stock[, ac(hist.yrs)]),
                    discards.n =discards.n(ane.stock[, ac(hist.yrs)]),
                    discards =discards(ane.stock[, ac(hist.yrs)]),
                    discards.wt =discards.wt(ane.stock[, ac(hist.yrs)]))

catches <- FLCatchesExt(catch)
names(catches) <- SEINE.ALL.stks

m <- FLMetierExt(catches = catches,name = SEINE.mets)

# 2. Variables externas de la flota
# Estas no vienen de SS3 

# Variable cost
seine.vcost <-read.csv(seine_vcost_csv) %>%filter(year <= maxyear)
seine.vcost$data <- 1
seine.vcost.flq <-propagate(as.FLQuant(seine.vcost),ni)

# Effort share
seine.met1_effshare <-read.csv(seine.met1_effshare_csv) %>%filter(year <= maxyear)
seine.met1_effshare.flq <-propagate(as.FLQuant(seine.met1_effshare),ni)

m@effshare <- seine.met1_effshare.flq
m@vcost    <- seine.vcost.flq
metiers <- FLMetiersExt(m)
names(metiers) <- SEINE.mets

# Componentes a nivel de flota

# Effort
seine.effort <-read.csv(seine_effort_csv) %>%filter(year <= maxyear)
seine.effort$data <- 1
seine.effort.flq <-propagate(as.FLQuant(seine.effort),ni)

# Capacity
seine.capacity <-read.csv(seine_capacity_csv) %>%filter(year <= maxyear)
seine.capacity$data <- 5000
seine.capacity.flq <-propagate(as.FLQuant(seine.capacity),ni)

# Crew share
seine.crewshare <-read.csv(seine_crewshare_csv) %>%filter(year <= maxyear)
seine.crewshare$data <- 0
seine.crewshare.flq <-propagate(as.FLQuant(seine.crewshare),ni)

# Fixed cost
seine.fcost <-read.csv(seine_fcost_csv) %>%filter(year <= maxyear)
seine.fcost$data <- 0
seine.fcost.flq <-propagate(as.FLQuant(seine.fcost),ni)

# 3. Construir FLFleetExt
fleet <- FLFleetExt(
  metiers   = metiers,
  name      = fls,
  effort    = seine.effort.flq,
  fcost     = seine.fcost.flq,
  capacity  = seine.capacity.flq,
  crewshare = seine.crewshare.flq)

fleets <- FLFleetsExt(fleet)
names(fleets) <- "SEINE"

cobj <- fleets$SEINE@metiers$ALL@catches$ANE
wt_fleet_hist_diff <- max(abs(as.numeric(landings.wt(cobj)[,ac(hist.yrs)]-wt(biols$ANE)[,ac(hist.yrs)])),na.rm=TRUE)
stopifnot(wt_fleet_hist_diff < 1e-10)

# 4. Unidades y selectividad landings/discards
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("landings.n", "discards.n")] <-units(biols[[stks]])$n
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("landings.wt", "discards.wt")] <-units(biols[[stks]])$wt
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("landings", "discards")] <-units(biols[[stks]]@n *wt(biols[[stks]]))
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("alpha", "beta", "catch.q")] <- "1"
fleets[[fle]]@metiers[[met]]@catches[[stks]]@landings.sel[] <- 1
fleets[[fle]]@metiers[[met]]@catches[[stks]]@discards.sel[] <- 0
fleets[[fle]]@metiers[[met]]@catches[[stks]]@discards.wt <-fleets[[fle]]@metiers[[met]]@catches[[stks]]@landings.wt

# 5. Extender hasta 2054 y proyectar pesos de captura
fleets <- window(fleets,start = first.yr,end   = last.yr)
myrs <- ac((maxyear - 2):maxyear)

ANEcwa.mean <-yearMeans(landings.wt(fleets$SEINE@metiers$ALL@catches$ANE)[, myrs, ])
landings.wt(fleets$SEINE@metiers$ALL@catches$ANE)[, ac(proj.yrs), ] <- ANEcwa.mean
ANEdwa.mean <-yearMeans(discards.wt(fleets$SEINE@metiers$ALL@catches$ANE)[, myrs, ])
discards.wt(fleets$SEINE@metiers$ALL@catches$ANE)[, ac(proj.yrs), ] <- ANEdwa.mean
# Variable cost
vcost.mean <- yearMeans(fleets$SEINE@metiers$ALL@vcost[, myrs, ])
fleets$SEINE@metiers$ALL@vcost[, ac(proj.yrs), ] <- vcost.mean
# Fixed cost
fcost.mean <- yearMeans(fleets$SEINE@fcost[, myrs, ])
fleets$SEINE@fcost[, ac(proj.yrs), ] <- fcost.mean
# Crew share
crewshare.mean <- yearMeans(fleets$SEINE@crewshare[, myrs, ])
fleets$SEINE@crewshare[, ac(proj.yrs), ] <- crewshare.mean
# ============================================================
# FLEETS CONTROLS
# ============================================================

n.fls.stks         <- 1
fls.stksnames      <- "ANE"
effort.models      <- "SMFB"
effort.restr.SEINE <- "ANE"
restriction.SEINE  <- "catch"
catch.models       <- "Baranov"
capital.models     <- "fixedCapital"


flq.ANE <- FLQuant(dimnames = list(
    age    = "all",
    year   = first.yr:last.yr,
    unit   = ANE.unit,
    season = 1:ns,
    iter   = 1:ni))


fleets.ctrl <- create.fleets.ctrl(
  fls                = fls,
  n.fls.stks         = n.fls.stks,
  fls.stksnames      = fls.stksnames,
  effort.models      = effort.models,
  catch.models       = catch.models,
  capital.models     = capital.models,
  flq                = flq.ANE,
  effort.restr.SEINE = effort.restr.SEINE,
  restriction.SEINE  = restriction.SEINE)

fleets.ctrl$SEINE$ANE$discard.TAC.OS <- FALSE


# 7. Calcular q
# Historical alpha
fleets[[fle]]@metiers[[met]]@catches[[stks]]@alpha[, ac(hist.yrs), ] <- 1
mean.yrs <-hist.yrs[(length(hist.yrs) - 5 + 1):length(hist.yrs)]

fleets <- calculate.q.sel.flrObjs(
  biols       = biols,
  fleets      = fleets,
  fleets.ctrl = fleets.ctrl,
  mean.yrs    = mean.yrs,
  sim.yrs     = proj.yrs)

cq <-catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])

cq[is.na(cq)]       <- 0
cq[is.infinite(cq)] <- 0

catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]]) <- cq

# 8. Completar futuro 

# Alpha
alpha.mean <-yearMeans(fleets$SEINE@metiers$ALL@catches$ANE@alpha[, myrs, ])
fleets$SEINE@metiers$ALL@catches$ANE@alpha[, ac(proj.yrs), ] <- alpha.mean

# Beta
beta.mean <-yearMeans(fleets$SEINE@metiers$ALL@catches$ANE@beta[, myrs, ])
fleets$SEINE@metiers$ALL@catches$ANE@beta[, ac(proj.yrs), ] <- beta.mean

# Effort share
effshare.mean <-yearMeans(fleets$SEINE@metiers$ALL@effshare[, myrs, ])
fleets$SEINE@metiers$ALL@effshare[, ac(proj.yrs), ] <- effshare.mean

# Capacity
fleetcapacity.mean <-yearMeans(fleets$SEINE@capacity[, myrs, ])
fleets$SEINE@capacity[, ac(proj.yrs), ] <-fleetcapacity.mean

# 9. Seasonal share
# ============================================================
# SEASONAL CATCH SHARE — REFERENCE
# ============================================================

c1s.yr <-quantSums(catchWStock(fleets,stock = "ANE"))[,,,1,]
c2s.yr <-quantSums(catchWStock(fleets,stock = "ANE"))[,,,2,]
c3s.yr <-quantSums(catchWStock(fleets,stock = "ANE"))[,,,3,]
c4s.yr <-quantSums(catchWStock(fleets,stock = "ANE"))[,,,4,]

total.catch <-seasonSums(quantSums(catchWStock(fleets,stock = "ANE")))

fleets.ctrl$seasonal.share[[1]][,,,1,] <-c1s.yr / total.catch
fleets.ctrl$seasonal.share[[1]][,,,2,] <-c2s.yr / total.catch
fleets.ctrl$seasonal.share[[1]][,,,3,] <-c3s.yr / total.catch
fleets.ctrl$seasonal.share[[1]][,,,4,] <-c4s.yr / total.catch

# para el futuro
myrs.seasonal <- ac((ass.yr - 9):ass.yr)

for(s in 1:4) {
  fleets.ctrl$seasonal.share[[1]][,ac(proj.yrs),,s,] <-
    yearMeans(fleets.ctrl$seasonal.share[[1]][,myrs.seasonal,,s,])
}

# 10. Inicializar los slot futuros
landings.n(fleets$SEINE@metiers$ALL@catches$ANE)[, ac(proj.yrs), ] <- 0
discards.n(fleets$SEINE@metiers$ALL@catches$ANE)[, ac(proj.yrs), ] <- 0
fleets$SEINE@metiers$ALL@catches$ANE@beta[, ac(proj.yrs), ] <- 1
fleets$SEINE@effort[, ac(proj.yrs), ] <- 1

cobj <- fleets$SEINE@metiers$ALL@catches$ANE
wt_final_hist_diff <- max(abs(as.numeric(landings.wt(cobj)[,ac(hist.yrs)]-wt(biols$ANE)[,ac(hist.yrs)])),na.rm=TRUE)
wt_final_proj_diff <- max(abs(as.numeric(landings.wt(cobj)[,ac(proj.yrs)]-wt(biols$ANE)[,ac(proj.yrs)])),na.rm=TRUE)
stopifnot(wt_final_hist_diff < 1e-10,wt_final_proj_diff < 1e-10)
message("Fleet/biol weight consistency: OK")

# ============================================================
# WORKING DOCUMENT OUTPUTS
# ============================================================
# Reference fishery-conditioning configuration.
# Documentation output only; no OM objects are modified.

fishery_tab_dir <- file.path(here(), "outputs", "fishery_conditioning", "tables")

dir.create(fishery_tab_dir, recursive = TRUE, showWarnings = FALSE)

fishery_conditioning_configuration <- data.frame(
  component = c(
    "Historical period",
    "Projection period",
    "Fleet",
    "Metier",
    "Stock",
    "Number of iterations",
    "Number of seasons",
    "Effort model",
    "Catch model",
    "Capital model",
    "Effort restriction",
    "Catchability reference period",
    "Catch-weight reference period",
    "Fleet-parameter reference period",
    "Seasonal-share reference period",
    "Landings selectivity multiplier",
    "Discards selectivity multiplier",
    "Future landings numbers",
    "Future discards numbers"
  ),
  value = c(
    "1989-2024",
    "2025-2054",
    "SEINE",
    "ALL",
    "ANE",
    "1000",
    "4",
    "SMFB",
    "Baranov",
    "fixedCapital",
    "catch",
    "2020-2024",
    "2022-2024",
    "2022-2024",
    "2015-2024",
    "1",
    "0",
    "Initialised at 0",
    "Initialised at 0"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  fishery_conditioning_configuration,
  file.path(fishery_tab_dir, "fishery_conditioning_configuration.csv"),
  row.names = FALSE
)

message("Fishery-conditioning WD configuration saved in: ", fishery_tab_dir)

# ============================================================
# SAVE
# ============================================================

saveRDS(ane.stock,file.path(rds_boot_dir,"ane_stock_conditioned_1000iter.rds"))
saveRDS(fleets,file.path(rds_boot_dir,"fleets_conditioned_1000iter.rds"))
saveRDS(fleets.ctrl,file.path(rds_boot_dir,"fleets_ctrl_reference_1000iter.rds"))


