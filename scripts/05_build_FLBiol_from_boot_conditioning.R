# ============================================================
# 05_build_FLBiol_from_bootstrap_conditioning.R
# Build biological OM from bootstrap-conditioned historical stock
#
# Purpose:
#   Build the FLBiol object used as the biological component of
#   the operating model from the 100 bootstrap-conditioned
#   historical SS3 trajectories.
#
#   The script:
#     - Loads the 100-iteration conditioned FLStock.
#     - Configures spawning timing and maturity.
#     - Converts the conditioned FLStock into an FLBiol.
#     - Extends biological quantities through 2025-2054.
#     - Projects weight-at-age using the 2022-2024 mean.
#     - Scales future natural mortality among OMs using
#       Lorenzen weight scaling.
#     - Performs consistency checks on historical and
#       projected biological quantities.
#
# INPUTS:
#   data/Rdata/
#     - stk_conditioned_100iter.rds
#     - iteration_map.csv
#
#   boot/data/
#     - stk_ane9aS.rds
#
# OUTPUTS:
#   data/Rdata/
#     - ane_stock_conditioned_100iter.rds
#     - biols_conditioned_100iter.rds
#     - M_future_summary.rds
#     - weight_projection_bootstrap.rds
#     - M_projection_bootstrap.rds
#
# Simulation period:
#   Historical: 1989-2024
#   Projection: 2025-2054
#
# Notes:
#   - 100 bootstrap-conditioned historical OMs are retained.
#   - Projected weight-at-age is fixed within each OM using
#     the mean of the last three historical years (2022-2024).
#   - Future M differs among OMs according to projected weight
#     and Lorenzen weight scaling, but remains constant through
#     the projection period within each OM.
# ============================================================

rm(list = ls())

# ============================================================
# LIBRARIES
# ============================================================

library(FLBEIA)
library(FLCore)
library(dplyr)
library(tidyr)
library(purrr)
library(here)

# ============================================================
# WORKING DIRECTORY
# ============================================================

wd <- here()
setwd(wd)

rds_boot_dir <- file.path(getwd(),"data","Rdata")

# ============================================================
# SIMULATION PARAMETERS
# ============================================================

first.yr    <- 1989
last.obs.yr <- 2024
proj.yr     <- 2025
proj.nyr    <- 30

hist.yrs <- first.yr:last.obs.yr
last.yr  <- proj.yr + proj.nyr - 1
proj.yrs <- proj.yr:last.yr
ass.yr   <- last.obs.yr

ns <- 4

# Actualmente estamos trabajando con
# 100 condiciones históricas distintas
n_boot <- 100
ni     <- n_boot

# ============================================================
# LOAD CONDITIONED STOCK
# ============================================================
stk <- readRDS(file.path(rds_boot_dir,"stk_conditioned_100iter.rds"))
iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map.csv"))
# ============================================================
# STOCK CONFIGURATION
# ============================================================

name(stk) <- "ANE"

desc(stk) <- "100 bootstrap-conditioned SS3 historical trajectories"

stk@harvest@units <- "f"

range(stk, "minfbar") <- 3
range(stk, "maxfbar") <- 3

# ============================================================
# SPAWNING TIMING
# ============================================================

# Toda la F de season 1 ocurre antes del desove
fspwn <- harvest.spwn(stk)
fspwn[] <- 0
fspwn[,,, "1", ] <- 1
harvest.spwn(stk) <- fspwn

# Toda la M de season 1 ocurre antes del desove
mspwn <- m.spwn(stk)
mspwn[] <- 0
mspwn[,,, "1", ] <- 1
m.spwn(stk) <- mspwn

# ============================================================
# MATURITY
# ============================================================

mat(stk)[ac(0), ] <- 0

# SSB sólo en season 2
mat(stk)[,,,1] <- 0
mat(stk)[,,,3] <- 0
mat(stk)[,,,4] <- 0

ane.stock <- stk

rm(stk)

# ============================================================
# BIOLOGICAL DATA
# ============================================================

stks <- "ANE"

ANE.age.min <- 0
ANE.age.max <- 3
ANE.unit    <- 1

ANE_range.plusgroup <- ANE.age.max
ANE_range.minyear   <- first.yr

ane <- FLBiol(
  n    = stock.n(ane.stock),
  wt   = stock.wt(ane.stock),
  m    = m(ane.stock),
  spwn = m.spwn(ane.stock),
  mat  = mat(ane.stock),
  fec  = predictModel(FLQuants(fec = ane.stock@mat * 0 + 1),model = ~fec),
  name  = stks,
  desc  = ane.stock@desc,
  range = ane.stock@range)

units(ane)$m <- units(fec(ane)) <- units(mat(ane)) <- ""

# Se amplia hasta el final de la proyección 
ane <- window(ane,start = first.yr,end   = last.yr)
ane@desc <- "100 bootstrap-conditioned SS3 historical trajectories"
biols <- FLBiols(ANE = ane)


# ============================================================
# PROJECTION PERIOD
# ============================================================

naver <- 3 - 1
maxyear <- ane.stock@range["maxyear"]
myrs <- ac((maxyear - naver):maxyear)
myrs
# Natural mortality
m(biols$ANE)[, ac(proj.yrs), ] <-yearMeans(m(biols$ANE[, myrs]))
# Weight-at-age
wt(biols$ANE)[, ac(proj.yrs), ] <-yearMeans(wt(biols$ANE[, myrs]))

# ============================================================
# FUTURE NATURAL MORTALITY UNCERTAINTY
# Lorenzen weight-scaling
# ============================================================

lorenzen_b <- 0.288

# Annual benchmark M by age
M_base_annual <- c("0" = 2.97,"1" = 1.60,"2" = 2.48,"3" = 2.48)
# FLBEIA operates with 4 seasonal time steps.
# Therefore annual M must be expressed on the seasonal scale.
M_base_season <- M_base_annual / ns
M_base_season
# ------------------------------------------------------------
# 1. Base-model reference weight
# ------------------------------------------------------------
stk_base <- readRDS(file.path(getwd(), "boot", "data", "stk_ane9aS.rds"))
wt_base_proj <- yearMeans(stock.wt(stk_base[, myrs]))
# ------------------------------------------------------------
# 2. Projected weight of bootstrap OMs
# ------------------------------------------------------------
wt_proj_boot <- wt(biols$ANE)[, ac(proj.yr), , , ,]

# Check historical weights transferred exactly from FLStock to FLBiol
wt_hist_diff <- max(abs(as.numeric(wt(biols$ANE)[,ac(hist.yrs)]-stock.wt(ane.stock)[,ac(hist.yrs)])),na.rm=TRUE)
print(wt_hist_diff)
stopifnot(wt_hist_diff < 1e-10)
# Check projected weights equal the 2022-2024 mean
wt_proj_expected <- yearMeans(wt(biols$ANE[,myrs]))
wt_proj_diff <- max(abs(as.numeric(wt(biols$ANE)[,ac(proj.yr)]-wt_proj_expected)),na.rm=TRUE)
print(wt_proj_diff)
stopifnot(wt_proj_diff < 1e-10)

# ------------------------------------------------------------
# 3. Build future seasonal M by bootstrap OM
# ------------------------------------------------------------

M_proj_boot <- m(biols$ANE)[, ac(proj.yr), , , ,]

for(i in seq_len(ni)) {
  for(s in seq_len(ns)) {
    for(a in names(M_base_season)) {
    
      W_boot <- as.numeric(wt_proj_boot[a, ac(proj.yr), , ac(s), , i])
      W_base <- as.numeric(wt_base_proj[a, , , ac(s), , 1])
      
      if(is.finite(W_boot) &&is.finite(W_base) &&W_boot > 0 &&W_base > 0) {
        
        M_proj_boot[a, ac(proj.yr), , ac(s), , i] <-M_base_season[a] *(W_base / W_boot)^lorenzen_b

      } else {
        
        # Age 0 in Q1-Q2 has weight = 0.
        # Keep benchmark seasonal M.
        M_proj_boot[a, ac(proj.yr), , ac(s), , i] <- M_base_season[a]
      }
    }
  }
}

# ------------------------------------------------------------
# 4. Aplicar la misma M a todo el período de proyección
# ------------------------------------------------------------
for(y in ac(proj.yrs)) {m(biols$ANE)[,y,,,,] <-M_proj_boot[,ac(proj.yr),,,,]}


# Maturity
mat(biols$ANE)[, ac(proj.yrs), ] <-yearMeans(mat(biols$ANE[, myrs]))
# Fecundity
fec(biols$ANE)[, ac(proj.yrs), ] <-yearMeans(fec(biols$ANE[, myrs]))
# Spawning
spwn(biols$ANE)[, ac(proj.yrs), ] <-yearMeans(spwn(biols$ANE[, myrs]))

# Comproar que los pesos varian entre iteraciones
# Los pesos deben variar entre iteraciones (tiene que ser FALSE)
all.equal(as.numeric(wt(biols$ANE)[, ac(proj.yr), , , , 1]),as.numeric(wt(biols$ANE)[, ac(proj.yr), , , , 2]))
# M futura debe variar entre iteraciones después del escalamiento por peso
# La comparación NO debe devolver TRUE
all.equal(as.numeric(m(biols$ANE)[, ac(proj.yr), , , , 1]), as.numeric( m(biols$ANE)[, ac(proj.yr), , , , 2]))
all.equal(as.numeric(m(ane.stock)[, ac(2024), , , , 1]),as.numeric(m(biols$ANE)[, ac(2024), , , , 1]))

# M futura constante temporalmente dentro de cada OM
stopifnot(isTRUE( all.equal( as.numeric(m(biols$ANE)[, ac(2025), , , , 1]),
                             as.numeric(m(biols$ANE)[, ac(2054), , , , 1]))))

M_future_summary <- as.data.frame(m(biols$ANE)[,ac(proj.yr),,,,]) %>%
  group_by(age, season) %>%
  summarise(
    median = median(data, na.rm = TRUE),
    p05 = quantile(data, 0.05, na.rm = TRUE),
    p95 = quantile(data, 0.95, na.rm = TRUE),
    .groups = "drop")

M_future_summary

# ============================================================
# SAVE REFERENCE CONDITIONING OBJECTS
# ============================================================

saveRDS(ane.stock,file.path(rds_boot_dir,"ane_stock_conditioned_100iter.rds"))
saveRDS(biols,file.path(rds_boot_dir,"biols_conditioned_100iter.rds"))
saveRDS(M_future_summary,file.path(rds_boot_dir,"M_future_summary.rds"))

# También guardamos los pesos y M proyectados por bootstrap
saveRDS(wt_proj_boot,file.path(rds_boot_dir,"weight_projection_bootstrap.rds"))
saveRDS(M_proj_boot,file.path(rds_boot_dir,"M_projection_bootstrap.rds"))

message("Reference biological conditioning saved in: ", rds_boot_dir)
