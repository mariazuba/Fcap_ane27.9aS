################################################################################
#
# OM anchovy 9aS
#
################################################################################



rm(list=ls())

#'*===========================================================================*
# LOAD LIBRARIES AND FUNCTIONS                                             ----
#'*===========================================================================*
library(FLBEIA)
library(ggplotFL)
library(r4ss)
library(icesTAF)
library(reshape2)
library(dplyr)
library(callr)
library(here)


#'*===========================================================================*
# WORKING DIRECTORY                                                        ----
#'*===========================================================================*

wd <-here() # main directory
setwd(wd)

data_om<-paste0("data")
mkdir(data_om)

# Data requerida 
stk_ane9aS_rds          <- "boot/data/stk_ane9aS.rds"
ss3_ane9aS.rds          <- "boot/data/ss3_ane9aS.rds"
seine_vcost_csv         <- "boot/data/seine_vcost.csv"
seine.met1_effshare_csv <- "boot/data/seine.met1_effshare.csv"
seine_effort_csv        <- "boot/data/seine_effort.csv"
seine_capacity_csv      <- "boot/data/seine_capacity.csv"
seine_crewshare_csv     <- "boot/data/seine_crewshare.csv"
seine_fcost_csv         <- "boot/data/seine_fcost.csv"

dir.create("data/mse", recursive = TRUE, showWarnings = FALSE)
dir.create("config", recursive = TRUE, showWarnings = FALSE)

#'*===========================================================================*
# SIMULATION PARAMETERS                                                    ----
#'*===========================================================================*

first.yr    <- 1989  # First year of the historic data.
last.obs.yr <- 2024 # último año del assessment / último año observado
proj.yr     <- 2025  # first year of projection
proj.nyr    <- 30   # Number of years in the projection period

#'---- periods ----

hist.yrs <- first.yr:(proj.yr-1)   # historical period
last.yr  <- proj.yr + (proj.nyr-1) # Last year of projection
proj.yrs <- proj.yr:last.yr        # projection period
ass.yr   <- proj.yr-1              # assessment year
  

# seasons and iterations
ni <- 1000
ns <- 4

seed <- 123

experiment_name <- paste0("test_", ni, "iter_", proj.nyr, "year")

#'*===========================================================================*
# LOAD HISTORICAL DATA (assessment + indices)                              ----
#'*===========================================================================*
ss3       <- readRDS(ss3_ane9aS.rds)
stk0      <- readRDS(stk_ane9aS_rds)
stk       <- propagate(stk0, ni) # Extend to the number of iterations

name(stk)<- "ANE"
desc(stk)<-"WGHANSA2025 assessment output"
dimnames(stk)
stk@harvest@units <-'f'

# el desove se asume al inicio de la season 2:
fspwn <- harvest.spwn(stk) 
fspwn[] <- 0
fspwn[,,, "1", ] <- 1       # toda la F de la season 1 es "antes del desove"
harvest.spwn(stk) <- fspwn

# Repite para la mortalidad natural (M):
mspwn <- m.spwn(stk)
mspwn[] <- 0
mspwn[,,, "1", ] <- 1       # toda la M de la season 1 es "antes del desove"
m.spwn(stk) <- mspwn

mat(stk)[ac(0),]<-0 # edad 0=0 maduros
mat(stk)[,,,1]<-0 # maturity in the 1st season is 0, no ssb in 1st season
mat(stk)[,,,3]<-0 # maturity in the 3 season is 0, no ssb in 3 season
mat(stk)[,,,4]<-0 # maturity in the 4 season is 0, no ssb in 4 season

range(stk, "minfbar") <- 3
range(stk, "maxfbar") <- 3

#change name of FLR object
ane.stock <- stk
rm(stk)

#'*===========================================================================*
# BIOLOGICAL DATA ----                                                         ----
#'*===========================================================================*
#'---- stocks ----

stks <- "ANE"

# stock ANE
ANE.age.min  <- 0
ANE.age.max  <- 3
ANE.unit     <- 1
ANE_range.plusgroup <- ANE.age.max
ANE_range.minyear   <- first.yr

ane <- FLBiol( n=stock.n(ane.stock), 
               wt=stock.wt(ane.stock), 
               m=m(ane.stock), 
               spwn=m.spwn(ane.stock),
               mat=mat(ane.stock),
               fec=predictModel(FLQuants(fec=ane.stock@mat*0+1), model=~fec),
               name=stks, 
               desc=ane.stock@desc, 
               range=ane.stock@range)

units(ane)$m <- units(fec(ane)) <- units(mat(ane)) <- ""

ane <- window(ane, start=first.yr, end=last.yr ) 

ane@desc <- "WGHANSA2025 assessment output"


#'---- FLBiols object ---- 

biols <- FLBiols(ANE = ane)
#'*===========================================================================*
## Projection period ----
#'*===========================================================================*

# Values set as the mean of the last 63years (mean.yrs = hl6.yr)
naver=3-1
maxyear<-ane.stock@range["maxyear"]
myrs <- ac((maxyear-naver):maxyear)

# - Natural mortality 
m(biols$ANE)[,ac(proj.yrs),]<-yearMeans(m(biols$ANE[,myrs]))
# - Weights at age in the stock
wt(biols$ANE)[,ac(proj.yrs),]<-yearMeans(wt(biols$ANE[,myrs]))
# - Maturity at age
mat(biols$ANE)[,ac(proj.yrs),]<-yearMeans(mat(biols$ANE[,myrs]))
# - Fecundity
fec(biols$ANE)[,ac(proj.yrs),]<-yearMeans(fec(biols$ANE[,myrs]))
# - Spawning
spwn(biols$ANE)[,ac(proj.yrs),]<-yearMeans(spwn(biols$ANE[,myrs]))
units(biols$ANE@n*wt(biols$ANE))

#'*===========================================================================*
## Reference points ----
#'*===========================================================================*
# reference points for ANE (tons) - WGHANSA2024
Fpa   <- NA 
Flim  <- NA 
Bpa   <- 6561
Blim  <- 4721 
Blow  <- NA
Bloss <- NA
Fmsy  <- NA
Floss <- NA 
Flow  <- NA 
Bmsy  <- NA

ANE_ref.pts <- c( Fpa=Fpa, Flim=Flim, Floss= Floss, Flow=Flow, Bpa=Bpa, Blim=Blim, Fmsy=Fmsy, Bmsy=Bmsy, Bloss=Bloss)

#'*===========================================================================*
# BIOLS CONTROLS ----                                                          
#'*===========================================================================*

growth.model <- c("ASPG_Baranov")
biols.ctrl   <- create.biols.ctrl (stksnames=stks,growth.model=growth.model)

#'*===========================================================================*
## SRs data and model----
#'*===========================================================================*

ss.rec <- 3   # recruitment season
ss.ssb <- 2   # spawning season

rec_h <- rec(stk0)[,,,ss.rec, drop = FALSE]
ssb_h <- ssb(stk0)[,,,ss.ssb, drop = FALSE]

# Years used to condition the SR relationship
fit_yrs <- first.yr:last.obs.yr

# -----------------------------------------------------------------------------
# 1. Fit alternative SR models
# -----------------------------------------------------------------------------

# Base model: Beverton-Holt
mod_bh <- FLSR(rec = rec_h, ssb = ssb_h, model = bevholt)
fit_bh <- fmle(window(mod_bh, start = first.yr, end = last.obs.yr))
# Alternative 1: segmented regression
mod_seg <- FLSR(rec = rec_h, ssb = ssb_h, model = "segreg")
fit_seg <- fmle(window(mod_seg, start = first.yr, end = last.obs.yr))
# Alternative 2: segmented regression with breakpoint fixed at Blim
mod_segBlim   <- FLSR(rec = rec_h, ssb = ssb_h, model = "segreg")
fit_segBlim_a <- fmle(window(mod_segBlim, start = first.yr, end = last.obs.yr),
                 fixed = list(b = ANE_ref.pts[["Blim"]]), method = "Brent", lower = 35, upper = 10000)
fit_segBlim   <- fmle(mod_segBlim,fixed = list(a = fit_segBlim_a@params[1], b = ANE_ref.pts[["Blim"]]))
# Alternative 3: segmented regression with breakpoint fixed at Bpa
mod_segBpa   <- FLSR(rec = rec_h, ssb = ssb_h, model = "segreg")
fit_segBpa_a <- fmle(window(mod_segBpa, start = first.yr, end = last.obs.yr),
                fixed = list(b = ANE_ref.pts[["Bpa"]]), method = "Brent", lower = 35, upper = 10000)
fit_segBpa   <- fmle(mod_segBpa, fixed = list(a = fit_segBpa_a@params[1], b = ANE_ref.pts[["Bpa"]]))

save(fit_bh, fit_seg, fit_segBlim, fit_segBpa, file = "data/mse/SRs_models.RData")


#'*===========================================================================*
## FLBEIA input object: SRs ----
#'*===========================================================================*
# Helper function to create FLBEIA SR object
make_SRs <- function(fit, model_name, stk, stock_name = "ANE", ss.rec = 3, ss.ssb = 2) {
  
  SRs <- list()
  SRs[[stock_name]] <- FLSRsim(name=stock_name, model=model_name, rec=rec(stk), ssb=ssb(stk))
  SRs[[stock_name]]@params[] <- params(fit)
  
  # Recruitment enters in season 3
  SRs[[stock_name]]@proportion[] <- 0
  SRs[[stock_name]]@proportion[,,,ss.rec,,] <- 1
  
  # Spawning in season 2 produces recruitment in season 3 of the same year
  SRs[[stock_name]]@timelag["year", ]   <- 0
  SRs[[stock_name]]@timelag["season", ] <- ss.ssb
  
  SRs[[stock_name]]@uncertainty[] <- 1
  
  return(SRs)
}

SRs_bh      <- make_SRs(fit_bh,      "bevholt", ane, ss.rec = ss.rec, ss.ssb = ss.ssb)
SRs_seg     <- make_SRs(fit_seg,     "segreg",  ane, ss.rec = ss.rec, ss.ssb = ss.ssb)
SRs_segBlim <- make_SRs(fit_segBlim, "segreg",  ane, ss.rec = ss.rec, ss.ssb = ss.ssb)
SRs_segBpa  <- make_SRs(fit_segBpa,  "segreg",  ane, ss.rec = ss.rec, ss.ssb = ss.ssb)


#'*===========================================================================*
## Recruitment uncertainty scenarios ----
#'*===========================================================================*

# Historical residuals from base BH model
rec_dev_hist <- exp(residuals(fit_bh))

tmp_hist <- SRs_bh$ANE@uncertainty[, ac(fit_yrs), , ss.rec, , ]
tmp_hist[] <- as.numeric(rec_dev_hist)

SRs_bh$ANE@uncertainty[, ac(fit_yrs), , ss.rec, , ] <- tmp_hist
# Historical recruitment variability estimated from BH residuals
sigmaR_hist <- sd(as.numeric(residuals(fit_bh)), na.rm = TRUE)
# sigmaR_hist <- sqrt(var(log(SRs_bh$ANE@uncertainty[, ac(fit_yrs), , ss.rec, ]),na.rm = TRUE))

sigmaR_scenarios <- c(
  hist = sigmaR_hist,
  SS3  = 0.33,
  s05  = 0.50,
  s07  = 0.70
)


# La relación stock-reclutamiento de Beverton-Holt se mantuvo como modelo base 
# porque es coherente con el modelo de evaluación y mostró un ajuste comparable al 
# modelo segmentado libre. Aunque los modelos segmentados con punto de quiebre fijo
# en Blim y Bpa presentaron valores de AIC más bajos, se trataron como escenarios 
# de robustez estructural porque sus puntos de quiebre fueron impuestos en niveles 
# de biomasa de referencia y no estimados libremente desde los datos.

# Extreme recruitment event based on the historical recruitment distribution
rec_hist_values <- as.numeric(rec_h)

extreme_prob <- 0.05

extreme_multiplier_hist <- as.numeric(
  quantile(rec_hist_values, probs = extreme_prob, na.rm = TRUE) /
    median(rec_hist_values, na.rm = TRUE)
)


#'*===========================================================================*
## Function to simulate future recruitment deviations ----
#'*===========================================================================*

add_rec_uncertainty <- function(SRs,
                                proj.yrs,
                                sigmaR,
                                stock_name = "ANE",
                                ss.rec = 3,
                                seed = 123,
                                extreme_event = FALSE,
                                extreme_prob = 0.05,
                                extreme_multiplier = 0.25) {
  
  set.seed(seed)
  
  target   <- SRs[[stock_name]]@uncertainty[, ac(proj.yrs), , ss.rec, , ]
  rec_devs <- exp(rnorm(length(target), mean = 0, sd = sigmaR))
  
  if (extreme_event) {
    is_extreme <- rbinom(length(target), size = 1, prob = extreme_prob)
    rec_devs   <- ifelse(is_extreme == 1,rec_devs * extreme_multiplier,rec_devs)
  }
  
  target[] <- rec_devs
  SRs[[stock_name]]@uncertainty[, ac(proj.yrs), , ss.rec, , ] <- target
  
  return(SRs)
}



#'*===========================================================================*
## Build SR uncertainty scenarios ----
#'*===========================================================================*

SRs_bh_hist <- add_rec_uncertainty(
  SRs = SRs_bh,
  proj.yrs = proj.yrs,
  sigmaR = sigmaR_scenarios["hist"],
  ss.rec = ss.rec,
  seed = 123
)

SRs_bh_SS3 <- add_rec_uncertainty(
  SRs = SRs_bh,
  proj.yrs = proj.yrs,
  sigmaR = sigmaR_scenarios["SS3"],
  ss.rec = ss.rec,
  seed = 123
)

SRs_bh_s05 <- add_rec_uncertainty(
  SRs = SRs_bh,
  proj.yrs = proj.yrs,
  sigmaR = sigmaR_scenarios["s05"],
  ss.rec = ss.rec,
  seed = 123
)

SRs_bh_s07 <- add_rec_uncertainty(
  SRs = SRs_bh,
  proj.yrs = proj.yrs,
  sigmaR = sigmaR_scenarios["s07"],
  ss.rec = ss.rec,
  seed = 123
)

SRs_bh_extreme <- add_rec_uncertainty(
  SRs                = SRs_bh,
  proj.yrs           = proj.yrs,
  sigmaR             = sigmaR_scenarios["hist"],
  ss.rec             = ss.rec,
  seed               = 123,
  extreme_event      = TRUE,
  extreme_prob       = extreme_prob,
  extreme_multiplier = extreme_multiplier_hist
)

#'*===========================================================================*
## Save FLBEIA SR objects ----
#'*===========================================================================*

save(
  SRs_bh_hist,
  SRs_bh_SS3,
  SRs_bh_s05,
  SRs_bh_s07,
  SRs_bh_extreme,
  SRs_seg,
  SRs_segBlim,
  SRs_segBpa,
  sigmaR_scenarios,
  file = "data/mse/SRs_uncertainty_scenarios.RData"
)

#'*===========================================================================*
# MAIN CONTROLS ----                                                            ----
#'*===========================================================================*

main.ctrl           <- list()
#main.ctrl$sim.years <- c(initial = proj.yr, final = last.yr)
main.ctrl$sim.years <- c(initial = ass.yr, final = last.yr)
main.ctrl   
#'*===========================================================================*
# Biomass dynamics (for stocks in biomass) ----
#'*===========================================================================*
BDs <- NULL
#'*===========================================================================*
# FLEET DATA ----
#'*===========================================================================*
#Note: we didn't use the function create.fleets.array

###---- fleets ----
fls <- "SEINE"

###---- metiers ----
SEINE.mets <- "ALL"

###---- metiers * stocks ----
SEINE.ALL.stks <- "ANE"

#---- FLFleets object ----
# Populate the FLCatchExt object with catch related data taken from the FLStock object
catch <- FLCatchExt(name = names(biols),
                    landings.n   = landings.n(ane.stock[, ac(hist.yrs)]), 
                    landings     = landings(ane.stock[, ac(hist.yrs)]), 
                    landings.wt  = landings.wt(ane.stock[, ac(hist.yrs)]),
                    discards.n   = discards.n(ane.stock[, ac(hist.yrs)]), 
                    discards     = discards(ane.stock[, ac(hist.yrs)]),
                    discards.wt  = discards.wt(ane.stock[, ac(hist.yrs)]))

# Create an FLCatchesExt() object. Since there is only one fleet with one metier, it is a simple one.
catches <- FLCatchesExt(catch)
names(catches) <- SEINE.ALL.stks

# Next create a FLMetierExt object.
m <- FLMetierExt(catches = catches,name = SEINE.mets)


# vcost for metier
seine.vcost     <- read.csv(file = seine_vcost_csv) %>% filter(year<=maxyear)
seine.vcost.flq <- as.FLQuant(seine.vcost)
seine.vcost.flq <- propagate(seine.vcost.flq, ni)

# effort share for metier
seine.met1_effshare     <- read.csv(file = seine.met1_effshare_csv) %>% filter(year<=maxyear)
seine.met1_effshare.flq <- as.FLQuant(seine.met1_effshare)
seine.met1_effshare.flq <- propagate(seine.met1_effshare.flq, ni)

# Set the effort share and variable cost for the metier as defined in the flquants
m@effshare  <-seine.met1_effshare.flq 
m@vcost <- seine.vcost.flq

# m <- window(m,start=first.yr,end=last.yr)
metiers <- FLMetiersExt(m)
names(metiers) <- SEINE.mets

# Next define some elements ascribed to the entire fleet
# fleet effort
seine.effort     <- read.csv(file = seine_effort_csv) %>% filter(year<=maxyear)
seine.effort$data<-1
seine.effort.flq <- as.FLQuant(seine.effort)
seine.effort.flq <- propagate(seine.effort.flq, ni)

#fleet capacity, assume a high value
seine.capacity     <- read.csv(file = seine_capacity_csv) %>% filter(year<=maxyear)
seine.capacity$data<-5000
seine.capacity.flq <- as.FLQuant(seine.capacity)
seine.capacity.flq <- propagate(seine.capacity.flq, ni)

# fleet crewshare
seine.crewshare     <- read.csv(file = seine_crewshare_csv) %>% filter(year<=maxyear)
seine.crewshare.flq <- as.FLQuant(seine.crewshare)
seine.crewshare.flq <- propagate(seine.crewshare.flq, ni)

#fcost
seine.fcost     <- read.csv(file = seine_fcost_csv) %>% filter(year<=maxyear)
seine.fcost.flq <- as.FLQuant(seine.fcost)
seine.fcost.flq <- propagate(seine.fcost.flq, ni)

# Create the FLFleetExt and FLFleetsExt objects
fleet <- FLFleetExt(metiers   = metiers,
                    name      = fls,
                    effort    = seine.effort.flq, 
                    fcost     = seine.fcost.flq,
                    capacity  = seine.capacity.flq, 
                    crewshare = seine.crewshare.flq)

fleets <- FLFleetsExt(fleet)
names(fleets) <- "SEINE"

# Asignar nombres por claridad
fle <- "SEINE"
met <- "ALL"

# --- Asegurar que las unidades de medida sean coherentes con biols ---
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("landings.n", "discards.n")]   <- units(biols[[stks]])$n
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("landings.wt", "discards.wt")] <- units(biols[[stks]])$wt
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("landings", "discards")]       <- units(biols[[stks]]@n * wt(biols[[stks]]))
units(fleets[[fle]]@metiers[[met]]@catches[[stks]])[c("alpha", "beta", "catch.q")]   <- "1"

##################################################################-
# Definition of landings and discards selectivity
#'* Here I assume that selectivity is zero at age 0 all over the year, and discards selectivity is 0 for all ages, so no discards
#' Esta selectividad no se refiere a la maniobra de pesca, sino que se usa para desagregar las capturas en descargas y descartes
fleets[[fle]]@metiers[[met]]@catches[[stks]]@landings.sel[] <- 1
fleets[[fle]]@metiers[[met]]@catches[[stks]]@discards.sel[] <- 0
fleets[[fle]]@metiers[[met]]@catches[[stks]]@discards.wt <- fleets[[fle]]@metiers[[met]]@catches[[stks]]@landings.wt

#############################################################-
# Extend the FLFleetsExt up to last.yr 
fleets <- window(fleets, start=first.yr, end=last.yr )

#============================================================================-
# Initiate population of the projection period for the fleet data
#============================================================================-

# mean of last 3 years
myrs <- ac((maxyear-2):maxyear)

#'1. Peso medio en los desembarques (`landings.wt`)
ANEcwa.mean <- yearMeans(landings.wt(fleets$SEINE@metiers$ALL@catches$ANE)[,myrs,])
landings.wt(fleets$SEINE@metiers$ALL@catches$ANE)[,ac(proj.yrs),] <- ANEcwa.mean
#'2. Peso medio en descartes (`discards.wt`)
ANEdwa.mean <- yearMeans(discards.wt(fleets$SEINE@metiers$ALL@catches$ANE)[,myrs,])
discards.wt(fleets$SEINE@metiers$ALL@catches$ANE)[,ac(proj.yrs),] <- ANEdwa.mean

#'*===========================================================================*
# FLEETS CONTROLS                                                          ----
#'*===========================================================================*

n.fls.stks         <- 1                # Define que hay 1 stock por flota (es decir, cada flota captura solo un stock).
fls.stksnames      <- "ANE"            # Define el nombre del stock capturado por la flota.
effort.models      <- "SMFB"           #  "SMFB" si modelas comportamiento
effort.restr.SEINE <- "ANE"            # Define qué stock restringe el esfuerzo en la flota
restriction.SEINE  <- "catch"          # el esfuerzo depende del nivel de captura deseado o permitido.
catch.models       <- "Baranov"# "CobbDouglasAge" # Define el modelo de captura para el stock. "CobbDouglasAge" es un modelo en el que la captura depende de
capital.models     <- "fixedCapital"   # "fixedCapital" implica que la capacidad de la flota es constante en el tiempo (no hay inversión/desinversión).

#'8. FLQuant vacío (`flq`)
#Crea un objeto FLQuant vacío con dimensiones adecuadas para las variables que dependen del tiempo y temporada, pero sin valores aún (NA).
flq.ANE <- FLQuant(dimnames = list(age    = "all",
                                   year   = first.yr:last.yr,
                                   unit   = ANE.unit,
                                   season = 1:ns,
                                   iter   = 1:ni))
#'*===========================================================================*
# FLBEIA input object: fleets.ctrl
#'*===========================================================================*
fleets.ctrl <- create.fleets.ctrl( fls                = fls,          # nombres de las flotas
                                   n.fls.stks         = n.fls.stks,   # número de stocks por flota
                                   fls.stksnames      = fls.stksnames ,# nombre del stock capturado por la flota 
                                   effort.models      = effort.models , # modelo de esfuerzo, "fixedEffort" o  "SMFB" si modelas comportamiento
                                   catch.models       = catch.models,     # modelo de captura,
                                   capital.models     = capital.models,       # modelo de capital
                                   flq                = flq.ANE,              # un FLQuant con la dimensión adecuada (se usa para inicializar effort, capacity, etc.)
                                   effort.restr.SEINE = effort.restr.SEINE,                # indica qué stock restringe el esfuerzo de la flota 
                                   restriction.SEINE  = restriction.SEINE  )             # tipo de restricción ("catch", "landings", etc.) 

fleets.ctrl$SEINE$ANE$discard.TAC.OS <- FALSE
# fleets.ctrl$seine$ane9as$catch.model <- ifelse(biols.ctrl$ane9as$growth.model=="ASPG", "CobbDouglasAge", "Baranov") # de alfonso

###################################################################################################-
# Continuación con la configuración del objeto fleet una vez definido el catch.model en el fleet.ctrl

# Fix alpha and beta params of Cobb-Douglas function to 1. If we have selected the Baranov catch equation, we only need to set the alpha param as 1.
fleets[[fle]]@metiers[[met]]@catches[[stks]]@alpha[, ac(hist.yrs), ] <- 1
#fleets[[fle]]@metiers[[met]]@catches[[stks]]@beta[,  ac(hist.yrs), ] <- 1

# Estimate catchability using the function calculate.q.sel.flrObjs, which, 
# besides estimating q using fleet and stock data in the historic period, 
# populates q in the projection period based in the mean q at age in the period specified, as accomplished below.
mean.yrs <- hist.yrs[(length(hist.yrs)-5+1):length(hist.yrs)] # Define the number of years to estimate the mean catchability to be used in the projection
fleets   <- calculate.q.sel.flrObjs(biols=biols, fleets=fleets, fleets.ctrl=fleets.ctrl, mean.yrs=mean.yrs, sim.yrs=proj.yrs)  # apply function calculate.q.sel.flrObjs() to estimate and project catchability

# Dado que algunos valores de capturas y abudnancias son 0, se obtiene q=NA. Para evitar problemas posteriores estos NA se transforman a 0
# catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])[
#   is.infinite(catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]]))] <- 0

cq <- catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])
cq[is.na(cq)] <- 0
cq[is.infinite(cq)] <- 0
catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]]) <- cq

# plot of mean catchability by year and iteration, at age in seasons 1 and 2
q1 <- as.vector(yearMeans(iterMeans(catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])[,,,1,,])))
q2 <- as.vector(yearMeans(iterMeans(catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])[,,,2,,])))
q3 <- as.vector(yearMeans(iterMeans(catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])[,,,3,,])))
q4 <- as.vector(yearMeans(iterMeans(catch.q(fleets[[fle]]@metiers[[met]]@catches[[stks]])[,,,4,,])))
ages   <- rep(c(ANE.age.min:ANE.age.max),ns)
season <- rep(c(1:ns), each=length(unique(ages)))
qdat   <- data.frame(season=season, age=ages, q=c(q1,q2,q3,q4))

#============================================================================-
# Continue  populating the projection period for the fleet data
#============================================================================-

# Parameter alpha of Cobb-Doublas catch equation
alpha.mean <- yearMeans(fleets$SEINE@metiers$ALL@catches$ANE@alpha[,myrs,])
fleets$SEINE@metiers$ALL@catches$ANE@alpha[,ac(proj.yrs),] <- alpha.mean

# Parameter beta of Cobb-Doublas catch equation
beta.mean <- yearMeans(fleets$SEINE@metiers$ALL@catches$ANE@beta[,myrs,])
fleets$SEINE@metiers$ALL@catches$ANE@beta[,ac(proj.yrs),] <- beta.mean

# Effort share
effshare.mean <- yearMeans(fleets$SEINE@metiers$ALL@effshare[,myrs,])
fleets$SEINE@metiers$ALL@effshare[,ac(proj.yrs),] <- effshare.mean

# Fleet capacity
fleetcapacity.mean <- yearMeans(fleets$SEINE@capacity[,myrs,])
fleets$SEINE@capacity[,ac(proj.yrs),] <- fleetcapacity.mean

### In this case the fleets.ctrl0$seasonal.share is modified. 
#The percentage of the catch taken on each season is estimated from the last 10 years in the historic period
c1s.yr <- quantSums(catchWStock(fleets, stock = "ANE"))[,,,1,]
c2s.yr <- quantSums(catchWStock(fleets, stock = "ANE"))[,,,2,]
c3s.yr <- quantSums(catchWStock(fleets, stock = "ANE"))[,,,3,]
c4s.yr <- quantSums(catchWStock(fleets, stock = "ANE"))[,,,4,]

fleets.ctrl$seasonal.share[[1]][,,,1,]   <- c1s.yr/ seasonSums(quantSums(catchWStock(fleets, stock = "ANE")))
fleets.ctrl$seasonal.share[[1]][,,,2,]   <- c2s.yr/ seasonSums(quantSums(catchWStock(fleets, stock = "ANE")))
fleets.ctrl$seasonal.share[[1]][,,,3,]   <- c3s.yr/ seasonSums(quantSums(catchWStock(fleets, stock = "ANE")))
fleets.ctrl$seasonal.share[[1]][,,,4,]   <- c4s.yr/ seasonSums(quantSums(catchWStock(fleets, stock = "ANE")))

# # Consider the past myrs to estimate the average percentage of the catches taken in season 1 and 4
naver=10-1
init.yr <- ass.yr #main.ctrl$sim.years[["initial"]] - 1
myrs <- ac((init.yr-naver):init.yr)

fleets.ctrl$seasonal.share[[1]][,ac(proj.yrs),,1,] <- yearMeans(fleets.ctrl$seasonal.share[[1]][,myrs,,1,])
fleets.ctrl$seasonal.share[[1]][,ac(proj.yrs),,2,] <- yearMeans(fleets.ctrl$seasonal.share[[1]][,myrs,,2,])
fleets.ctrl$seasonal.share[[1]][,ac(proj.yrs),,3,] <- yearMeans(fleets.ctrl$seasonal.share[[1]][,myrs,,3,])
fleets.ctrl$seasonal.share[[1]][,ac(proj.yrs),,4,] <- yearMeans(fleets.ctrl$seasonal.share[[1]][,myrs,,4,])

#=============================================================================
# Catch seasonal share scenarios
#=============================================================================

make_seasonal_share <- function(fleets.ctrl, fleets, stock = "ANE", proj.yrs, scenario = "recent10", ass.yr) {
  
  total_catch <- seasonSums(quantSums(catchWStock(fleets, stock = stock)))
  
  # Historical annual seasonal proportions
  prop_hist <- fleets.ctrl$seasonal.share[[1]]
  
  if (scenario == "recent10") {
    
    yrs <- ac((ass.yr - 9):ass.yr)
    for (s in 1:4) {prop_hist[, ac(proj.yrs), , s, ] <- yearMeans(prop_hist[, yrs, , s, ])}
    
  } else if (scenario == "recent3") {
    
    yrs <- ac((ass.yr - 2):ass.yr)
    for (s in 1:4) {prop_hist[, ac(proj.yrs), , s, ] <- yearMeans(prop_hist[, yrs, , s, ])}
    
  } else if (scenario == "historical") {
    
    yrs <- ac(first.yr:ass.yr)
    for (s in 1:4) {prop_hist[, ac(proj.yrs), , s, ] <- yearMeans(prop_hist[, yrs, , s, ])}
    
  } else if (scenario == "high_Q2Q3") {
    
    prop_hist[, ac(proj.yrs), , 1, ] <- 0.10
    prop_hist[, ac(proj.yrs), , 2, ] <- 0.45
    prop_hist[, ac(proj.yrs), , 3, ] <- 0.35
    prop_hist[, ac(proj.yrs), , 4, ] <- 0.10
    
  } else {
    stop("Unknown catch_prop scenario: ", scenario)
  }
  
  fleets.ctrl$seasonal.share[[1]] <- prop_hist
  
  return(fleets.ctrl)
}

catch_prop_base <- "recent10"

fleets.ctrl <- make_seasonal_share(
  fleets.ctrl = fleets.ctrl,
  fleets      = fleets,
  stock       = "ANE",
  proj.yrs    = proj.yrs,
  scenario    = catch_prop_base,
  ass.yr      = ass.yr
) 


#' =============================================================================
# Completar slots de fleet en años de proyección para evitar NA estructurales
#' =============================================================================

landings.n(fleets$SEINE@metiers$ALL@catches$ANE)[, ac(proj.yrs), ] <- 0
discards.n(fleets$SEINE@metiers$ALL@catches$ANE)[, ac(proj.yrs), ] <- 0

fleets$SEINE@metiers$ALL@catches$ANE@beta[, ac(proj.yrs), ] <- 1

fleets$SEINE@effort[, ac(proj.yrs), ] <- 1

#'*===========================================================================*
# COVARS DATA                                                              ----
#'*===========================================================================*
# Aquí se pueden incluir variables externas (covariables) que afectan el sistema, 
# como variables ambientales, económicas o de modelo de reclutamiento.
covars <- NULL

#'*===========================================================================*
# COVARS CONTROLS                                                          ----
#'*===========================================================================*

covars.ctrl <- NULL

#'*===========================================================================*
# INDICES DATA                                                           ----
#'*===========================================================================*

indices<-NULL

#'*===========================================================================*
#'*OBSERVATION CONTROLS*                                                     ----
#'*===========================================================================*
obs_Fcap.ctrl <- create.obs.ctrl(stksnames     = stks,
                                flq.ANE       = flq.ANE,
                                stkObs.models = "perfectObs")
obs_Fcap.ctrl$ANE$obs.curryr <- TRUE

#'*===========================================================================*
#'*ASSESSMENT CONTROLS*                                                    ----
#'*===========================================================================*
assess.models <- "NoAssessment"
assess.ctrl   <- create.assess.ctrl(stksnames = stks, assess.models = assess.models)
assess.ctrl$ANE$ass.curryr <- TRUE

#'*===========================================================================*
#'*ADVICE DATA*                                                              ----
#'*===========================================================================*
advice_TACNA <- list()

advice_TACNA$TAC <- FLQuant(
  NA,
  dimnames = list(
    stock  = stks,
    year   = first.yr:last.yr,
    unit   = "unique",
    season = "all",
    area   = "unique",
    iter   = 1:ni
  )
)

advice_TACNA$TAC[, ac(2019:last.obs.yr), ] <- c(5278, 8856, 9459, 4383, 1892, 1733)
advice_TACNA$TAC[, ac(proj.yr), ] <- NA
units(advice_TACNA$TAC) <- "t"

advice_TACNA$TAE <- NULL

advice_TACNA$quota.share <- list()

advice_TACNA$quota.share[["ANE"]] <- FLQuant(
  1,
  dimnames = list(
    fleet  = fls,
    year   = first.yr:last.yr,
    unit   = "unique",
    season = "all",
    area   = "unique",
    iter   = 1:ni
  ))

#' ============================================================
# advice.ctrl para usar con FcapBpaHCR_ane ----
#' ============================================================

# Bescapement base
Besc_base <- Bpa

# Proporción estacional
seasonal_prop <- round(as.numeric(fleets.ctrl$seasonal.share[[1]][, ac(proj.yr), , , , 1]),5)

seasonal_prop <- seasonal_prop / sum(seasonal_prop)

# Reference points matrix
ref.pts.mat <- matrix(
  rep(ANE_ref.pts, ni),
  nrow = length(ANE_ref.pts),
  ncol = ni,
  dimnames = list(names(ANE_ref.pts), 1:ni))

advice_Fcap.ctrl <- list()

advice_Fcap.ctrl[["ANE"]] <- list(
  HCR.model       = "annualTAC",
  nyears          = 1,
  wts.nyears      = 3,
  fbar.nyears     = 3,
  f.rescale       = TRUE,
  spawn.season    = 2,
  rec.season      = 3,
  f.search.season = 1,
  Fcap            = 2,
  Besc            = Besc_base,
  propf           = seasonal_prop,
  ref.pts         = ref.pts.mat,
  sr              = fit_bh
)

advice_Fcap.ctrl$ANE$adv.year   <- NULL
advice_Fcap.ctrl$ANE$adv.season <- NULL


#===============================================================================
# GUARDAR INPUTS MSE MODULARES
#===============================================================================

dir.create("data/mse", recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------------------------
# 1. BASE OBJECTS
#------------------------------------------------------------------------------
save(
  ane.stock,
  biols,
  fleets,
  main.ctrl,
  biols.ctrl,
  fleets.ctrl,
  advice_TACNA,
  ANE_ref.pts,
  Blim,
  Bpa,
  first.yr,
  last.obs.yr,
  proj.yr,
  proj.nyr,
  last.yr,
  proj.yrs,
  ass.yr,
  ni,
  ns,
  seed,
  file = file.path("data/mse", paste0(experiment_name, ".RData")),
  compress = "xz"
)

#------------------------------------------------------------------------------
# 2. STOCK-RECRUITMENT SCENARIO
#------------------------------------------------------------------------------
bh_params <- params(fit_bh)

bh_params_mat <- rbind(
  a = as.numeric(bh_params["a", ]),
  b = as.numeric(bh_params["b", ])
)

colnames(bh_params_mat) <- seq_len(ncol(bh_params_mat))


SRs <- SRs_bh_hist
save(SRs,fit_bh,bh_params_mat,sigmaR_scenarios,
  file = "data/mse/SR_bh_hist.RData", compress = "xz")

SRs <- SRs_bh_SS3
save(SRs, fit_bh, bh_params_mat, sigmaR_scenarios,
     file = "data/mse/SR_bh_SS3.RData", compress = "xz")

SRs <- SRs_bh_s05
save(SRs, fit_bh, bh_params_mat, sigmaR_scenarios,
     file = "data/mse/SR_bh_s05.RData", compress = "xz")

SRs <- SRs_bh_s07
save(SRs, fit_bh, bh_params_mat, sigmaR_scenarios,
     file = "data/mse/SR_bh_s07.RData", compress = "xz")

SRs <- SRs_bh_extreme
save(SRs, fit_bh, bh_params_mat, sigmaR_scenarios,
     extreme_prob, extreme_multiplier_hist,
     file = "data/mse/SR_bh_extreme.RData", compress = "xz")

SRs <- SRs_seg
save(SRs, fit_seg,
     file = "data/mse/SR_seg.RData", compress = "xz")

SRs <- SRs_segBlim
save(SRs, fit_segBlim,
     file = "data/mse/SR_segBlim.RData", compress = "xz")

SRs <- SRs_segBpa
save(SRs, fit_segBpa,
     file = "data/mse/SR_segBpa.RData", compress = "xz")

#------------------------------------------------------------------------------
# Seasonal catch allocation scenarios
#------------------------------------------------------------------------------

catch_prop_values <- tibble::tribble(
  ~catch_prop, ~Q1,    ~Q2,    ~Q3,    ~Q4,
  "recent10",  NA,     NA,     NA,     NA,
  "cluster_1", 0.271,  0.534,  0.164,  0.031,
  "cluster_2", 0.238,  0.315,  0.288,  0.159,
  "cluster_3", 0.107,  0.423,  0.370,  0.100
)

catch_prop_grid_all <- c(
  "recent10",
  "cluster_1",
  "cluster_2",
  "cluster_3"
)

catch_prop_grid_run <- c(
  "cluster_1",
  "cluster_2",
  "cluster_3"
)

make_seasonal_share <- function(fleets.ctrl,
                                proj.yrs,
                                scenario,
                                ass.yr,
                                catch_prop_values) {
  
  prop_hist <- fleets.ctrl$seasonal.share[[1]]
  
  if (scenario == "recent10") {
    
    yrs <- ac((ass.yr - 9):ass.yr)
    
    for (s in 1:4) {
      prop_hist[, ac(proj.yrs), , s, ] <-
        yearMeans(prop_hist[, yrs, , s, ])
    }
    
  } else {
    
    if (!scenario %in% catch_prop_values$catch_prop) {
      stop("Unknown catch_prop scenario: ", scenario)
    }
    
    prop_vec <- catch_prop_values |>
      dplyr::filter(catch_prop == scenario) |>
      dplyr::select(Q1, Q2, Q3, Q4) |>
      unlist(use.names = FALSE)
    
    if (any(is.na(prop_vec))) {
      stop("Missing seasonal proportions for scenario: ", scenario)
    }
    
    prop_vec <- prop_vec / sum(prop_vec)
    
    for (s in 1:4) {
      prop_hist[, ac(proj.yrs), , s, ] <- prop_vec[s]
    }
  }
  
  fleets.ctrl$seasonal.share[[1]] <- prop_hist
  
  return(fleets.ctrl)
}

# Save all seasonal allocation objects, including recent10
dir.create("data/mse/catch_prop", recursive = TRUE, showWarnings = FALSE)

for (cp in catch_prop_grid_all) {
  
  fleets.ctrl.cp <- make_seasonal_share(
    fleets.ctrl       = fleets.ctrl,
    proj.yrs          = proj.yrs,
    scenario          = cp,
    ass.yr            = ass.yr,
    catch_prop_values = catch_prop_values
  )
  
  save(
    fleets.ctrl.cp,
    file = file.path("data/mse/catch_prop", paste0("fleets_ctrl_", cp, ".RData")),
    compress = "xz"
  )
}

# Keep recent10 as the base case in the objects already saved above
fleets.ctrl <- make_seasonal_share(
  fleets.ctrl       = fleets.ctrl,
  proj.yrs          = proj.yrs,
  scenario          = "recent10",
  ass.yr            = ass.yr,
  catch_prop_values = catch_prop_values
)

#------------------------------------------------------------------------------
# Observation scenario
#------------------------------------------------------------------------------

obs.ctrl <- obs_Fcap.ctrl

save(
  obs.ctrl,
  file = "data/mse/obs_perfect.RData",
  compress = "xz"
)

#------------------------------------------------------------------------------
# Assessment scenario
#------------------------------------------------------------------------------

save(
  assess.ctrl,
  file = "data/mse/assess_none.RData",
  compress = "xz"
)

#------------------------------------------------------------------------------
# Advice scenario
#------------------------------------------------------------------------------

advice.ctrl <- advice_Fcap.ctrl

save(
  advice.ctrl,
  file = "data/mse/adv_fcap.RData",
  compress = "xz"
)

message("MSE input objects saved in data/mse/")

#------------------------------------------------------------------------------
# MSE scenarios
#------------------------------------------------------------------------------

Fcap_grid <- c(0.75, 1.00, 1.25, 1.50, 2.00)

Besc_grid <- c(5000, Bpa, 8000, 10000)

SR_grid <- c(
  "bh_hist",
  "bh_SS3",
  "bh_s05",
  "bh_s07",
  "bh_extreme"
)

scenarios <- expand.grid(
  SR         = SR_grid,
  Fcap       = Fcap_grid,
  Besc       = Besc_grid,
  catch_prop = catch_prop_grid_all,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

scenarios <- scenarios |>
  dplyr::mutate(
    scenario_id = dplyr::row_number(),
    obs         = "perfect",
    assess      = "none",
    advice      = "fcap",
    ni          = ni,
    proj_nyr    = proj.nyr
  ) |>
  dplyr::select(
    scenario_id,
    SR,
    obs,
    assess,
    advice,
    Fcap,
    Besc,
    catch_prop,
    ni,
    proj_nyr
  )

stopifnot(nrow(scenarios) == length(SR_grid) *
            length(Fcap_grid) *
            length(Besc_grid) *
            length(catch_prop_grid_all))

stopifnot(all(scenarios$catch_prop %in% catch_prop_grid_all))

message("Number of scenarios to run: ", nrow(scenarios))
print(table(scenarios$catch_prop))

dir.create("config", recursive = TRUE, showWarnings = FALSE)

write.csv(
  scenarios,
  file = "config/scenarios.csv",
  row.names = FALSE
)

