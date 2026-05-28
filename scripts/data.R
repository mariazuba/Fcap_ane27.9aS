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
seine_vcost_csv         <- "boot/data/seine_vcost.csv"
seine.met1_effshare_csv <- "boot/data/seine.met1_effshare.csv"
seine_effort_csv        <- "boot/data/seine_effort.csv"
seine_capacity_csv      <- "boot/data/seine_capacity.csv"
seine_crewshare_csv     <- "boot/data/seine_crewshare.csv"
seine_fcost_csv         <- "boot/data/seine_fcost.csv"



#'*===========================================================================*
# SIMULATION PARAMETERS                                                    ----
#'*===========================================================================*

first.yr <- 1989  # First year of the historic data.
last.obs.yr <- 2024 # último año del assessment / último año observado
proj.yr  <- 2025  # first year of projection
proj.nyr <- 5   # Number of years in the projection period

#'---- periods ----

hist.yrs <- first.yr:(proj.yr-1)   # historical period
last.yr  <- proj.yr + (proj.nyr-1) # Last year of projection
proj.yrs <- proj.yr:last.yr        # projection period
ass.yr   <- proj.yr-1              # assessment year
  

# seasons and iterations
ni <- 5 #! coded for 1 iteration
ns <- 4

#'*===========================================================================*
# LOAD HISTORICAL DATA (assessment + indices)                              ----
#'*===========================================================================*
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
ss.rec <- 3
ss.ssb <- 2

rec_h <- rec(stk0)[,,,ss.rec, drop = FALSE]
ssb_h <- ssb(stk0)[,,,ss.ssb, drop = FALSE]

mod2_bh <- FLSR(
  rec   = rec_h,
  ssb   = ssb_h,
  model = bevholt
)

fit2_bh <- fmle(
  window(mod2_bh,
         start = first.yr,
         end   = last.obs.yr)
)

SRs_bh <- list(
  ANE = FLSRsim(
    name  = "ANE",
    model = "bevholt",
    rec   = rec(ane),
    ssb   = ssb(ane)
  )
)

SRs_bh$ANE@params[] <- params(fit2_bh)

SRs_bh$ANE@proportion[,,,1,,] <- 0
SRs_bh$ANE@proportion[,,,2,,] <- 0
SRs_bh$ANE@proportion[,,,3,,] <- 1
SRs_bh$ANE@proportion[,,,4,,] <- 0

SRs_bh$ANE@timelag["year", ]   <- 0
SRs_bh$ANE@timelag["season", ] <- 2

SRs_bh$ANE@uncertainty[, ac(first.yr:last.obs.yr), , ss.rec, ] <-
  exp(residuals(fit2_bh))

residsd_bh <- sqrt(
  var(
    log(SRs_bh$ANE@uncertainty[, ac(first.yr:last.obs.yr), , ss.rec, ]),
    na.rm = TRUE
  )
)

seed <- 123
set.seed(seed)

SRs_bh$ANE@uncertainty[, ac(proj.yrs), , ss.rec, , ] <-
  exp(rnorm(length(proj.yrs) * ni, 0, residsd_bh))

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
bh_params <- params(fit2_bh)

bh_params_mat <- rbind(
  a = as.numeric(bh_params["a", ]),
  b = as.numeric(bh_params["b", ]))

colnames(bh_params_mat) <- seq_len(ncol(bh_params_mat))

bh_params_mat

seasonal_prop <- round(as.numeric(fleets.ctrl$seasonal.share[[1]][, ac("2025"), , , , 1]), 5)
seasonal_prop 

advice_Fcap.ctrl <- list()

advice_Fcap.ctrl[["ANE"]] <- list(
  HCR.model     = "annualTAC", #"FcapBpaHCR_ane",
  # Configuración del short-term forecast
  nyears        = 1,
  wts.nyears    = 3,
  fbar.nyears   = 3,
  f.rescale     = TRUE,
  # Estructura temporal biológica
  spawn.season  = 2,
  rec.season    = 3,
  # Límite técnico de búsqueda de F, equivalente al Fmax de SS3
  Fcap          = 2,
  # Proporción estacional de F
  propf         = seasonal_prop,
  # Modelo stock-reclutamiento
  sr.model      = "bevholt",
  sr.params     = bh_params_mat,
  # Puntos de referencia
  ref.pts       = matrix(
    rep(ANE_ref.pts, ni),
    nrow = length(ANE_ref.pts),
    ncol = ni,
    dimnames = list(names(ANE_ref.pts), 1:ni)))


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
  file = "data/mse/base.RData",
  compress = "xz"
)

#------------------------------------------------------------------------------
# 2. STOCK-RECRUITMENT SCENARIO
#------------------------------------------------------------------------------

SRs <- SRs_bh

save(
  SRs,
  fit2_bh,
  bh_params_mat,
  file = "data/mse/SR_bh.RData",
  compress = "xz"
)

#------------------------------------------------------------------------------
# 3. OBSERVATION SCENARIO
#------------------------------------------------------------------------------

obs.ctrl <- obs_Fcap.ctrl

save(
  obs.ctrl,
  file = "data/mse/obs_perfect.RData",
  compress = "xz"
)

#------------------------------------------------------------------------------
# 4. ASSESSMENT SCENARIO
#------------------------------------------------------------------------------

save(
  assess.ctrl,
  file = "data/mse/assess_none.RData",
  compress = "xz"
)

#------------------------------------------------------------------------------
# 5. ADVICE SCENARIO
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

Fcap_grid <- c(
  1.00,
  1.25,
  1.50,
  2.00,
  3.00
)

Besc_grid <- c(6561)

scenarios <- data.frame(
  scenario_id = seq_along(Fcap_grid),
  SR          = "bh",
  obs         = "perfect",
  assess      = "none",
  advice      = "fcap",
  Fcap        = Fcap_grid,
  Besc        = Besc_grid,
  ni          = ni,
  proj_nyr    = proj.nyr
)

write.csv(scenarios,"config/scenarios.csv", row.names = FALSE)


