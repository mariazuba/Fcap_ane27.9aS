#-------------------------------------------------------------------------------
# perfectObs4seas
#
# Perfect observation model for seasonal FLBEIA operating models
#
# Creates an FLStock object representing the stock perceived by the Management
# Procedure under perfect observation, while preserving the seasonal structure
# of the Operating Model.
#
# The standard FLBEIA perfectObs() implementation is primarily designed for
# annual observations. This modified version retains the full seasonal
# dimensions of the biological and fishery objects, allowing quarterly
# dynamics to be propagated consistently to the assessment/advice component.
#
# Biological quantities (stock.n, stock.wt, m, maturity and spawning
# information) are obtained directly from the FLBiol object. Fishery quantities
# (landings, discards, catch numbers and weights) are reconstructed from the
# FLFleet objects for each season and iteration.
#
# Fishing mortality stored in the returned FLStock must represent the realised
# mortality experienced by the Operating Model. Historical fishing mortality
# can be obtained from the conditioned stock, whereas projected fishing
# mortality should be reconstructed from realised catches rather than copied
# from the final historical year.
#
# Arguments:
#   biol      FLBiol object containing the true biological state of the OM.
#   fleets    FLFleets object containing realised fishery dynamics.
#   covars    FLBEIA covariates object (kept for compatibility with perfectObs).
#   obs.ctrl  Observation-control object (kept for FLBEIA compatibility).
#   year      Current observation/assessment year.
#   season    Current season; retained for compatibility with FLBEIA calls.
#   stknm     Stock name. If NULL, biol@name is used.
#   ...       Additional arguments passed by FLBEIA.
#
# Returns:
#   FLStock containing the observed stock history up to 'year', preserving
#   seasons and iterations and suitable for subsequent assessment/advice steps.
#
# Notes:
#   - No observation error is introduced by this function.
#   - Assessment error, when required in the shortcut MSE, is applied later
#     within the Management Procedure to the perceived stock only.
#   - The true Operating Model is never modified by this function.
#
# Shortcut MSE - Anchovy 9a South
#-------------------------------------------------------------------------------

perfectObs4seas <- function(biol, fleets, covars, obs.ctrl,
                            year=1, season=NULL, stknm=NULL,historical.stock=NULL,...) {

  
  st <- if (!is.null(stknm)) stknm else biol@name
  it <- dim(biol@n)[6]
  
  if(is.null(historical.stock) && exists("ane.stock",envir=.GlobalEnv)) 
    historical.stock <- get("ane.stock",envir=.GlobalEnv)
  
  yrs <- dimnames(biol@n)$year
  yrs.num <- suppressWarnings(as.numeric(yrs))
  
  if (!is.na(as.numeric(year)) && as.numeric(year) %in% yrs.num) {
    yrs.keep <- yrs[yrs.num <= as.numeric(year)]
  } else {
    year.pos <- as.integer(year)
    year.pos <- max(1, min(year.pos, length(yrs)))
    yrs.keep <- yrs[seq_len(year.pos)]
  }
  
  if (length(yrs.keep) < 1) {
    stop("perfectObs4seas: no hay años disponibles para observar.")
  }
  
  # Crear FLStock conservando seasons
  stk <- as(biol, "FLStock")
  stk <- stk[, yrs.keep, , , , ]
  stk <- propagate(iter(stk, 1), it, fill.iter = TRUE)
  
  #' ------------------------------------------------------------
  # Biología, conservando seasons
  #' ------------------------------------------------------------
  stock.n(stk)  <- unitSums(biol@n)[, yrs.keep, , , , ]
  stock.wt(stk) <- unitSums(biol@wt)[, yrs.keep, , , , ]
  m(stk)        <- unitSums(biol@m)[, yrs.keep, , , , ]
  
  mat_biol <- predict(biol@mat)
  mat(stk) <- unitSums(mat_biol)[, yrs.keep, , , , ]
  
  harvest.spwn(stk) <- unitSums(biol@spwn)[, yrs.keep, , , , ]
  m.spwn(stk)       <- unitSums(biol@spwn)[, yrs.keep, , , , ]
  
  #' ------------------------------------------------------------
  # Capturas desde fleets, conservando seasons
  #' ------------------------------------------------------------
  landings.n(stk) <- landStock(fleets, st)[, yrs.keep, 1, , , ]
  discards.n(stk) <- discStock(fleets, st)[, yrs.keep, 1, , , ]
  catch.n(stk)    <- landings.n(stk) + discards.n(stk)
  landings.wt(stk) <- wtalStock(fleets, st)[, yrs.keep, 1, , , ]
  discards.wt(stk) <- wtadStock(fleets, st)[, yrs.keep, 1, , , ]
  # catch.wt ponderado por capturas en número
  # evitando división por cero cuando catch.n = 0
  # cw <- landings.wt(stk)
  # idx <- catch.n(stk) > 0
  # cw[idx] <- (
  #   landings.n(stk)[idx] * landings.wt(stk)[idx] +
  #     discards.n(stk)[idx] * discards.wt(stk)[idx]
  # ) / catch.n(stk)[idx]
  # 
  # cw[is.na(cw)] <- 0
  # cw[is.infinite(cw)] <- 0
  # catch.wt(stk) <- cw
  
  # catch.wt ponderado por capturas en número
  num <- landings.n(stk) * landings.wt(stk) + discards.n(stk) * discards.wt(stk)
  cw <- landings.wt(stk)
  cw[] <- ifelse(as.array(catch.n(stk)) > 0,as.array(num) / as.array(catch.n(stk)),0)
  catch.wt(stk) <- cw
  
  # Biomasa por season
  landings(stk) <- quantSums(landings.n(stk) * landings.wt(stk))
  discards(stk) <- quantSums(discards.n(stk) * discards.wt(stk))
  catch(stk)    <- landings(stk) + discards(stk)
  stock(stk) <- quantSums(stock.n(stk) * stock.wt(stk))
  
  #' ------------------------------------------------------------
  # Realised fishing mortality
  #' ------------------------------------------------------------
  
  harvest(stk)[] <- 0
  
  if (!is.null(historical.stock)) {
    harv_hist <- harvest(historical.stock)
    yrs.hist <- intersect(yrs.keep,dimnames(harv_hist)$year)
    if (length(yrs.hist)>0) harvest(stk)[,yrs.hist,,,,] <- harv_hist[,yrs.hist,,,,]
  } else {
    yrs.hist <- character(0)
  }
  
  yrs.proj <- setdiff(yrs.keep,yrs.hist)
  
  
  # Reconstruct realised F from Baranov catch equation
  baranov_F <- function(N,M,C) {
    if (!is.finite(N) || !is.finite(M) || !is.finite(C)) return(NA_real_)
    if (N<=0 || C<=0) return(0)
    fobj <- function(F) F/(F+M)*(1-exp(-(F+M)))*N-C
    out <- try(uniroot(fobj,c(0,1e6))$root,silent=TRUE)
    if (inherits(out,"try-error")) NA_real_ else out
  }
  
  if (length(yrs.proj)>0) {
    ages <- dimnames(stock.n(stk))$age
    seasons <- dimnames(stock.n(stk))$season
    
    for (yy in yrs.proj) {
      for (a in ages) {
        for (s in seasons) {
          for (i in seq_len(it)) {
            N <- as.numeric(stock.n(stk)[a,yy,,s,,i])
            M <- as.numeric(m(stk)[a,yy,,s,,i])
            C <- as.numeric(catch.n(stk)[a,yy,,s,,i])
            harvest(stk)[a,yy,,s,,i] <- length(seasons)*baranov_F(N,M,C)
          }
        }
      }
    }
  }
  
   harvest(stk)[is.na(harvest(stk))] <- 0
  harvest(stk)[is.infinite(harvest(stk))] <- 0
  units(harvest(stk)) <- "f"
  return(stk)
}

environment(perfectObs4seas) <- asNamespace("FLBEIA")
assignInNamespace(x     = "perfectObs", value = perfectObs4seas, ns    = "FLBEIA")


