perfectObs4seas <- function(biol, fleets, covars, obs.ctrl,
                            year = 1, season = NULL, stknm = NULL, ...) {
  
  st <- if (!is.null(stknm)) stknm else biol@name
  it <- dim(biol@n)[6]
  
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
  cw <- landings.wt(stk)
  idx <- catch.n(stk) > 0
  cw[idx] <- (
    landings.n(stk)[idx] * landings.wt(stk)[idx] +
      discards.n(stk)[idx] * discards.wt(stk)[idx]
  ) / catch.n(stk)[idx]
  
  cw[is.na(cw)] <- 0
  cw[is.infinite(cw)] <- 0
  catch.wt(stk) <- cw
  # Biomasa por season
  landings(stk) <- quantSums(landings.n(stk) * landings.wt(stk))
  discards(stk) <- quantSums(discards.n(stk) * discards.wt(stk))
  catch(stk)    <- landings(stk) + discards(stk)
  stock(stk) <- quantSums(stock.n(stk) * stock.wt(stk))
  
  #' ------------------------------------------------------------
  # Harvest observado
  #' ------------------------------------------------------------
  if (exists("ane.stock", envir = .GlobalEnv)) {
    
    harv0 <- harvest(get("ane.stock", envir = .GlobalEnv))
    yrs.harv <- dimnames(harv0)$year
    
    yrs.harv.keep <- intersect(yrs.keep, yrs.harv)
    yrs.missing   <- setdiff(yrs.keep, yrs.harv)
    
    harvest(stk)[, yrs.harv.keep, , , , ] <- harv0[, yrs.harv.keep, , , , ]
    
    # if (length(yrs.missing) > 0) {harvest(stk)[, yrs.missing, , , , ] <- 0 }
    if (length(yrs.missing) > 0) {
      
      last.harv.yr <- tail(yrs.harv.keep, 1)
      
      for (yy in yrs.missing) {
        harvest(stk)[, yy, , , , ] <- harvest(stk)[, last.harv.yr, , , , ]
      }
    }
    
  } else {
    warning("perfectObs4seas: no existe ane.stock. Se mantiene harvest desde as(biol, 'FLStock').")
  }
  
  harvest(stk)[is.na(harvest(stk))] <- 0
  harvest(stk)[is.infinite(harvest(stk))] <- 0
  units(harvest(stk)) <- "f"
  return(stk)
}

environment(perfectObs4seas) <- asNamespace("FLBEIA")
assignInNamespace(x     = "perfectObs", value = perfectObs4seas, ns    = "FLBEIA")


