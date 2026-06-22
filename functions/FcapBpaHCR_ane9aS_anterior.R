FcapBpaHCR_ane <- function(stocks, advice, advice.ctrl, year, stknm, ...) {
  
  run_flasher_fwd <- function(stk, ctrl, sr, maxF = 1000) {
    callr::r(
      func = function(stk, ctrl, sr, maxF) {
        library(FLCore)
        library(FLFishery)
        library(FLasher)
        FLasher::fwd(stk, control = ctrl, sr = sr, maxF = maxF)
      },
      args = list(stk = stk, ctrl = ctrl, sr = sr, maxF = maxF)
    )
  }
  
  get_ctrl <- function(x, default) { if (is.null(x)) default else x}
  
  ctrl_stk <- advice.ctrl[[stknm]]
  
  nyears      <- get_ctrl(ctrl_stk$nyears, 1)
  wts.nyears  <- get_ctrl(ctrl_stk$wts.nyears, 3)
  fbar.nyears <- get_ctrl(ctrl_stk$fbar.nyears, 3)
  f.rescale   <- get_ctrl(ctrl_stk$f.rescale, TRUE)
  
  spawn.season <- get_ctrl(ctrl_stk$spawn.season, 2)
  rec.season   <- get_ctrl(ctrl_stk$rec.season, 3)
  Fcap         <- get_ctrl(ctrl_stk$Fcap, 8)
  propf        <- get_ctrl(ctrl_stk$propf, rep(1/4, 4))
  ref.pts      <- ctrl_stk$ref.pts
  sr.params    <- ctrl_stk$sr.params
  
  if (inherits(stocks, "FLStock")) {
    stk <- stocks
  } else {
    stk <- stocks[[stknm]]
  }
  
  ns <- dims(stk)$season
  ni <- dims(stk)$iter
  
  if (length(propf) != ns) {
    stop("FcapBpaHCR_ane: length(propf) debe ser igual al número de seasons.")
  }
  
  propf <- propf / sum(propf)
  
  yrs.stk <- dimnames(stk)$year
  
  if (year < 1 || year > length(yrs.stk)) {
    return(advice)
  }
  
  ass.yr.name <- yrs.stk[year]
  ass.yr      <- as.numeric(ass.yr.name)
  
  tac.yr      <- ass.yr + 1
  tac.yr.name <- as.character(tac.yr)
  
  if (!tac.yr.name %in% dimnames(advice$TAC)$year) {
    return(advice)
  }
  
  stk <- window(stk, end = ass.yr)
  
  stk <- stf(
    stk,
    nyears      = nyears,
    wts.nyears  = wts.nyears,
    fbar.nyears = fbar.nyears,
    f.rescale   = f.rescale
  )
  
  if (!tac.yr.name %in% dimnames(stk)$year) {
    stop("FcapBpaHCR_ane: stf() no generó el año TAC ", tac.yr.name)
  }
  
  # ----------------------------------------------------------
  # Asegurar patrón de F no nulo en el año TAC
  # ----------------------------------------------------------
  
  harv_ref_yr <- ass.yr.name
  
  if (all(harvest(stk)[, harv_ref_yr, , , , ] == 0, na.rm = TRUE)) {
    
    yrs_h <- dimnames(stk)$year
    
    cand <- yrs_h[as.numeric(yrs_h) < as.numeric(harv_ref_yr)]
    
    cand <- cand[sapply(cand, function(y) {
      any(harvest(stk)[, y, , , , ] > 0, na.rm = TRUE)
    })]
    
    if (length(cand) == 0) {
      stop("FcapBpaHCR_ane: no hay año histórico con harvest > 0.")
    }
    
    harv_ref_yr <- tail(cand, 1)
  }
  
  harvest(stk)[, tac.yr.name, , , , ] <- harvest(stk)[, harv_ref_yr, , , , ]
  
  # ----------------------------------------------------------
  # Escapement target: SSB = Bpa
  # ----------------------------------------------------------
  
  Bescapement <- as.numeric(ref.pts["Bpa", 1])
  
  a <- as.numeric(sr.params["a", ])
  b <- as.numeric(sr.params["b", ])
  
  sr <- rec(stk)[, , , rec.season, drop = FALSE]
  sr[, tac.yr.name] <- (a * Bescapement) / (b + Bescapement)
  
  # Fbase es la F anual de decisión.
  # La F aplicada en season 1 será Fbase * propf[1].
  # Por eso el límite superior para Fbase es Fcap.
  
  F_upper <- Fcap
  
  ssb_at_Fbase <- function(Fbase_value) {
    
    df_tmp <- data.frame(
                  year   = rep(tac.yr, ns),
                  quant  = rep("f", ns),
                  value  = Fbase_value * propf,
                  season = seq_len(ns))
    
    ctrl_tmp <- FLasher::fwdControl(df_tmp, iters = ni)
    
    stk_tmp <- run_flasher_fwd(
                    stk  = stk,
                    ctrl = ctrl_tmp,
                    sr   = sr,
                    maxF = 1000)
    
    as.numeric(ssb(stk_tmp)[, tac.yr.name, , spawn.season])
  }
  
  Fbase <- numeric(ni)
  
  for (i in seq_len(ni)) {
    
f_obj <- function(Fbase_value) {
      ssb_at_Fbase(Fbase_value)[i] - Bescapement
    }
    
    ssb_low <- f_obj(0)
    ssb_up  <- f_obj(F_upper)
    
    if (ssb_low <= 0) {
      Fbase[i] <- 0
    } else if (ssb_up >= 0) {
      Fbase[i] <- F_upper
    } else {
      Fbase[i] <- uniroot(
        f_obj,
        lower = 0,
        upper = F_upper
      )$root
    }
  }
  
  # ----------------------------------------------------------
  # Forward final con Fbase distribuida por seasons
  # ----------------------------------------------------------
  
  df_final <- data.frame(
               year   = rep(tac.yr, ns),
               quant  = rep("f", ns),
               value  = NA,
               season = seq_len(ns))
  
  ctrl_final <- FLasher::fwdControl(df_final, iters = ni)
  
  for (s in seq_len(ns)) {
    ctrl_final@iters[s, "value", ] <- Fbase * propf[s]
  }
  
  stk_fwd <- run_flasher_fwd(
            stk  = stk,
            ctrl = ctrl_final,
            sr   = sr,
            maxF = 1000)
  
  res_catch <- seasonSums(catch(stk_fwd))[, tac.yr.name]
  
  advice$TAC[stknm, tac.yr.name] <- res_catch
  
  # ----------------------------------------------------------
  # Guardar F de decisión comparable con Fcap
  # ----------------------------------------------------------
  
  if (is.null(advice$Fadv)) {
    advice$Fadv <- advice$TAC
    advice$Fadv[] <- NA
    units(advice$Fadv) <- "f"
  }
  
  advice$Fadv[stknm, tac.yr.name] <- Fbase
  
  return(advice)
}

environment(FcapBpaHCR_ane) <- asNamespace("FLBEIA")
assignInNamespace(x     = "annualTAC",value = FcapBpaHCR_ane, ns    = "FLBEIA")
