FcapBpaHCR_ane <- function(stocks, advice, advice.ctrl, year, stknm, ...) {
  
  run_flasher_fwd <- function(stk, ctrl, sr, maxF = 1000) {
    callr::r(
      func = function(stk, ctrl, sr, maxF) {
        library(methods)
        library(FLCore)
        library(FLFishery)
        library(FLasher)
        
        FLasher::fwd(stk, control = ctrl, sr = sr, maxF = maxF)
      },
      args = list(stk  = stk, ctrl = ctrl, sr   = sr, maxF = maxF))
  }
  
  get_ctrl <- function(x, default = NULL) {
    if (is.null(x)) return(default) 
    x
  }
  
  #------------------------------------------------------------
  # 1. Extract control block
  #------------------------------------------------------------
  
  ctrl_stk <- advice.ctrl[[stknm]]
  
  if (is.null(ctrl_stk)) {
    stop("FcapBpaHCR_ane: no control block found for stock ", stknm)
  }
  
  nyears      <- get_ctrl(ctrl_stk$nyears, 1)
  wts.nyears  <- get_ctrl(ctrl_stk$wts.nyears, 3)
  fbar.nyears <- get_ctrl(ctrl_stk$fbar.nyears, 3)
  f.rescale   <- get_ctrl(ctrl_stk$f.rescale, TRUE)
  
  spawn.season <- get_ctrl(ctrl_stk$spawn.season, 2)
  rec.season   <- get_ctrl(ctrl_stk$rec.season, 3)
  
  Fcap  <- get_ctrl(ctrl_stk$Fcap, NA)
  propf <- get_ctrl(ctrl_stk$propf, rep(1 / 4, 4))
  
  ref.pts <- ctrl_stk$ref.pts
  sr      <- ctrl_stk$sr
  
  n_search    <- get_ctrl(ctrl_stk$n_search, 8)
  diagnostics <- get_ctrl(ctrl_stk$diagnostics, FALSE)
  
  if (is.null(sr)) {
    stop("FcapBpaHCR_ane: sr object is missing in advice.ctrl.")
  }
  
  if (is.na(Fcap)) {
    stop("FcapBpaHCR_ane: Fcap must be provided.")
  }
  
  #------------------------------------------------------------
  # 2. Extract stock
  #------------------------------------------------------------
  
  if (inherits(stocks, "FLStock")) {
    stk <- stocks
  } else {
    stk <- stocks[[stknm]]
  }
  
  if (is.null(stk)) {
    stop("FcapBpaHCR_ane: stock ", stknm, " not found.")
  }
  
  ns <- dims(stk)$season
  ni <- dims(stk)$iter
  
  if (length(propf) != ns) {
    stop("FcapBpaHCR_ane: length(propf) must equal number of seasons.")
  }
  
  propf <- propf / sum(propf)
  
  if (spawn.season > ns) {
    stop("FcapBpaHCR_ane: spawn.season is larger than number of seasons.")
  }
  
  if (rec.season > ns) {
    stop("FcapBpaHCR_ane: rec.season is larger than number of seasons.")
  }
  
  #------------------------------------------------------------
  # 3. Reference points
  #------------------------------------------------------------
  
  if (is.null(ref.pts)) {
    stop("FcapBpaHCR_ane: ref.pts is missing.")
  }
  
  if (!"Bpa" %in% rownames(ref.pts)) {
    stop("FcapBpaHCR_ane: Bpa not found in ref.pts.")
  }
  
  if (!"Blim" %in% rownames(ref.pts)) {
    stop("FcapBpaHCR_ane: Blim not found in ref.pts.")
  }
  
  # Bpa is only the default. Alternative Besc values can be passed in ctrl_stk$Besc
  Bescapement <- as.numeric(get_ctrl(ctrl_stk$Besc, ref.pts["Bpa", 1]))
  Blim        <- as.numeric(ref.pts["Blim", 1])
  
  #------------------------------------------------------------
  # 4. Year handling
  #------------------------------------------------------------
  
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
  
  #------------------------------------------------------------
  # 5. Window and short-term extension
  #------------------------------------------------------------
  
  stk <- window(stk, end = ass.yr)
  
  stk <- FLasher::stf(
    stk,
    nyears      = nyears,
    wts.nyears  = wts.nyears,
    fbar.nyears = fbar.nyears,
    f.rescale   = f.rescale
  )
  
  if (!tac.yr.name %in% dimnames(stk)$year) {
    stop("FcapBpaHCR_ane: stf() did not generate TAC year ", tac.yr.name)
  }
  
  #------------------------------------------------------------
  # 6. Ensure non-zero harvest pattern in TAC year
  #------------------------------------------------------------
  
  harv_ref_yr <- ass.yr.name
  
  if (all(harvest(stk)[, harv_ref_yr, , , , ] == 0, na.rm = TRUE)) {
    
    yrs_h <- dimnames(stk)$year
    cand <- yrs_h[as.numeric(yrs_h) < as.numeric(harv_ref_yr)]
    
    cand <- cand[sapply(cand, function(y) {
      any(harvest(stk)[, y, , , , ] > 0, na.rm = TRUE)
    })]
    
    if (length(cand) == 0) {
      stop("FcapBpaHCR_ane: no historical year with harvest > 0.")
    }
    
    harv_ref_yr <- tail(cand, 1)
  }
  
  harvest(stk)[, tac.yr.name, , , , ] <- harvest(stk)[, harv_ref_yr, , , , ]
  

  
  #------------------------------------------------------------
  # 8. Helper: run one vector of Fbase values
  #------------------------------------------------------------
  
  run_one_Fvector <- function(Fbase_vec) {
    
    if (length(Fbase_vec) != ni) {
      stop("run_one_Fvector: length(Fbase_vec) must equal number of iterations.")
    }
    
    df_tmp <- data.frame(
      year   = rep(tac.yr, ns),
      quant  = rep("f", ns),
      value  = NA,
      season = seq_len(ns)
    )
    
    ctrl_tmp <- FLasher::fwdControl(df_tmp, iters = ni)
    
    for (s in seq_len(ns)) {
      ctrl_tmp@iters[s, "value", ] <- Fbase_vec * propf[s]
    }
    
    stk_tmp <- run_flasher_fwd(
      stk  = stk,
      ctrl = ctrl_tmp,
      sr   = sr,
      maxF = 1000
    )
    
    as.numeric(ssb(stk_tmp)[, tac.yr.name, , spawn.season])
  }
  
  #------------------------------------------------------------
  # 9. Vectorised binary search
  #    Objective: highest Fbase <= Fcap with SSB >= Bescapement
  #------------------------------------------------------------
  
  F_low  <- rep(0, ni)
  F_high <- rep(Fcap, ni)
  
  SSB_low  <- run_one_Fvector(F_low)
  SSB_high <- run_one_Fvector(F_high)
  
  Fbase <- numeric(ni)
  
  closed_iter <- SSB_low < Bescapement
  fcap_iter   <- SSB_high >= Bescapement
  search_iter <- !closed_iter & !fcap_iter
  
  Fbase[closed_iter] <- 0
  Fbase[fcap_iter]   <- Fcap
  
  lower <- F_low
  upper <- F_high
  
  if (any(search_iter)) {
    
    for (k in seq_len(n_search)) {
      
      mid <- (lower + upper) / 2
      
      SSB_mid <- run_one_Fvector(mid)
      
      ok <- SSB_mid >= Bescapement
      
      lower[ok]  <- mid[ok]
      upper[!ok] <- mid[!ok]
    }
    
    Fbase[search_iter] <- lower[search_iter]
  }
  
  Fbase <- pmin(Fbase, Fcap)
  
  #------------------------------------------------------------
  # 10. Final forward with Fbase distributed by season
  #------------------------------------------------------------
  
  df_final <- data.frame(
    year   = rep(tac.yr, ns),
    quant  = rep("f", ns),
    value  = NA,
    season = seq_len(ns)
  )
  
  ctrl_final <- FLasher::fwdControl(df_final, iters = ni)
  
  for (s in seq_len(ns)) {
    ctrl_final@iters[s, "value", ] <- Fbase * propf[s]
  }
  
  stk_fwd <- run_flasher_fwd(
    stk  = stk,
    ctrl = ctrl_final,
    sr   = sr,
    maxF = 1000
  )
  
  res_catch <- seasonSums(catch(stk_fwd))[, tac.yr.name]
  
  
  #------------------------------------------------------------
  # 11. Write advice
  #------------------------------------------------------------
  
  advice$TAC[stknm, tac.yr.name] <- res_catch
  
  if (is.null(advice$Fadv)) {
    advice$Fadv <- advice$TAC
    advice$Fadv[] <- NA
    units(advice$Fadv) <- "f"
  }
  
  advice$Fadv[stknm, tac.yr.name] <- Fbase
  
  if (diagnostics) {
    return(list(
      advice           = advice,
      stk              = stk,
      stk_fwd          = stk_fwd,
      ass.yr.name      = ass.yr.name,
      tac.yr.name      = tac.yr.name,
      Fcap             = Fcap,
      Fbase            = Fbase,
      propf            = propf,
      Blim             = Blim,
      Bescapement      = Bescapement,
      SSB_low          = SSB_low,
      SSB_high         = SSB_high,
      SSB_fwd          = ssb(stk_fwd)[, tac.yr.name, , spawn.season],
      res_catch        = res_catch,
      closed_iter      = closed_iter,
      fcap_iter        = fcap_iter,
      search_iter      = search_iter,
      n_search         = n_search,
      sr               = sr,
      df_final         = df_final
    ))
  }
  
  return(advice)
}
environment(FcapBpaHCR_ane) <- asNamespace("FLBEIA")
assignInNamespace(x     = "annualTAC",value = FcapBpaHCR_ane, ns    = "FLBEIA")

