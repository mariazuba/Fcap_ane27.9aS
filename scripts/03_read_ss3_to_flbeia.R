# ============================================================
# 03_read_ss3_to_flbeia.R
# Convert SS3 output to FLBEIA objects
#
# Purpose:
#   Read an executed SS3 model and convert its outputs into
#   FLR/FLBEIA objects for operating-model conditioning.
#
#   The function:
#     - Reads SS3 model output.
#     - Creates FLStock, FLSR and FLIndex objects.
#     - Replaces FLStock harvest with SS3 F-at-age.
#     - Reconstructs survey selectivity patterns from SS3.
#     - Sets the Fbar age range to age 3.
#     - Saves the resulting objects as RDS files.
#
# REQUIRED PACKAGES:
#   - r4ss
#   - ss3om
#   - FLCore
#   - dplyr
#   - tidyr
#
# FUNCTION:
#   read_ss3_to_flbeia(wd,output_name,rds_dir)
#
# ARGUMENTS:
#   wd          Directory containing an executed SS3 model.
#   output_name Identifier used for output filenames.
#   rds_dir     Directory where converted objects are saved.
#
# INPUT:
#   A complete executed SS3 model directory containing the
#   input and output files required by r4ss and ss3om.
#
# OUTPUTS:
#   ss3_<output_name>.rds   - SS3 output
#   stk_<output_name>.rds   - FLStock
#   sr_<output_name>.rds    - FLSR
#   idxs_<output_name>.rds  - FLIndex objects
#
# Returns invisibly:
#   list(output,stk,sr,idxs)
# ============================================================
`%||%` <- function(x, y) if (is.null(x)) y else x


read_ss3_to_flbeia <- function(wd,output_name,rds_dir) {
                      stopifnot(dir.exists(wd))
                      if (!dir.exists(rds_dir)) {dir.create(rds_dir,recursive = TRUE,showWarnings = FALSE)
  }
  
  # ==========================================================
  # Helper para selectividad de índices
  # ==========================================================
  
  build_selpattern_from_ss3 <- function(index_obj, sel_df) {
    
    valid_pairs <- as.data.frame(index_obj@index) %>%
      dplyr::filter(iter == 1,unit == "unique",area == "unique",!is.na(data)) %>%
      dplyr::transmute(Yr = as.integer(as.character(year)),
                       Seas = as.integer(as.character(season))) %>%
      dplyr::distinct()
    
    age_cols <- grep("^\\d+$",names(sel_df),value = TRUE)
    ages_in <- sort(as.integer(age_cols))
    
    sel_df2 <- sel_df %>%dplyr::left_join(valid_pairs %>%
          dplyr::mutate(valid = TRUE),by = c("Yr", "Seas")) %>%
          dplyr::mutate(dplyr::across(dplyr::all_of(age_cols),~ ifelse(!is.na(valid),.x,NA_real_))) %>%
          dplyr::select(-valid)
    
    dn_sp <- tryCatch(dimnames(index_obj@sel.pattern),error = function(e) NULL)
    dn_idx <- dimnames(index_obj@index)
    ages_tgt <- if (!is.null(dn_sp)) {as.integer(dn_sp$age)} else {ages_in}
    years_tgt <- if (!is.null(dn_sp)) {as.integer(dn_sp$year)} else {as.integer(dn_idx$year)}
    units_tgt <- (if (!is.null(dn_sp)) dn_sp$unit else dn_idx$unit) %||% "1"
    seasons_tgt <- if (!is.null(dn_sp)) {as.integer(dn_sp$season)} else {as.integer(dn_idx$season)}
    area_tgt <- (if (!is.null(dn_sp)) dn_sp$area else dn_idx$area) %||% "unique"
    iters_tgt <- (if (!is.null(dn_sp)) dn_sp$iter else dn_idx$iter) %||% "1"
    df_long <- sel_df2 %>% dplyr::select(Yr,Seas,dplyr::all_of(as.character(ages_in))) %>%
      tidyr::pivot_longer(dplyr::all_of(as.character(ages_in)),names_to = "age",values_to = "val") %>%
      dplyr::mutate(age = as.integer(age),Yr = as.integer(Yr),Seas = as.integer(Seas))
    arr <- array(NA_real_,
      dim = c(
        length(ages_tgt),
        length(years_tgt),
        length(units_tgt),
        length(seasons_tgt),
        length(area_tgt),
        length(iters_tgt)
      ),
      dimnames = list(
        age = as.character(ages_tgt),
        year = as.character(years_tgt),
        unit = units_tgt,
        season = as.character(seasons_tgt),
        area = area_tgt,
        iter = iters_tgt))
    
    a <- match(df_long$age, ages_tgt)
    y <- match(df_long$Yr, years_tgt)
    s <- match(df_long$Seas, seasons_tgt)
    
    keep <- !is.na(a) &
      !is.na(y) &
      !is.na(s)
    
    if (any(keep)) {arr[cbind(a[keep],y[keep],1,s[keep],1,1)] <- df_long$val[keep]}
    
    FLQuant(arr)
  }
  
  
  # ==========================================================
  # LEER SS3 YA EJECUTADO
  # ==========================================================
  
  message("Reading ",basename(wd))
  
  output <- r4ss::SS_output(dir = wd,forecast = FALSE,covar = FALSE,verbose = FALSE)
  sr <- ss3om::readFLSRss3(dir = wd)
  idxs <- ss3om::readFLIBss3(dir = wd)
  stk <- ss3om::readFLSss3(dir = wd)
  
  # ==========================================================
  # REEMPLAZAR harvest(stk) POR F-AT-AGE SS3
  # ==========================================================
  
  fatage_ss3 <- subset(output$fatage,Era == "TIME")
  age_cols <- intersect(colnames(fatage_ss3),as.character(dimnames(harvest(stk))$age))
  fatage_ss3 <- fatage_ss3[rowSums(fatage_ss3[,age_cols,drop = FALSE]) != 0,c("Yr","Seas",age_cols)]
  fq <- harvest(stk)
  dn <- dimnames(fq)
  fat_long <- do.call(rbind,lapply(age_cols,
      function(a) {
        data.frame(year =as.character(fatage_ss3$Yr),
          season =as.character(fatage_ss3$Seas),age = a,F =fatage_ss3[[a]],stringsAsFactors = FALSE)
      }
    ))
  
  fat_ok <- subset(fat_long,year %in% dn$year &season %in% dn$season & age %in% dn$age)
  newF <- fq
  newF[] <- NA_real_
  
  if (nrow(fat_ok)) {
    a <- match(fat_ok$age,dn$age)
    y <- match(fat_ok$year,dn$year)
    s <- match(fat_ok$season,dn$season)
    newF[cbind(a,y,1,s,1,1)] <- fat_ok$F
  }
  
  harvest(stk) <- newF
  
  
  # ==========================================================
  # SELECTIVIDAD DE CAMPAÑAS
  # ==========================================================
  
  if ("ageselex" %in% names(output)) {
    
    fleets_map <- list(PELAGO = 5,ECOCADIZ = 6,BOCADEVA = 7,ECORECLUTAS = 8)
    
    for (nm in names(fleets_map)) {
      if (nm %in% names(idxs)) {
        sel_df <- subset(output$ageselex,Fleet == fleets_map[[nm]] &Factor == "Asel2")
        idxs[[nm]]@sel.pattern <-build_selpattern_from_ss3(idxs[[nm]],sel_df)
      }}}
  
  # ==========================================================
  # Fbar
  # ==========================================================
  
  range(stk, "minfbar") <- 3
  range(stk, "maxfbar") <- 3
  
  # ==========================================================
  # GUARDAR
  # ==========================================================
  
  saveRDS(output,file.path(rds_dir,paste0("ss3_",output_name, ".rds")))
  saveRDS(stk,file.path(rds_dir,paste0("stk_",output_name,".rds")))
  saveRDS(sr,file.path(rds_dir,paste0("sr_",output_name,".rds")))
  saveRDS(idxs,file.path(rds_dir,paste0("idxs_",output_name,".rds")))
  invisible(list(output = output,stk = stk,sr = sr,idxs = idxs))
}

