##==============================================================================
## RESUMEN MÍNIMO PARA ESCENARIOS FCAP
##==============================================================================

rm(list = ls())
graphics.off()

library(FLBEIA)
library(dplyr)
library(tidyr)
library(here)

wd <- here()
setwd(wd)

dir.create("output", showWarnings = FALSE)

Blim <- 4721
Bpa  <- 6561

get_iters      <- function(obj) {
  as.integer(dim(obj$biols[[1]]@n)[6])
}
proj_year_from <- function(obj) {
  y0 <- obj$main.ctrl$sim.years[["initial"]]
  as.integer(y0 - 1)
}
safe_name      <- function(x) {
  gsub("^flbeia_|\\.rds$", "", x)
}
bio_sum_one    <- function(obj, sc_name) {
  bioSum(obj,scenario = sc_name,ssb_season = 2,long = TRUE)
}
risk_sum_one <- function(obj, sc_name, stock_name = "ANE", ssb_season = "2") {
  biol_i <- obj$biols[[stock_name]]
  ssb_all <- quantSums(biol_i@n *biol_i@wt *predict(biol_i@mat))
  ssb_s2 <- ssb_all[, , , ssb_season, , ]
  ssb_df <- as.data.frame(ssb_s2)
  names(ssb_df)[names(ssb_df) == "data"] <- "SSB"
  ssb_df %>% mutate(
    scenario = sc_name,
    unit = stock_name,
    year = as.numeric(as.character(year)) ) %>%
    group_by(year, unit, scenario) %>%
    summarise(
      pBlim = mean(SSB < Blim, na.rm = TRUE),
      pBpa  = mean(SSB < Bpa,  na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(
      cols = c(pBlim, pBpa),
      names_to = "indicator",
      values_to = "value"
    )
}

mod_files <- list.files("model",pattern = "^flbeia_Fcap_.*\\.rds$",full.names = TRUE)

if (!length(mod_files)) {stop("No hay archivos model/flbeia_Fcap_*.rds")}

tmp_obj <- readRDS(mod_files[1])
NITER   <- get_iters(tmp_obj)
proj.yr <- proj_year_from(tmp_obj)
rm(tmp_obj)

bio_all  <- list()
risk_all <- list()

for (f in mod_files) {
  
  sc_name <- safe_name(basename(f))
  message("Cargando ", sc_name)
  
  obj <- readRDS(f)
  
  bio_all[[sc_name]] <- bio_sum_one(obj, sc_name)
  
  risk_all[[sc_name]] <- tryCatch(
    risk_sum_one(obj, sc_name),
    error = function(e) NULL)
  
  rm(obj)
  gc()
}

bio_all  <- bind_rows(bio_all)
bioQ_all <- bioSumQ(bio_all)
risk_all <- bind_rows(Filter(Negate(is.null), risk_all))

summary_rds <- file.path("output",paste0("summary_flbeia_Fcap_it", NITER, ".rds"))

saveRDS(list(
  bio  = bio_all,
  bioQ = bioQ_all,
  risk = risk_all,
  meta = list(
    created   = Sys.time(),
    n_models  = length(mod_files),
    files     = basename(mod_files),
    iters     = NITER,
    proj_year = proj.yr,
    refs      = list(Blim = Blim, Bpa = Bpa))),
  summary_rds)

message("Resumen guardado en: ", summary_rds)

