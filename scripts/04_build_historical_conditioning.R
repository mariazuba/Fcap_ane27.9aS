# ============================================================
# 04_build_historical_conditioning.R
# Build historical conditioning across 100 OM iterations
#
# Purpose:
#   Combine the 100 validated SS3 bootstrap FLStock objects
#   into a single FLStock using the iteration dimension.
#
#   The script:
#     - Loads the 100 validated historical FLStock conditions.
#     - Combines them into a single 100-iteration FLStock.
#     - Verifies that each iteration matches its source bootstrap.
#     - Checks consistency among stock, catch and landings weights.
#     - Saves the conditioned historical stock for OM construction.
#
# INPUTS:
#   data/Rdata/
#     - stk_boot_100.rds
#     - iteration_map.csv
#
# OUTPUTS:
#   data/Rdata/
#     - stk_conditioned_100iter.rds
#     - iteration_map.csv
#
# Notes:
#   - Each iteration corresponds to one validated SS3 bootstrap.
#   - The bootstrap-to-iteration correspondence is preserved
#     through iteration_map.csv.
#   - Number of historical conditions: 100.
# ============================================================

rm(list = ls())

# ============================================================
# LOAD LIBRARIES
# ============================================================

library(FLBEIA)
library(ggplotFL)
library(r4ss)
library(icesTAF)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(callr)
library(here)

# ============================================================
# LOAD 100 HISTORICAL CONDITIONS
# ============================================================

n_boot <- 100
ni <- n_boot

rds_boot_dir <- file.path(getwd(),"data","Rdata")
stk_boot <- readRDS(file.path(rds_boot_dir,"stk_boot_100.rds"))
iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map.csv"))

# ============================================================
# CHECK INPUTS
# ============================================================

stopifnot(length(stk_boot) == n_boot,nrow(iteration_map) == n_boot)
stopifnot(identical(names(stk_boot),iteration_map$bootstrap))

names(stk_boot)[1:10]
head(iteration_map)

# ============================================================
# COMBINE 100 FLStock INTO ITER DIMENSION
# ============================================================

fill_FLStock_iterations <- function(stk_list) {
  
  n <- length(stk_list)
  out <- propagate(stk_list[[1]],n)
  
  for(i in seq_len(n)) {
    
    x <- stk_list[[i]]
    stock.n(out)[,,,,,i]      <- stock.n(x)[,,,,,1]
    stock.wt(out)[,,,,,i]     <- stock.wt(x)[,,,,,1]
    m(out)[,,,,,i]            <- m(x)[,,,,,1]
    mat(out)[,,,,,i]          <- mat(x)[,,,,,1]
    harvest(out)[,,,,,i]      <- harvest(x)[,,,,,1]
    harvest.spwn(out)[,,,,,i] <- harvest.spwn(x)[,,,,,1]
    m.spwn(out)[,,,,,i]       <- m.spwn(x)[,,,,,1]
    catch.n(out)[,,,,,i]      <- catch.n(x)[,,,,,1]
    catch.wt(out)[,,,,,i]     <- catch.wt(x)[,,,,,1]
    landings.n(out)[,,,,,i]   <- landings.n(x)[,,,,,1]
    landings.wt(out)[,,,,,i]  <- landings.wt(x)[,,,,,1]
    discards.n(out)[,,,,,i]   <- discards.n(x)[,,,,,1]
    discards.wt(out)[,,,,,i]  <- discards.wt(x)[,,,,,1]
  }
  
  out
}

stk <- fill_FLStock_iterations(
  stk_boot
)

# ============================================================
# BASIC CHECKS
# ============================================================

stopifnot(length(dimnames(stk)$iter) == n_boot)
range(stk)
dimnames(stk)$iter[1:10]

# ============================================================
# VERIFY THAT ITERATIONS MATCH ORIGINAL FLStock
# ============================================================
check_iter <- function(i) {
  
  boot_id <- iteration_map$bootstrap[i]
  
  tibble(
    iter = i,
    bootstrap = boot_id,
    stock_n = isTRUE(all.equal(as.numeric(stock.n(stk)[,,,,,i]),as.numeric(stock.n(stk_boot[[boot_id]])[,,,,,1]))),
    stock_wt = isTRUE(all.equal(as.numeric(stock.wt(stk)[,,,,,i]),as.numeric(stock.wt(stk_boot[[boot_id]])[,,,,,1]))),
    harvest = isTRUE(all.equal(as.numeric(harvest(stk)[,,,,,i]),as.numeric(harvest(stk_boot[[boot_id]])[,,,,,1]))),
    catch_n = isTRUE(all.equal(as.numeric(catch.n(stk)[,,,,,i]),as.numeric(catch.n(stk_boot[[boot_id]])[,,,,,1]))),
    catch_wt = isTRUE(all.equal(as.numeric(catch.wt(stk)[,,,,,i]),as.numeric(catch.wt(stk_boot[[boot_id]])[,,,,,1]))),
    landings_n = isTRUE(all.equal(as.numeric(landings.n(stk)[,,,,,i]),as.numeric(landings.n(stk_boot[[boot_id]])[,,,,,1]))),
    landings_wt = isTRUE(all.equal(as.numeric(landings.wt(stk)[,,,,,i]),as.numeric(landings.wt(stk_boot[[boot_id]])[,,,,,1])))
  )
}

check_all <- map_dfr(seq_len(n_boot),check_iter)
print(check_all)

check_summary <- check_all %>% summarise(across(c(stock_n,stock_wt,harvest,catch_n,catch_wt,landings_n,landings_wt),all))
print(check_summary)

stopifnot(all(unlist(check_summary)))

boot_wt_check <- map_dfr(seq_along(stk_boot),function(i) {
  
  x <- stk_boot[[i]]
  d_catch <- as.numeric(catch.wt(x) - stock.wt(x))
  d_land <- as.numeric(landings.wt(x) - stock.wt(x))
  
  tibble(
    iter = i,
    bootstrap = names(stk_boot)[i],
    maxdiff_stock_catch = max(abs(d_catch),na.rm=TRUE),
    maxdiff_stock_landings = max(abs(d_land),na.rm=TRUE),
    stock_catch_equal = all(abs(d_catch) < 1e-10,na.rm=TRUE),
    stock_landings_equal = all(abs(d_land) < 1e-10,na.rm=TRUE)
  )
})

print(boot_wt_check,n=Inf)

boot_wt_summary <- boot_wt_check %>% summarise(n_boot=n(),stock_catch_equal=sum(stock_catch_equal),stock_landings_equal=sum(stock_landings_equal),maxdiff_stock_catch=max(maxdiff_stock_catch),maxdiff_stock_landings=max(maxdiff_stock_landings))
print(boot_wt_summary)

# ============================================================
# SAVE COMBINED HISTORICAL STOCK
# ============================================================

saveRDS(stk,file.path(rds_boot_dir,"stk_conditioned_100iter.rds"))
write.csv(iteration_map,file.path(rds_boot_dir,"iteration_map.csv"),row.names = FALSE)

