# ============================================================
# 06_fit_SR_bootstrap.R
# Stock-recruitment conditioning across bootstrap OMs
#
# Stock / analysis:
#   European anchovy, ane.27.9aS (Gulf of Cádiz)
#   Beverton-Holt stock-recruitment conditioning across the
#   100 bootstrap-conditioned historical operating models.
#
# Role in the workflow:
#   This script follows historical OM conditioning. It extracts
#   Q2 SSB and Q3 recruitment from each of the 100 conditioned OMs
#   and fits one stock-recruitment relationship per OM. The positive
#   Beverton-Holt summaries are used by the subsequent SR-reference
#   construction step.
#
# Purpose:
#   - Extract historical SSB and recruitment for 1989-2024.
#   - Preserve the OM iteration to SS3-bootstrap correspondence.
#   - Fit unconstrained Beverton-Holt relationships as diagnostics.
#   - Fit positive Beverton-Holt relationships (a > 0, b > 0).
#   - Estimate OM-specific recruitment residual variability (sigmaR).
#   - Compare unconstrained and positive fits.
#
# Historical period:
#   1989-2024
#
# Biological timing:
#   SSB         : Q2
#   Recruitment : Q3
#
# Inputs:
#   data/Rdata/ane_stock_conditioned_100iter.rds
#   data/Rdata/iteration_map.csv
#
# Outputs:
#   data/Rdata/SR_historical_bootstrap.rds
#   data/Rdata/fit_bh_boot_unconstrained_100.rds
#   data/Rdata/bh_summary_unconstrained.csv
#   data/Rdata/fit_bh_positive_100.rds
#   data/Rdata/bh_summary_positive.csv
#   data/Rdata/bh_free_vs_positive.csv
#
# Working Document outputs:
#   The CSV summaries provide the OM-specific BH parameters and
#   sigmaR values used to document variability among conditioned OMs.
#   This script does not generate figures.
#
# Important:
#   The positive BH fit constrains a > 0 and b > 0 through log-scale
#   parameterisation. No recruitment process error is generated here.
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(tidyr)
library(purrr)
library(here)

wd <- here()
setwd(wd)

rds_boot_dir <- file.path(getwd(), "data", "Rdata")

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

first.yr <- 1989
last.obs.yr <- 2024

ss.rec <- 3   # Recruitment enters in Q3
ss.ssb <- 2   # Spawning biomass in Q2

n_boot <- 100

# ============================================================
# 2. LOAD HISTORICAL CONDITIONING
# ============================================================

ane.stock <- readRDS(file.path(rds_boot_dir,"ane_stock_conditioned_100iter.rds"))
iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map.csv"))


# ============================================================
# 3. EXTRACT HISTORICAL SSB AND RECRUITMENT
# ============================================================

rec_h <- rec(ane.stock)
ssb_h <- ssb(ane.stock)


# ============================================================
# 4. BUILD HISTORICAL STOCK-RECRUITMENT TABLE
# ============================================================

sr_hist <- map_dfr(seq_len(n_boot),
  function(i) {
    tibble(iter      = i,
           bootstrap = iteration_map$bootstrap[i],
           year      = first.yr:last.obs.yr,
           SSB = as.numeric(ssb_h[,ac(first.yr:last.obs.yr),,ac(ss.ssb),,i]),
           R = as.numeric(rec_h[,ac(first.yr:last.obs.yr),,ac(ss.rec),,i]))
  })

# Check integrity of the historical SR table
sr_check <- sr_hist %>%
  dplyr::summarise(
    n          = dplyr::n(),
    n_boot     = dplyr::n_distinct(bootstrap),
    n_year     = dplyr::n_distinct(year),
    n_na_SSB   = sum(is.na(SSB)),
    n_na_R     = sum(is.na(R)),
    n_zero_SSB = sum(SSB <= 0, na.rm = TRUE),
    n_zero_R   = sum(R <= 0, na.rm = TRUE))

print(sr_check)

# ============================================================
# 5. CHECK HISTORICAL STOCK-RECRUITMENT DATA
# ============================================================
# The summary above checks dimensions, missing values and non-positive
# SSB or recruitment values before model fitting.

# ============================================================
# 6. FIT UNCONSTRAINED BEVERTON-HOLT RELATIONSHIPS
# ============================================================

fit_bh_boot <- vector("list", n_boot)
names(fit_bh_boot) <- iteration_map$bootstrap

bh_summary <- vector("list", n_boot)

for(i in seq_len(n_boot)) {
  
  boot_id <- iteration_map$bootstrap[i]
  
  cat("Fitting BH:",boot_id,"(",i,"/",n_boot,")\n")
  
  # Extraer una sola historia
  rec_i <- rec_h[,ac(first.yr:last.obs.yr),,ac(ss.rec),,i,drop = FALSE]
  ssb_i <- ssb_h[,ac(first.yr:last.obs.yr),,ac(ss.ssb),,i,drop = FALSE]
  # Ajuste BH
  mod_i <- FLSR(rec = rec_i,ssb = ssb_i,model = bevholt)
  fit_i <- tryCatch(fmle(mod_i),error = function(e) NULL)
  fit_bh_boot[[i]] <- fit_i
  if (!is.null(fit_i)) {res_i <- as.numeric(residuals(fit_i))
    res_i <- res_i[is.finite(res_i)]
    bh_summary[[i]] <- tibble(iter = i,
                              bootstrap = boot_id,
                              fit_ok = TRUE,
                              a = as.numeric(params(fit_i)["a", ]),
                              b = as.numeric(params(fit_i)["b", ]),
                              sigmaR = sd(res_i,na.rm = TRUE),
                              mean_residual = mean(res_i,na.rm = TRUE),
                              n = length(res_i))} 
  else {
    bh_summary[[i]] <- tibble(iter = i,
                              bootstrap = boot_id,
                              fit_ok = FALSE,
                              a = NA_real_, 
                              b = NA_real_,
                              sigmaR = NA_real_,
                              mean_residual = NA_real_,
                              n = NA_integer_)
  }
}

bh_summary <- dplyr::bind_rows(bh_summary)

# ============================================================
# 7. DIAGNOSE UNCONSTRAINED FITS
# ============================================================
# Inspect fit success, parameter distributions, residual variability
# and the occurrence of non-positive b estimates.

table(bh_summary$fit_ok)
summary(bh_summary$a)
summary(bh_summary$b)
summary(bh_summary$sigmaR)

bh_summary %>%dplyr::arrange(desc(sigmaR)) %>%
  dplyr::select(iter,bootstrap,a,b,sigmaR,mean_residual,fit_ok) %>%
  print(n = Inf)

bh_summary %>%dplyr::summarise(
    n = dplyr::n(),
    n_b_negative = sum(b <= 0),
    prop_b_negative = mean(b <= 0))


# ============================================================
# 8. SAVE UNCONSTRAINED FIT OUTPUTS
# ============================================================

saveRDS(fit_bh_boot,file.path(rds_boot_dir, "fit_bh_boot_unconstrained_100.rds"))
write.csv(bh_summary,file.path(rds_boot_dir, "bh_summary_unconstrained.csv"),row.names = FALSE)

# ============================================================
# 9. DEFINE POSITIVE BEVERTON-HOLT OBJECTIVE FUNCTION
# a > 0, b > 0
# ============================================================

bh_ssq_pos <- function(par, SSB, R) {
  a <- exp(par[1])
  b <- exp(par[2])
  Rpred <- a * SSB / (b + SSB)
  res <- log(R) - log(Rpred)
  sum(res^2)
}

# ============================================================
# 10. FIT POSITIVE BEVERTON-HOLT RELATIONSHIPS
# ============================================================

fit_bh_pos <- vector("list",n_boot)
names(fit_bh_pos) <-iteration_map$bootstrap
bh_summary_pos <- vector("list",n_boot)

for(i in seq_len(n_boot)) {
  
  boot_id <- iteration_map$bootstrap[i]
  cat("Fitting positive BH:",boot_id,"(",i,"/",n_boot,")\n")
  
  dat_i <- sr_hist %>%dplyr::filter(iter == i)
  fit_i <- optim(par = c(log(median(dat_i$R)),log(median(dat_i$SSB))),
                 fn = bh_ssq_pos,
                 SSB = dat_i$SSB,
                 R   = dat_i$R,
                 method = "Nelder-Mead",
                 control = list(maxit = 10000,reltol = 1e-12))
  
  a_i <- exp(fit_i$par[1])
  b_i <- exp(fit_i$par[2])
  Rpred_i <- a_i * dat_i$SSB /(b_i + dat_i$SSB)
  res_i <-log(dat_i$R) -log(Rpred_i)
  sigmaR_i <- sd(res_i,na.rm = TRUE)
  fit_bh_pos[[i]] <- fit_i
  bh_summary_pos[[i]] <- tibble(
    iter        = i,
    bootstrap   = boot_id,
    convergence = fit_i$convergence,
    a           = a_i,
    b           = b_i,
    sigmaR      = sigmaR_i,
    objective   = fit_i$value)
}

bh_summary_pos <- dplyr::bind_rows(bh_summary_pos)

# ============================================================
# 11. VALIDATE POSITIVE BH FITS
# ============================================================

stopifnot(nrow(bh_summary_pos) == n_boot)
stopifnot(all(bh_summary_pos$convergence == 0))
stopifnot(all(is.finite(bh_summary_pos$a)))
stopifnot(all(is.finite(bh_summary_pos$b)))
stopifnot(all(is.finite(bh_summary_pos$sigmaR)))

# ============================================================
# 12. DIAGNOSE POSITIVE FITS
# ============================================================
# Summarise optimiser convergence, fitted parameter distributions
# and the frequency of small b estimates.

table(bh_summary_pos$convergence)
summary(bh_summary_pos$a)
summary(bh_summary_pos$b)
summary(bh_summary_pos$sigmaR)

bh_summary_pos %>%dplyr::summarise(n = dplyr::n(),
                                   n_converged =sum(convergence == 0),
                                   n_b_lt_1 =sum(b < 1),
                                   n_b_lt_100 =sum(b < 100),
                                   n_b_lt_1000 =sum(b < 1000))


# ============================================================
# 13. SUMMARISE POSITIVE BH PARAMETER DISTRIBUTIONS
# ============================================================
# Calculate median and 5th-95th percentile summaries for a, b and
# sigmaR across the 100 positive BH fits.

bh_diag <- bh_summary_pos %>%
  summarise(
    n = dplyr::n(),
    a_median = median(a),
    a_p05 = quantile(a, 0.05),
    a_p95 = quantile(a, 0.95),
    b_median = median(b),
    b_p05 = quantile(b, 0.05),
    b_p95 = quantile(b, 0.95),
    sigmaR_median = median(sigmaR),
    sigmaR_p05 = quantile(sigmaR, 0.05),
    sigmaR_p95 = quantile(sigmaR, 0.95),
    n_b_near_zero = sum(b < 1)
  )

bh_diag

# ============================================================
# 14. COMPARE UNCONSTRAINED AND POSITIVE FITS
# ============================================================
# Join both summaries by bootstrap identifier for direct comparison.

bh_compare <- bh_summary %>%
  dplyr::select(bootstrap,a_free = a,b_free = b,sigmaR_free = sigmaR) %>%
  dplyr::left_join(
    bh_summary_pos %>%
      dplyr::select(bootstrap,convergence,a_pos = a,b_pos = b,sigmaR_pos = sigmaR),by = "bootstrap")

bh_compare %>%dplyr::filter(b_free <= 0)

# ============================================================
# 15. SAVE POSITIVE FIT AND COMPARISON OUTPUTS
# ============================================================

saveRDS(fit_bh_pos,file.path(rds_boot_dir,"fit_bh_positive_100.rds"))
write.csv(bh_summary_pos,file.path(rds_boot_dir,"bh_summary_positive.csv"),row.names = FALSE)
write.csv(bh_compare,file.path(rds_boot_dir,"bh_free_vs_positive.csv"),row.names = FALSE)
saveRDS(sr_hist,file.path(rds_boot_dir,"SR_historical_bootstrap.rds"))
write.csv(bh_diag,file.path(rds_boot_dir,"bh_positive_diagnostics.csv"),row.names = FALSE)
