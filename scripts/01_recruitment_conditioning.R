# ============================================================
# 01_recruitment_conditioning.R
# Recruitment conditioning diagnostics for the Operating Model
# Anchovy 9a South
#
# Purpose:
#   1. Compare the OM Beverton-Holt SR relationship with SS3.
#   2. Diagnose recruitment residuals.
#   3. Produce figures and tables for the MSE Working Document.
#
# Outputs:
#   outputs/recruitment_conditioning/
#     - tables/*.csv
#     - figures/*.png
#   data/processed/recruitment_conditioning_OM_vs_SS3.RData
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "FLBEIA", "FLCore", "dplyr", "tibble", "ggplot2",
  "stringr", "here", "readr",
  "changepoint", "strucchange", "broom",
  "zoo", "randtests"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),"\nInstall them before running this script.")
}

library(FLBEIA)
library(FLCore)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(here)
library(readr)
library(zoo)
library(randtests)

set.seed(123)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

first_yr    <- 1989
last_obs_yr <- 2024

# Seasonal timing used for the SR relationship
# Recruitment occurs in Q3; SSB is evaluated in Q2
ss_rec <- 3
ss_ssb <- 2

# Reference points
Blim <- 4721
Bpa  <- 6561

# Input files
stk_ane9aS_rds <- here("boot/data/stk_ane9aS.rds")
ss3_rds        <- here("boot/data/ss3_ane9aS.rds")

# Output folders
out_dir  <- here("outputs/recruitment_conditioning")
fig_dir  <- file.path(out_dir, "figures")
tab_dir  <- file.path(out_dir, "tables")
data_dir <- here("data/processed")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

check_file_exists <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
}

find_value <- function(df, patterns, label_col = "Label", value_col = "Value") {
  if (!all(c(label_col, value_col) %in% names(df))) {
    stop("Columns ", label_col, " and/or ", value_col, " not found.")
  }
  pattern <- paste(patterns, collapse = "|")
  out <- df %>%
    filter(str_detect(.data[[label_col]], regex(pattern, ignore_case = TRUE)))
  if (nrow(out) == 0) return(NA_real_)
  as.numeric(out[[value_col]][1])
}

ss3_bh <- function(SSB, R0, h, SSB0) {
  (4 * h * R0 * SSB) /(SSB0 * (1 - h) + SSB * (5 * h - 1))
}

om_bh <- function(SSB, a, b) {a * SSB / (b + SSB)}

make_lognormal_bands <- function(rec_med, sigmaR) {
  tibble(
    rec_lower_95 = rec_med * exp(-0.5 * sigmaR^2 - 2.00 * sigmaR),
    rec_upper_95 = rec_med * exp(-0.5 * sigmaR^2 + 2.00 * sigmaR),
    rec_lower_80 = rec_med * exp(-0.5 * sigmaR^2 - 1.28 * sigmaR),
    rec_upper_80 = rec_med * exp(-0.5 * sigmaR^2 + 1.28 * sigmaR)
  )
}

save_plot <- function(plot, filename, width = 8, height = 5) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

safe_box_test <- function(x, lag) {
  out <- Box.test(x, lag = lag, type = "Ljung-Box")
  tibble(
    test = paste0("Ljung-Box lag ", lag),
    statistic = unname(out$statistic),
    p_value = out$p.value
  )
}

calc_skewness <- function(x) {
  x <- x[is.finite(x)]
  mean((x - mean(x))^3) / sd(x)^3}

calc_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  mean((x - mean(x))^4) / sd(x)^4}

# ------------------------------------------------------------
# 3. Load stock and define biological timing
# ------------------------------------------------------------

check_file_exists(stk_ane9aS_rds)
check_file_exists(ss3_rds)

ane_stock <- readRDS(stk_ane9aS_rds)

name(ane_stock) <- "ANE"
desc(ane_stock) <- "WGHANSA2025 assessment output"
ane_stock@harvest@units <- "f"

# Fishing and natural mortality timing
fspwn <- harvest.spwn(ane_stock)
fspwn[] <- 0
fspwn[, , , "1", ] <- 1
harvest.spwn(ane_stock) <- fspwn

mspwn <- m.spwn(ane_stock)
mspwn[] <- 0
mspwn[, , , "1", ] <- 1
m.spwn(ane_stock) <- mspwn

# Maturity timing: spawning in Q2, age 0 not mature
mat(ane_stock)[ac(0), ] <- 0
mat(ane_stock)[, , , 1] <- 0
mat(ane_stock)[, , , 3] <- 0
mat(ane_stock)[, , , 4] <- 0

range(ane_stock, "minfbar") <- 3
range(ane_stock, "maxfbar") <- 3

# Create FLBiol only for historical period
ane_biol <- FLBiol(
  n = stock.n(ane_stock),
  wt = stock.wt(ane_stock),
  m = m(ane_stock),
  spwn = m.spwn(ane_stock),
  mat = mat(ane_stock),
  fec = predictModel(
    FLQuants(fec = ane_stock@mat * 0 + 1),
    model = ~fec
  ),
  name = "ANE",
  desc = ane_stock@desc,
  range = ane_stock@range
)

units(ane_biol)$m <- units(fec(ane_biol)) <- units(mat(ane_biol)) <- ""
ane_biol <- window(ane_biol, start = first_yr, end = last_obs_yr)

# ------------------------------------------------------------
# 4. Extract historical recruitment and SSB
# ------------------------------------------------------------

rec_h <- rec(ane_biol)[, ac(first_yr:last_obs_yr), , ss_rec, drop = FALSE]
ssb_h <- ssb(ane_biol)[, ac(first_yr:last_obs_yr), , ss_ssb, drop = FALSE]

# Use one iteration only if the object contains iterations
rec_fit <- iter(rec_h, 1)
ssb_fit <- iter(ssb_h, 1)

sr_data <- tibble(
  year = as.numeric(dimnames(rec_fit)$year),
  ssb = as.numeric(ssb_fit),
  rec = as.numeric(rec_fit)
) %>%
  filter(is.finite(ssb), is.finite(rec), ssb > 0, rec > 0)

if (nrow(sr_data) != length(first_yr:last_obs_yr)) {
  warning(
    "Expected ", length(first_yr:last_obs_yr),
    " years but found ", nrow(sr_data),
    " valid SR observations."
  )
}

write_csv(sr_data, file.path(tab_dir, "00_historical_ssb_recruitment.csv"))

# ------------------------------------------------------------
# 5. Fit OM Beverton-Holt SR relationship
# ------------------------------------------------------------

mod_bh_om <- FLSR(rec = rec_fit, ssb = ssb_fit, model = bevholt)
fit_bh_om <- fmle(mod_bh_om)

hist_resid_om <- as.numeric(residuals(fit_bh_om))
hist_resid_om <- hist_resid_om[is.finite(hist_resid_om)]

sigmaR_om <- sd(hist_resid_om, na.rm = TRUE)

om_pars <- as.numeric(params(fit_bh_om))
names(om_pars) <- dimnames(params(fit_bh_om))$params

a_om <- om_pars["a"]
b_om <- om_pars["b"]

if (!is.finite(a_om) || !is.finite(b_om)) {
  stop("OM parameters a or b are not finite. Check params(fit_bh_om).")
}

# Segmented regresion

mod_seg <- FLSR(
  rec = rec_fit,
  ssb = ssb_fit,
  model = "segreg"
)

fit_seg <- fmle(mod_seg)

seg_pars <- as.numeric(params(fit_seg))
names(seg_pars) <- dimnames(params(fit_seg))$params

sr_model_comparison <- tibble(
  model = c("Beverton-Holt", "Segmented"),
  logLik = c(as.numeric(logLik(fit_bh_om)), as.numeric(logLik(fit_seg))),
  AIC = c(AIC(fit_bh_om), AIC(fit_seg)),
  delta_AIC = AIC - min(AIC, na.rm = TRUE)
)

write_csv(
  sr_model_comparison,
  file.path(tab_dir, "16_SR_model_comparison_BH_segmented.csv")
)
# ------------------------------------------------------------
# 6. Load SS3 SR parameters
# ------------------------------------------------------------

ss3rep <- readRDS(ss3_rds)
pars <- ss3rep$parameters
dq   <- ss3rep$derived_quants

lnR0_ss3 <- find_value(pars,
  c("SR_LN\\(R0\\)", "SR_Ln\\(R0\\)", "SR_log\\(R0\\)", "SR_LN_R0", "SR_lnR0"))

R0_ss3 <- exp(lnR0_ss3)

h_ss3 <- find_value(pars, c("SR_BH_steep", "SR_steep", "steep"))

sigmaR_ss3 <- find_value(pars,
  c("SR_sigmaR", "SR_sigma", "sigmaR", "SigmaR", "rec_sigma", "ln\\(R\\)_sigma"))

if (!is.finite(sigmaR_ss3)) {
  sigmaR_ss3 <- sigmaR_om
  sigmaR_source <- "OM residual SD"
} else {
  sigmaR_source <- "SS3"
}

SSB0_ss3 <- find_value(
  dq,
  c("SSB_unfished", "SSB_Virgin", "SSB_virgin", "SSB_Bzero", "SSB_0", "SSB_unf")
)

if (!is.finite(R0_ss3)) stop("R0 not found in SS3 parameters.")
if (!is.finite(h_ss3)) stop("Steepness not found in SS3 parameters.")
if (!is.finite(SSB0_ss3)) stop("SSB0 not found in SS3 derived quantities.")

sr_parameter_table <- tibble(
  parameter = c("a_OM", "b_OM", "R0_SS3", "h_SS3", "SSB0_SS3", "sigmaR_SS3_or_used", "sigmaR_OM_resid"),
  value = c(a_om, b_om, R0_ss3, h_ss3, SSB0_ss3, sigmaR_ss3, sigmaR_om),
  source = c("FLSR OM fit", "FLSR OM fit", "SS3", "SS3", "SS3", sigmaR_source, "OM residuals")
)

write_csv(sr_parameter_table, file.path(tab_dir, "01_SR_parameter_table.csv"))

# ------------------------------------------------------------
# 7. OM vs SS3 comparison at observed SSB values
# ------------------------------------------------------------

sr_obs <- sr_data %>%
  mutate(
    rec_pred_om  = om_bh(ssb, a_om, b_om),
    rec_pred_ss3 = ss3_bh(ssb, R0_ss3, h_ss3, SSB0_ss3),
    resid_om     = log(rec) - log(rec_pred_om),
    resid_ss3    = log(rec) - log(rec_pred_ss3),
    ratio_om_ss3 = rec_pred_om / rec_pred_ss3,
    diff_percent = 100 * (rec_pred_om - rec_pred_ss3) / rec_pred_ss3
  )

compare_summary <- sr_obs %>%
  summarise(
    n_years = n(),
    mean_ratio_om_ss3 = mean(ratio_om_ss3, na.rm = TRUE),
    min_ratio_om_ss3 = min(ratio_om_ss3, na.rm = TRUE),
    max_ratio_om_ss3 = max(ratio_om_ss3, na.rm = TRUE),
    mean_diff_percent = mean(diff_percent, na.rm = TRUE),
    mean_abs_percent_diff = mean(abs(diff_percent), na.rm = TRUE),
    cor_residuals_om_ss3 = cor(resid_om, resid_ss3, use = "complete.obs"),
    mean_abs_residual_difference = mean(abs(resid_om - resid_ss3), na.rm = TRUE)
  )

write_csv(sr_obs, file.path(tab_dir, "02_SR_observed_OM_vs_SS3.csv"))
write_csv(compare_summary, file.path(tab_dir, "03_OM_vs_SS3_summary.csv"))

# ------------------------------------------------------------
# 8. Recruitment residual diagnostics
# ------------------------------------------------------------

rec_diag <- sr_obs %>%
  transmute(
    year,
    ssb,
    rec_obs = rec,
    rec_pred_om,
    rec_pred_ss3,
    rec_resid_om = resid_om,
    rec_resid_ss3 = resid_ss3,
    rec_mult_om = exp(resid_om),
    rec_resid_std = as.numeric(scale(resid_om)),
    resid_percentile = percent_rank(resid_om),
    residual_class = case_when(
      resid_percentile <= 0.10 ~ "lower 10%",
      resid_percentile >= 0.90 ~ "upper 10%",
      TRUE ~ "central"
    ),
    sigmaR_lower_95 = -0.5 * sigmaR_ss3^2 - 2.00 * sigmaR_ss3,
    sigmaR_upper_95 = -0.5 * sigmaR_ss3^2 + 2.00 * sigmaR_ss3,
    sigmaR_lower_80 = -0.5 * sigmaR_ss3^2 - 1.28 * sigmaR_ss3,
    sigmaR_upper_80 = -0.5 * sigmaR_ss3^2 + 1.28 * sigmaR_ss3,
    below_sigmaR_95 = rec_resid_om < sigmaR_lower_95,
    above_sigmaR_95 = rec_resid_om > sigmaR_upper_95,
    inside_sigmaR_95 = !below_sigmaR_95 & !above_sigmaR_95
  )

resid_summary <- rec_diag %>%
  summarise(
    n_years = n(),
    mean_resid = mean(rec_resid_om, na.rm = TRUE),
    sd_resid = sd(rec_resid_om, na.rm = TRUE),
    sigmaR_used = sigmaR_ss3,
    sigmaR_source = sigmaR_source,
    skewness = calc_skewness(rec_resid_om),
    kurtosis = calc_kurtosis(rec_resid_om),
    q01 = quantile(rec_resid_om, 0.01, na.rm = TRUE),
    q05 = quantile(rec_resid_om, 0.05, na.rm = TRUE),
    q10 = quantile(rec_resid_om, 0.10, na.rm = TRUE),
    median = median(rec_resid_om, na.rm = TRUE),
    q90 = quantile(rec_resid_om, 0.90, na.rm = TRUE),
    q95 = quantile(rec_resid_om, 0.95, na.rm = TRUE),
    q99 = quantile(rec_resid_om, 0.99, na.rm = TRUE),
    min_resid = min(rec_resid_om, na.rm = TRUE),
    max_resid = max(rec_resid_om, na.rm = TRUE),
    prop_inside_sigmaR_95 = mean(inside_sigmaR_95, na.rm = TRUE),
    prop_below_sigmaR_95 = mean(below_sigmaR_95, na.rm = TRUE),
    prop_above_sigmaR_95 = mean(above_sigmaR_95, na.rm = TRUE)
  )

resid_tests <- bind_rows(
  safe_box_test(rec_diag$rec_resid_om, lag = 1),
  safe_box_test(rec_diag$rec_resid_om, lag = 3),
  safe_box_test(rec_diag$rec_resid_om, lag = 5),
  {
    sw <- shapiro.test(rec_diag$rec_resid_om)
    tibble(test = "Shapiro-Wilk normality", statistic = unname(sw$statistic), p_value = sw$p.value)
  }
)

tail_diagnostics <- tibble(
  probability = c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99),
  observed_residual = as.numeric(quantile(rec_diag$rec_resid_om, probability, na.rm = TRUE)),
  normal_expected = qnorm(probability, mean = -0.5 * sigmaR_ss3^2, sd = sigmaR_ss3)
) %>%
  mutate(difference_observed_minus_expected = observed_residual - normal_expected)

extreme_residual_years <- rec_diag %>%
  filter(residual_class != "central") %>%
  arrange(rec_resid_om) %>%
  select(
    year, ssb, rec_obs, rec_pred_om, rec_pred_ss3,
    rec_resid_om, rec_mult_om, resid_percentile,
    residual_class, inside_sigmaR_95
  )

# Are poor recruitment years adequately represented?
tail_severity <- tail_diagnostics %>%
  filter(probability %in% c(0.01,0.05,0.10)) %>%
  mutate(
    severity_ratio =
      abs(observed_residual) /
      abs(normal_expected)
  )

tail_severity <- tail_severity %>%
  mutate(
    interpretation =
      case_when(
        severity_ratio > 1.25 ~
          "historical tail more extreme",
        
        severity_ratio < 0.75 ~
          "historical tail less extreme",
        
        TRUE ~
          "similar to lognormal assumption"
      )
  )

write_csv(rec_diag, file.path(tab_dir, "04_recruitment_residual_diagnostics_by_year.csv"))
write_csv(resid_summary, file.path(tab_dir, "05_recruitment_residual_summary.csv"))
write_csv(resid_tests, file.path(tab_dir, "06_recruitment_residual_tests.csv"))
write_csv(tail_diagnostics, file.path(tab_dir, "07_tail_diagnostics.csv"))
write_csv(extreme_residual_years, file.path(tab_dir, "08_extreme_residual_years.csv"))
write_csv(tail_severity,file.path(tab_dir,"18_tail_severity.csv"))

# ------------------------------------------------------------
# 8b. Autocorrelation diagnostics
# ------------------------------------------------------------

resid_vec <- rec_diag$rec_resid_om
resid_vec <- resid_vec[is.finite(resid_vec)]

acf_obj <- acf(resid_vec, plot = FALSE)
pacf_obj <- pacf(resid_vec, plot = FALSE)

acf_table <- tibble(
  lag = as.numeric(acf_obj$lag[-1]),
  acf = as.numeric(acf_obj$acf[-1]),
  ci95 = 1.96 / sqrt(length(resid_vec)),
  significant = abs(acf) > ci95
)

pacf_table <- tibble(
  lag = as.numeric(pacf_obj$lag),
  pacf = as.numeric(pacf_obj$acf),
  ci95 = 1.96 / sqrt(length(resid_vec)),
  significant = abs(pacf) > ci95
)

ar1_fit <- arima(resid_vec, order = c(1, 0, 0), include.mean = TRUE)

ar1_table <- tibble(
  parameter = names(ar1_fit$coef),
  estimate = as.numeric(ar1_fit$coef),
  se = sqrt(diag(ar1_fit$var.coef)),
  z = estimate / se,
  p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
)

write_csv(acf_table, file.path(tab_dir, "10_ACF_recruitment_residuals.csv"))
write_csv(pacf_table, file.path(tab_dir, "11_PACF_recruitment_residuals.csv"))
write_csv(ar1_table, file.path(tab_dir, "12_AR1_recruitment_residuals.csv"))

# ------------------------------------------------------------
# 8c. Regime-shift / change-point diagnostics
# ------------------------------------------------------------

log_rec <- log(sr_data$rec)
years_rec <- sr_data$year

cpt_mean <- changepoint::cpt.mean(
  log_rec,
  method = "PELT",
  penalty = "MBIC"
)

cpt_meanvar <- changepoint::cpt.meanvar(
  log_rec,
  method = "PELT",
  penalty = "MBIC"
)

cpt_table <- tibble(
  method = c("mean", "mean_and_variance"),
  n_change_points = c(
    length(changepoint::cpts(cpt_mean)),
    length(changepoint::cpts(cpt_meanvar))
  ),
  change_point_index = c(
    paste(changepoint::cpts(cpt_mean), collapse = "; "),
    paste(changepoint::cpts(cpt_meanvar), collapse = "; ")
  ),
  change_point_year = c(
    paste(years_rec[changepoint::cpts(cpt_mean)], collapse = "; "),
    paste(years_rec[changepoint::cpts(cpt_meanvar)], collapse = "; ")
  )
)

bp <- strucchange::breakpoints(log_rec ~ 1)

bp_table <- tibble(
  method = "strucchange_breakpoints",
  break_index = bp$breakpoints,
  break_year = years_rec[bp$breakpoints]
) %>%
  filter(is.finite(break_index))

write_csv(cpt_table, file.path(tab_dir, "13_changepoint_recruitment.csv"))
write_csv(bp_table, file.path(tab_dir, "14_strucchange_breakpoints_recruitment.csv"))

# ------------------------------------------------------------
# 8d. Comparison with Thorson et al. 2014
# ------------------------------------------------------------

thorson_sigmaR_mean <- 0.74
thorson_sigmaR_sd   <- 0.35
thorson_rho_mean    <- 0.43
thorson_rho_sd      <- 0.28

rho_est <- ar1_table %>%
  filter(parameter == "ar1") %>%
  pull(estimate)

literature_comparison <- tibble(
  metric = c("sigmaR", "AR1 rho"),
  ane9aS_estimate = c(sigmaR_om, rho_est),
  thorson_mean = c(thorson_sigmaR_mean, thorson_rho_mean),
  thorson_sd = c(thorson_sigmaR_sd, thorson_rho_sd),
  interpretation = c(
    ifelse(sigmaR_om < thorson_sigmaR_mean,
           "Lower than global mean recruitment variability",
           "Similar to or higher than global mean recruitment variability"),
    ifelse(abs(rho_est) < 0.2,
           "Low temporal autocorrelation",
           "Moderate/high temporal autocorrelation")
  )
)

write_csv(literature_comparison, file.path(tab_dir, "15_literature_comparison_Thorson2014.csv"))


# ------------------------------------------------------------
# 8e. Historical vs simulated recruitment residuals
# ------------------------------------------------------------
#Does the OM reproduce historical recruitment variability?
n_sim <- 1000
n_years <- length(resid_vec)

sim_resid <- replicate(
  n_sim,
  rnorm(
    n_years,
    mean = -0.5 * sigmaR_ss3^2,
    sd = sigmaR_ss3
  )
)

sim_resid_vec <- as.vector(sim_resid)

hist_vs_sim_summary <- tibble(
  statistic = c(
    "mean","sd","q01","q05",
    "median","q95","q99"
  ),
  historical = c(
    mean(resid_vec),
    sd(resid_vec),
    quantile(resid_vec,0.01),
    quantile(resid_vec,0.05),
    median(resid_vec),
    quantile(resid_vec,0.95),
    quantile(resid_vec,0.99)
  ),
  simulated = c(
    mean(sim_resid_vec),
    sd(sim_resid_vec),
    quantile(sim_resid_vec,0.01),
    quantile(sim_resid_vec,0.05),
    median(sim_resid_vec),
    quantile(sim_resid_vec,0.95),
    quantile(sim_resid_vec,0.99)
  )
)

write_csv(
  hist_vs_sim_summary,
  file.path(tab_dir,
            "17_historical_vs_simulated_residuals.csv")
)

hist_vs_sim_plot <- bind_rows(
  tibble(residual = resid_vec,
         source = "Historical"),
  tibble(residual = sim_resid_vec,
         source = "Simulated")
)

p_hist_vs_sim <- ggplot(
  hist_vs_sim_plot,
  aes(residual,
      colour = source,
      fill = source)
) +
  geom_density(alpha = 0.2) +
  theme_bw() +
  labs(
    title = "Historical vs simulated recruitment residuals"
  )

save_plot(
  p_hist_vs_sim,
  "12_historical_vs_simulated_residuals.png"
)

# ------------------------------------------------------------
# 8f. Recruitment productivity diagnostics
# ------------------------------------------------------------

sr_data <- sr_data %>%
  mutate(
    productivity = rec / ssb,
    log_productivity = log(productivity),
    productivity_roll5 = zoo::rollmean(
      productivity,
      k = 5,
      fill = NA,
      align = "center"
    ),
    period_2015 = ifelse(year < 2015, "Historical", "Recent")
  )

p_prod <- ggplot(sr_data, aes(year, productivity)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.8) +
  geom_line(aes(y = productivity_roll5), linewidth = 1.1) +
  theme_bw() +
  labs(
    title = "Recruitment productivity",
    subtitle = "Recruitment per unit of SSB, with 5-year moving average",
    x = "Year",
    y = "R / SSB"
  )

save_plot(
  p_prod,
  "13_productivity_R_per_SSB_roll5.png"
)

cp_prod <- changepoint::cpt.meanvar(
  sr_data$log_productivity,
  method = "PELT",
  penalty = "MBIC"
)

prod_cp_table <- tibble(
  variable = "log(R/SSB)",
  method = "PELT mean and variance",
  penalty = "MBIC",
  n_changepoints = length(changepoint::cpts(cp_prod)),
  change_year = ifelse(
    length(changepoint::cpts(cp_prod)) == 0,
    NA_character_,
    paste(sr_data$year[changepoint::cpts(cp_prod)], collapse = "; ")
  )
)

write_csv(
  prod_cp_table,
  file.path(tab_dir, "19_productivity_changepoints.csv")
)

period_summary <- bind_rows(
  sr_data %>%
    summarise(
      period = paste0(min(year), "-", max(year)),
      mean_R = mean(rec, na.rm = TRUE),
      geomean_R = exp(mean(log(rec), na.rm = TRUE)),
      cv_R = sd(rec, na.rm = TRUE) / mean(rec, na.rm = TRUE),
      mean_prod = mean(productivity, na.rm = TRUE),
      geomean_prod = exp(mean(log_productivity, na.rm = TRUE)),
      cv_prod = sd(productivity, na.rm = TRUE) / mean(productivity, na.rm = TRUE),
      n_years = n()
    ),
  
  sr_data %>%
    filter(year >= 2005) %>%
    summarise(
      period = "2005-present",
      mean_R = mean(rec, na.rm = TRUE),
      geomean_R = exp(mean(log(rec), na.rm = TRUE)),
      cv_R = sd(rec, na.rm = TRUE) / mean(rec, na.rm = TRUE),
      mean_prod = mean(productivity, na.rm = TRUE),
      geomean_prod = exp(mean(log_productivity, na.rm = TRUE)),
      cv_prod = sd(productivity, na.rm = TRUE) / mean(productivity, na.rm = TRUE),
      n_years = n()
    ),
  
  sr_data %>%
    filter(year >= 2015) %>%
    summarise(
      period = "2015-present",
      mean_R = mean(rec, na.rm = TRUE),
      geomean_R = exp(mean(log(rec), na.rm = TRUE)),
      cv_R = sd(rec, na.rm = TRUE) / mean(rec, na.rm = TRUE),
      mean_prod = mean(productivity, na.rm = TRUE),
      geomean_prod = exp(mean(log_productivity, na.rm = TRUE)),
      cv_prod = sd(productivity, na.rm = TRUE) / mean(productivity, na.rm = TRUE),
      n_years = n()
    )
)

write_csv(
  period_summary,
  file.path(tab_dir, "20_productivity_by_period.csv")
)

prod_ttest <- t.test(
  productivity ~ period_2015,
  data = sr_data
)

prod_wilcox <- wilcox.test(
  productivity ~ period_2015,
  data = sr_data
)

productivity_tests <- tibble(
  test = c("Welch t-test", "Wilcoxon rank-sum test"),
  statistic = c(
    unname(prod_ttest$statistic),
    unname(prod_wilcox$statistic)
  ),
  p_value = c(
    prod_ttest$p.value,
    prod_wilcox$p.value
  ),
  interpretation = ifelse(
    p_value < 0.05,
    "Evidence of difference between historical and recent productivity",
    "No evidence of difference between historical and recent productivity"
  )
)

write_csv(
  productivity_tests,
  file.path(tab_dir, "21_productivity_historical_vs_recent_tests.csv")
)

# ------------------------------------------------------------
# 8h. Density-dependence / compensation diagnostics
# Following Hilborn & Walters: examine R/S and log(R/S) vs SSB
# ------------------------------------------------------------

sr_data <- sr_data %>%
  mutate(
    RPS = rec / ssb,
    log_RPS = log(rec / ssb),
    log_SSB = log(ssb)
  )

# Linear diagnostic: decreasing log(R/S) with SSB indicates compensation
comp_lm <- lm(log_RPS ~ ssb, data = sr_data)

comp_test <- broom::tidy(comp_lm) %>%
  mutate(
    interpretation = case_when(
      term == "ssb" & estimate < 0 & p.value < 0.05 ~
        "Evidence of compensatory density dependence",
      term == "ssb" & estimate < 0 & p.value >= 0.05 ~
        "Negative slope, but weak statistical evidence",
      term == "ssb" & estimate >= 0 ~
        "No evidence of compensatory density dependence",
      TRUE ~ NA_character_
    )
  )

write_csv(
  comp_test,
  file.path(tab_dir, "27_compensation_log_RPS_vs_SSB.csv")
)

p_rps_ssb <- ggplot(sr_data, aes(ssb, RPS)) +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Recruitment per spawner vs SSB",
    subtitle = "Diagnostic for compensatory recruitment dynamics",
    x = "SSB",
    y = "Recruitment / SSB"
  )

save_plot(
  p_rps_ssb,
  "15_RPS_vs_SSB_compensation.png"
)

p_log_rps_ssb <- ggplot(sr_data, aes(ssb, log_RPS)) +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "log(R/SSB) vs SSB",
    subtitle = "Negative slope is consistent with compensation",
    x = "SSB",
    y = "log(Recruitment / SSB)"
  )

save_plot(
  p_log_rps_ssb,
  "16_log_RPS_vs_SSB_compensation.png"
)


# ------------------------------------------------------------
# 8g. Bootstrap stability of Beverton-Holt SR parameters
# ------------------------------------------------------------

nboot <- 1000
boot_pars <- vector("list", nboot)

for(i in seq_len(nboot)) {
  
  idx <- sample(seq_len(nrow(sr_data)), replace = TRUE)
  
  tmp <- sr_data[idx, ] %>%
    arrange(year) %>%
    mutate(boot_year = seq_len(n()))
  
  rec_boot <- FLQuant(
    tmp$rec,
    dimnames = list(
      age = "all",
      year = as.character(tmp$boot_year),
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )
  
  ssb_boot <- FLQuant(
    tmp$ssb,
    dimnames = list(
      age = "all",
      year = as.character(tmp$boot_year),
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )
  
  fit <- try({
    mod <- FLSR(
      rec = rec_boot,
      ssb = ssb_boot,
      model = bevholt
    )
    fmle(mod)
  }, silent = TRUE)
  
  if(!inherits(fit, "try-error")) {
    
    pars <- as.numeric(params(fit))
    names(pars) <- dimnames(params(fit))$params
    
    boot_pars[[i]] <- tibble(
      iter = i,
      a = pars["a"],
      b = pars["b"],
      logLik = as.numeric(logLik(fit)),
      AIC = AIC(fit)
    )
  }
}

boot_pars <- bind_rows(boot_pars) %>%
  filter(
    is.finite(a),
    is.finite(b),
    a > 0,
    b > 0
  )

boot_summary <- boot_pars %>%
  summarise(
    n_success = n(),
    success_rate = n() / nboot,
    a_median = median(a, na.rm = TRUE),
    a_q05 = quantile(a, 0.05, na.rm = TRUE),
    a_q95 = quantile(a, 0.95, na.rm = TRUE),
    b_median = median(b, na.rm = TRUE),
    b_q05 = quantile(b, 0.05, na.rm = TRUE),
    b_q95 = quantile(b, 0.95, na.rm = TRUE),
    b_q25 = quantile(b, 0.25, na.rm = TRUE),
    b_q75 = quantile(b, 0.75, na.rm = TRUE)
  )

write_csv(
  boot_pars,
  file.path(tab_dir, "22_bootstrap_BH_parameters.csv")
)

write_csv(
  boot_summary,
  file.path(tab_dir, "23_bootstrap_BH_summary.csv")
)

p_boot_b <- ggplot(boot_pars, aes(x = b)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = b_om, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Bootstrap distribution of Beverton-Holt b parameter",
    subtitle = paste0("Successful fits: ", nrow(boot_pars), " / ", nboot),
    x = "BH b parameter",
    y = "Frequency"
  )

save_plot(
  p_boot_b,
  "14_bootstrap_BH_b_parameter.png"
)

# ------------------------------------------------------------
# 8i. Bootstrap correlation between Beverton-Holt parameters
# ------------------------------------------------------------

boot_ab_cor <- cor.test(
  boot_pars$a,
  boot_pars$b,
  method = "spearman",
  exact = FALSE
)

boot_ab_cor_table <- tibble(
  test = "Spearman correlation between BH a and b bootstrap estimates",
  estimate = unname(boot_ab_cor$estimate),
  statistic = unname(boot_ab_cor$statistic),
  p_value = boot_ab_cor$p.value,
  interpretation = case_when(
    abs(estimate) >= 0.7 ~ "Strong parameter correlation",
    abs(estimate) >= 0.4 ~ "Moderate parameter correlation",
    TRUE ~ "Weak parameter correlation"
  )
)

write_csv(
  boot_ab_cor_table,
  file.path(tab_dir, "28_bootstrap_BH_a_b_correlation.csv")
)

p_boot_ab <- ggplot(boot_pars, aes(a, b)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_point(aes(x = a_om, y = b_om), size = 3) +
  theme_bw() +
  labs(
    title = "Bootstrap correlation between Beverton-Holt parameters",
    subtitle = "Diagnostic for parameter uncertainty and identifiability",
    x = "BH a parameter",
    y = "BH b parameter"
  )

save_plot(
  p_boot_ab,
  "17_bootstrap_BH_a_b_correlation.png"
)

cor(boot_pars$a, boot_pars$b)

# ------------------------------------------------------------
# 8j. Contrast in spawning biomass
# Low contrast can weaken SRR identifiability
# ------------------------------------------------------------

ssb_contrast_table <- sr_data %>%
  summarise(
    min_SSB = min(ssb, na.rm = TRUE),
    q05_SSB = quantile(ssb, 0.05, na.rm = TRUE),
    median_SSB = median(ssb, na.rm = TRUE),
    mean_SSB = mean(ssb, na.rm = TRUE),
    q95_SSB = quantile(ssb, 0.95, na.rm = TRUE),
    max_SSB = max(ssb, na.rm = TRUE),
    contrast_max_min = max_SSB / min_SSB,
    contrast_q95_q05 = q95_SSB / q05_SSB,
    n_years = n(),
    interpretation = case_when(
      contrast_q95_q05 < 2 ~
        "Low SSB contrast; SRR parameters may be weakly identifiable",
      contrast_q95_q05 < 4 ~
        "Moderate SSB contrast",
      TRUE ~
        "High SSB contrast"
    )
  )

write_csv(
  ssb_contrast_table,
  file.path(tab_dir, "29_SSB_contrast_diagnostic.csv")
)

# ------------------------------------------------------------
# 8k. Conditioning decision table
# ------------------------------------------------------------

conditioning_decisions <- tibble(
  question = c(
    "Should a stock-recruitment relationship be used?",
    "Which SR relationship should be used in the base OM?",
    "Should recruitment deviations be lognormal?",
    "Should recruitment deviations be autocorrelated?",
    "Should productivity regimes be included?",
    "Should biomass-dependent recruitment variance be included?",
    "Should environmental covariates be included?",
    "Should SR parameter uncertainty be explored?"
  ),
  diagnostic = c(
    "Biological plausibility and SR model comparison",
    "BH vs segmented comparison, residual diagnostics",
    "Residual distribution, QQ plot, Shapiro-Wilk",
    "ACF, PACF, Ljung-Box, AR(1)",
    "Changepoint and productivity diagnostics",
    "Absolute residuals vs SSB, biomass-class tests",
    "No tested mechanistic environmental driver",
    "Bootstrap BH parameters"
  ),
  base_OM_decision = c(
    "Yes",
    "Beverton-Holt",
    "Yes",
    "No",
    "No",
    "No",
    "No",
    "Not in base OM; use sensitivity scenarios"
  )
)

write_csv(
  conditioning_decisions,
  file.path(tab_dir, "30_conditioning_decision_table.csv")
)
# ------------------------------------------------------------
# 9. Proposed OM scenarios for later MSE runs
# ------------------------------------------------------------
scenario_table <- tibble(
  OM_scenario = c(
    "base",
    "high_sigmaR_Thorson",
    "high_sigmaR_150",
    "low_productivity_075",
    "low_productivity_050",
    "AR1_recruitment_Thorson",
    "poor_recruitment_tail"
  ),
  SR_model = "Beverton-Holt",
  a_multiplier = c(1.00, 1.00, 1.00, 0.75, 0.50, 1.00, 1.00),
  b_multiplier = c(1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
  sigmaR_value = c(
    sigmaR_om,
    0.74,
    sigmaR_om * 1.5,
    sigmaR_om,
    sigmaR_om,
    sigmaR_om,
    sigmaR_om
  ),
  AR1_rho = c(0, 0, 0, 0, 0, 0.43, 0),
  rationale = c(
    "Reference OM based on fitted BH and observed recruitment residuals.",
    "High recruitment variability scenario using Thorson et al. 2014 global mean sigmaR.",
    "Strong increase in recruitment variability relative to ane9aS residuals.",
    "Moderate low-productivity regime.",
    "Severe low-productivity regime.",
    "Exploratory scenario using Thorson et al. 2014 mean autocorrelation.",
    "Stress-test based on observed lower-tail recruitment residuals."
  )
)

write_csv(scenario_table, file.path(tab_dir, "09_proposed_OM_scenarios.csv"))



# ------------------------------------------------------------
# 10. Curves for plotting
# ------------------------------------------------------------

ssb_grid <- seq(
  min(c(sr_obs$ssb, Blim, Bpa), na.rm = TRUE) * 0.5,
  max(sr_obs$ssb, na.rm = TRUE) * 1.5,
  length.out = 400
)

om_curve <- tibble(
  ssb = ssb_grid,
  rec_med = om_bh(ssb_grid, a_om, b_om),
  model = "OM FLSR BH"
)

ss3_curve <- tibble(
  ssb = ssb_grid,
  rec_med = ss3_bh(ssb_grid, R0_ss3, h_ss3, SSB0_ss3),
  model = "SS3 BH"
) %>%
  bind_cols(make_lognormal_bands(.$rec_med, sigmaR_ss3))

seg_curve <- tibble(
  ssb = ssb_grid,
  rec_med = as.numeric(predict(fit_seg, ssb = FLQuant(ssb_grid))),
  model = "OM segmented"
)

p_sr_models <- ggplot() +
  geom_point(data = sr_obs, aes(x = ssb, y = rec), alpha = 0.65, size = 1.8) +
  geom_line(data = om_curve, aes(x = ssb, y = rec_med, colour = "Beverton-Holt"), linewidth = 1.1) +
  geom_line(data = seg_curve, aes(x = ssb, y = rec_med, colour = "Segmented"), linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Alternative stock-recruitment models",
    subtitle = "Beverton-Holt and segmented relationships fitted to historical recruitment",
    x = "SSB",
    y = "Recruitment",
    colour = "SR model"
  )

save_plot(p_sr_models, "11_SR_models_BH_vs_segmented.png")

# ------------------------------------------------------------
# 11. Figures for WD
# ------------------------------------------------------------

p_rec_ts <- ggplot(sr_data, aes(year, rec)) +
  geom_line() +
  geom_point(size = 1.8) +
  theme_bw() +
  labs(
    title = "Historical recruitment",
    subtitle = paste0(first_yr, "-", last_obs_yr, ", recruitment season Q", ss_rec),
    x = "Year",
    y = "Recruitment"
  )

p_sr_compare <- ggplot() +
  geom_ribbon(data = ss3_curve, aes(x = ssb, ymin = rec_lower_95, ymax = rec_upper_95), alpha = 0.12) +
  geom_ribbon(data = ss3_curve, aes(x = ssb, ymin = rec_lower_80, ymax = rec_upper_80), alpha = 0.22) +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(data = sr_obs, aes(x = ssb, y = rec), alpha = 0.65, size = 1.8) +
  geom_line(data = ss3_curve, aes(x = ssb, y = rec_med, colour = "SS3 BH"), linewidth = 1.1) +
  geom_line(data = om_curve, aes(x = ssb, y = rec_med, colour = "OM FLSR BH"), linewidth = 1.1, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Stock-recruitment relationship: OM vs SS3",
    subtitle = paste0("Shaded bands: 80% and 95% recruitment envelopes using ", sigmaR_source, " sigmaR = ", round(sigmaR_ss3, 3), ". Vertical lines: Blim and Bpa."),
    x = "SSB",
    y = "Recruitment",
    colour = "Curve"
  )

p_resid_ts <- ggplot(rec_diag, aes(year, rec_resid_om)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(aes(yintercept = sigmaR_lower_95), linetype = "dotted") +
  geom_hline(aes(yintercept = sigmaR_upper_95), linetype = "dotted") +
  geom_line() +
  geom_point(aes(shape = residual_class), size = 2) +
  theme_bw() +
  labs(
    title = "OM recruitment residuals through time",
    subtitle = "Dotted lines represent approximate 95% limits using sigmaR",
    x = "Year",
    y = "Recruitment residual"
  )

p_resid_hist <- ggplot(rec_diag, aes(rec_resid_om)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12, fill = "grey80", colour = "black") +
  geom_density(linewidth = 1, adjust = 1.3) +
  stat_function(fun = dnorm, args = list(mean = -0.5 * sigmaR_ss3^2, sd = sigmaR_ss3), linewidth = 1, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Distribution of OM recruitment residuals",
    subtitle = "Dashed curve: expected normal distribution under the lognormal recruitment assumption",
    x = "Recruitment residual",
    y = "Density"
  )

p_qq <- ggplot(rec_diag, aes(sample = rec_resid_om)) +
  stat_qq() +
  stat_qq_line() +
  theme_bw() +
  labs(
    title = "QQ plot of OM recruitment residuals",
    x = "Theoretical quantiles",
    y = "Observed residual quantiles"
  )

p_resid_ssb <- ggplot(rec_diag, aes(ssb, rec_resid_om)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(aes(yintercept = sigmaR_lower_95), linetype = "dotted") +
  geom_hline(aes(yintercept = sigmaR_upper_95), linetype = "dotted") +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(aes(shape = residual_class), size = 2) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "OM recruitment residuals vs SSB",
    subtitle = "Used as a diagnostic for residual structure and potential low-SSB bias",
    x = "SSB",
    y = "Recruitment residual"
  )

p_tail <- ggplot(tail_diagnostics, aes(normal_expected, observed_residual)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_text(aes(label = probability), nudge_y = 0.05, size = 3) +
  theme_bw() +
  labs(
    title = "Tail diagnostics",
    subtitle = "Observed recruitment residual quantiles compared with expected normal quantiles",
    x = "Expected residual quantile",
    y = "Observed residual quantile"
  )


save_plot(p_rec_ts, "01_historical_recruitment_timeseries.png")
save_plot(p_sr_compare, "02_SR_curve_OM_vs_SS3_sigmaR.png")
save_plot(p_resid_ts, "03_OM_residuals_time_series_sigmaR.png")
save_plot(p_resid_hist, "04_OM_residual_distribution_vs_sigmaR.png")
save_plot(p_qq, "05_OM_residual_QQplot.png")
save_plot(p_resid_ssb, "06_OM_residuals_vs_SSB_sigmaR.png")
save_plot(p_tail, "07_tail_diagnostics.png")

png(file.path(fig_dir, "08_OM_residual_ACF.png"), width = 900, height = 700)
acf(rec_diag$rec_resid_om, main = "ACF of OM recruitment residuals")
dev.off()

png(file.path(fig_dir, "09_OM_residual_PACF.png"), width = 900, height = 700)
pacf(rec_diag$rec_resid_om, main = "PACF of OM recruitment residuals")
dev.off()


# ------------------------------------------------------------
# Residual magnitude vs SSB
# Diagnostic for stock-size-dependent recruitment variability
# ------------------------------------------------------------

rec_diag <- rec_diag %>%
  mutate(
    abs_rec_resid_om = abs(rec_resid_om),
    biomass_class = ifelse(
      ssb < median(ssb, na.rm = TRUE),
      "Low SSB",
      "High SSB"
    )
  )

p_abs_resid_ssb <- ggplot(rec_diag, aes(ssb, abs_rec_resid_om)) +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Absolute recruitment residuals vs SSB",
    subtitle = "Diagnostic for stock-size-dependent recruitment variability",
    x = "SSB",
    y = "|Recruitment residual|"
  )

save_plot(
  p_abs_resid_ssb,
  "06b_absolute_OM_residuals_vs_SSB.png"
)

# Spearman correlation between residual magnitude and SSB
abs_resid_spearman <- cor.test(
  rec_diag$abs_rec_resid_om,
  rec_diag$ssb,
  method = "spearman",
  exact = FALSE
)

abs_resid_spearman_table <- tibble(
  test = "Spearman correlation",
  variable_x = "|recruitment residual|",
  variable_y = "SSB",
  statistic = unname(abs_resid_spearman$statistic),
  estimate = unname(abs_resid_spearman$estimate),
  p_value = abs_resid_spearman$p.value
)

write_csv(
  abs_resid_spearman_table,
  file.path(tab_dir, "24_abs_residuals_vs_SSB_spearman.csv")
)

abs_resid_by_biomass <- rec_diag %>%
  group_by(biomass_class) %>%
  summarise(
    n = n(),
    mean_abs_resid = mean(abs_rec_resid_om, na.rm = TRUE),
    median_abs_resid = median(abs_rec_resid_om, na.rm = TRUE),
    sd_abs_resid = sd(abs_rec_resid_om, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  abs_resid_by_biomass,
  file.path(tab_dir, "25_abs_residuals_by_biomass_class.csv")
)

abs_resid_wilcox <- wilcox.test(
  abs_rec_resid_om ~ biomass_class,
  data = rec_diag
)

abs_resid_wilcox_table <- tibble(
  test = "Wilcoxon rank-sum test",
  comparison = "Absolute recruitment residuals by biomass class",
  statistic = unname(abs_resid_wilcox$statistic),
  p_value = abs_resid_wilcox$p.value
)

write_csv(
  abs_resid_wilcox_table,
  file.path(tab_dir, "26_abs_residuals_by_biomass_class_wilcox.csv")
)

p_abs_resid_box <- ggplot(
  rec_diag,
  aes(biomass_class, abs_rec_resid_om)
) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Absolute recruitment residuals by biomass class",
    x = "Biomass class",
    y = "|Recruitment residual|"
  )

save_plot(
  p_abs_resid_box,
  "06c_absolute_residuals_by_biomass_class.png"
)
#############################################################################

library(zoo)

sr_data <- sr_data %>%
  mutate(
    rec_roll5 = zoo::rollmean(
      rec,
      k = 5,
      fill = NA,
      align = "center"
    )
  )

p_roll <- ggplot(sr_data, aes(year)) +
  geom_line(aes(y = rec), alpha = 0.35) +
  geom_point(aes(y = rec), alpha = 0.5, size = 1.5) +
  geom_line(aes(y = rec_roll5), linewidth = 1.2) +
  theme_bw() +
  labs(
    title = "Historical recruitment with 5-year moving average",
    subtitle = "Visual diagnostic for persistent productivity changes",
    x = "Year",
    y = "Recruitment"
  )

save_plot(p_roll, "10_recruitment_rolling_mean_5yr.png")

library(randtests)

sign_series <- ifelse(rec_diag$rec_resid_om > 0, 1, 0)

runs_out <- runs.test(sign_series)

runs_table <- tibble(
  test = "Runs test",
  statistic = unname(runs_out$statistic),
  p_value = runs_out$p.value
)

write_csv(
  runs_table,
  file.path(tab_dir, "16_runs_test.csv")
)
# ------------------------------------------------------------
# 12. Save RData object
# ------------------------------------------------------------
save(
  fit_bh_om,
  fit_seg,
  mod_bh_om,
  mod_seg,
  sr_data,
  sr_obs,
  sr_parameter_table,
  sr_model_comparison,
  compare_summary,
  rec_diag,
  resid_summary,
  resid_tests,
  runs_table,
  tail_diagnostics,
  extreme_residual_years,
  scenario_table,
  om_curve,
  ss3_curve,
  seg_curve,
  prod_cp_table,
  period_summary,
  productivity_tests,
  boot_pars,
  boot_summary,
  cp_prod,
  R0_ss3,
  h_ss3,
  SSB0_ss3,
  sigmaR_ss3,
  sigmaR_om,
  sigmaR_source,
  a_om,
  b_om,
  file = file.path(data_dir, "recruitment_conditioning_OM_vs_SS3.RData")
)

message("Recruitment conditioning diagnostics completed.")
message("Tables saved in: ", tab_dir)
message("Figures saved in: ", fig_dir)
message("RData saved in: ", file.path(data_dir, "recruitment_conditioning_OM_vs_SS3.RData"))


sr_data <- sr_data %>%
  mutate(
    regime =
      ifelse(
        productivity >
          median(productivity),
        "High",
        "Low"
      )
  )

rle_regime <- rle(sr_data$regime)

tibble(
  regime = rle_regime$values,
  duration = rle_regime$lengths
)

tibble(
  regime = rle_regime$values,
  duration = rle_regime$lengths
) %>%
  group_by(regime) %>%
  summarise(
    n_runs = n(),
    mean_duration = mean(duration),
    max_duration = max(duration)
  )


# comparar BH vs Ricker vs Segmented

mod_ricker <- FLSR(
  rec = rec_fit,
  ssb = ssb_fit,
  model = ricker
)

fit_ricker <- fmle(mod_ricker)

sr_model_comparison <- tibble(
  model = c("Beverton-Holt", "Ricker", "Segmented"),
  logLik = c(
    as.numeric(logLik(fit_bh_om)),
    as.numeric(logLik(fit_ricker)),
    as.numeric(logLik(fit_seg))
  ),
  AIC = c(
    AIC(fit_bh_om),
    AIC(fit_ricker),
    AIC(fit_seg)
  )
) %>%
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE))

write_csv(
  sr_model_comparison,
  file.path(tab_dir, "31_SR_model_comparison_BH_Ricker_segmented.csv")
)

# Evaluarsi hay caída del reclutamiento total a SSB altas
sr_high_low <- sr_data %>%
  mutate(
    ssb_class = ifelse(
      ssb >= quantile(ssb, 0.75, na.rm = TRUE),
      "High SSB",
      "Other"
    )
  ) %>%
  group_by(ssb_class) %>%
  summarise(
    n = n(),
    mean_rec = mean(rec, na.rm = TRUE),
    median_rec = median(rec, na.rm = TRUE),
    geomean_rec = exp(mean(log(rec), na.rm = TRUE)),
    mean_RPS = mean(RPS, na.rm = TRUE),
    median_RPS = median(RPS, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  sr_high_low,
  file.path(tab_dir, "32_recruitment_at_high_SSB.csv")
)

ricker_curve <- tibble(
  ssb = ssb_grid,
  rec_med = as.numeric(predict(fit_ricker, ssb = FLQuant(ssb_grid))),
  model = "Ricker"
)

p_bh_ricker <- ggplot() +
  geom_point(data = sr_data, aes(ssb, rec), size = 2, alpha = 0.7) +
  geom_line(data = om_curve, aes(ssb, rec_med, colour = "Beverton-Holt"), linewidth = 1.1) +
  geom_line(data = ricker_curve, aes(ssb, rec_med, colour = "Ricker"), linewidth = 1.1) +
  geom_line(data = seg_curve, aes(ssb, rec_med, colour = "Segmented"), linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Compensation vs overcompensation",
    subtitle = "Comparison of Beverton-Holt, Ricker and segmented SR relationships",
    x = "SSB",
    y = "Recruitment",
    colour = "Model"
  )

save_plot(
  p_bh_ricker,
  "18_compensation_vs_overcompensation_BH_Ricker_segmented.png"
)


# ------------------------------------------------------------
# 8l. Quarterly spawner age-structure diagnostics
# SSB in Q2 -> recruitment in Q3, same year
# ------------------------------------------------------------
# ============================================================
# Helper: extract numeric FLQuant from FLBiol slot
# ============================================================

get_flq <- function(x, slot_name) {
  
  obj <- slot(x, slot_name)
  
  # If it is already an FLQuant
  if (inherits(obj, "FLQuant")) {
    return(obj)
  }
  
  # If it is an FLQuants or list-like object
  if (inherits(obj, "FLQuants") || is.list(obj)) {
    return(obj[[1]])
  }
  
  # If it has a slot called .Data
  if (isS4(obj) && ".Data" %in% slotNames(obj)) {
    return(obj@.Data[[1]])
  }
  
  stop("Could not extract FLQuant from slot: ", slot_name)
}

# ============================================================
# Quarterly spawner age structure
# SSB_Q2(year t) -> Recruitment_Q3(year t)
# ============================================================

n_all   <- get_flq(ane_biol, "n")
wt_all  <- get_flq(ane_biol, "wt")
mat_all <- get_flq(ane_biol, "mat")

n_q2 <- n_all[, ac(first_yr:last_obs_yr), , ss_ssb, drop = FALSE]
wt_q2 <- wt_all[, ac(first_yr:last_obs_yr), , ss_ssb, drop = FALSE]
mat_q2 <- mat_all[, ac(first_yr:last_obs_yr), , ss_ssb, drop = FALSE]

ssb_age_q2 <- n_q2 * wt_q2 * mat_q2

ssb_age_df <- as.data.frame(ssb_age_q2) %>%
  as_tibble() %>%
  rename(ssb_age_q2 = data) %>%
  mutate(
    age = as.numeric(as.character(age)),
    year = as.numeric(as.character(year))
  )

spawner_structure_q2 <- ssb_age_df %>%
  group_by(year) %>%
  mutate(
    ssb_q2_total = sum(ssb_age_q2, na.rm = TRUE),
    prop_ssb_q2 = ssb_age_q2 / ssb_q2_total
  ) %>%
  summarise(
    mean_spawner_age_q2 = sum(age * prop_ssb_q2, na.rm = TRUE),
    prop_ssb_q2_age1 = sum(prop_ssb_q2[age == 1], na.rm = TRUE),
    prop_ssb_q2_age2plus = sum(prop_ssb_q2[age >= 2], na.rm = TRUE),
    ssb_q2_total = first(ssb_q2_total),
    .groups = "drop"
  ) %>%
  left_join(
    rec_diag %>%
      select(year, rec_obs, rec_pred_om, rec_resid_om, rec_mult_om, ssb),
    by = "year"
  )

write_csv(
  spawner_structure_q2,
  file.path(tab_dir, "33_spawner_age_structure_Q2_by_year.csv")
)

spawner_structure_q2
# ============================================================
# Maternal / spawner age-structure diagnostics
# Q2 spawner structure -> Q3 recruitment residuals
# ============================================================

maternal_tests_q2 <- bind_rows(
  
  broom::tidy(
    lm(rec_resid_om ~ mean_spawner_age_q2,
       data = spawner_structure_q2)
  ) %>%
    mutate(model = "Q3 recruitment residual ~ Q2 mean spawner age"),
  
  broom::tidy(
    lm(rec_resid_om ~ prop_ssb_q2_age1,
       data = spawner_structure_q2)
  ) %>%
    mutate(model = "Q3 recruitment residual ~ Q2 proportion SSB age 1"),
  
  broom::tidy(
    lm(rec_resid_om ~ prop_ssb_q2_age2plus,
       data = spawner_structure_q2)
  ) %>%
    mutate(model = "Q3 recruitment residual ~ Q2 proportion SSB age 2+")
)

write_csv(
  maternal_tests_q2,
  file.path(tab_dir, "34_maternal_effects_Q2_spawners_Q3_recruitment.csv")
)

# Figure 1: mean spawner age
p_mean_age_resid_q2 <- ggplot(
  spawner_structure_q2,
  aes(mean_spawner_age_q2, rec_resid_om)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Q3 recruitment residuals vs Q2 mean spawner age",
    subtitle = "Quarterly diagnostic: SSB_Q2(year t) -> Recruitment_Q3(year t)",
    x = "Mean spawner age in Q2",
    y = "Q3 recruitment residual"
  )

save_plot(
  p_mean_age_resid_q2,
  "19_Q3_recruitment_residuals_vs_Q2_mean_spawner_age.png"
)

# Figure 2: proportion of SSB from age 1
p_prop_age1_resid_q2 <- ggplot(
  spawner_structure_q2,
  aes(prop_ssb_q2_age1, rec_resid_om)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Q3 recruitment residuals vs Q2 proportion of SSB from age 1",
    subtitle = "Quarterly diagnostic: young-spawner contribution in Q2",
    x = "Proportion of Q2 SSB from age 1",
    y = "Q3 recruitment residual"
  )

save_plot(
  p_prop_age1_resid_q2,
  "20_Q3_recruitment_residuals_vs_Q2_prop_SSB_age1.png"
)

# Figure 3: proportion of SSB from age 2+
p_prop_age2plus_resid_q2 <- ggplot(
  spawner_structure_q2,
  aes(prop_ssb_q2_age2plus, rec_resid_om)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Q3 recruitment residuals vs Q2 proportion of SSB from age 2+",
    subtitle = "Quarterly diagnostic: older-spawner contribution in Q2",
    x = "Proportion of Q2 SSB from age 2+",
    y = "Q3 recruitment residual"
  )

save_plot(
  p_prop_age2plus_resid_q2,
  "21_Q3_recruitment_residuals_vs_Q2_prop_SSB_age2plus.png"
)
# ------------------------------------------------------------
# 8n. Quarterly cohort strength diagnostics
# Cohort strength is defined using recruitment in Q3
# ------------------------------------------------------------

cohort_strength_q3 <- rec_diag %>%
  mutate(
    cohort = year,
    recruitment_season = ss_rec,
    spawning_season = ss_ssb,
    log_rec_q3 = log(rec_obs),
    log_rec_q3_anom = as.numeric(scale(log_rec_q3)),
    cohort_class = case_when(
      log_rec_q3_anom >= 1.0  ~ "Strong Q3 cohort",
      log_rec_q3_anom <= -1.0 ~ "Weak Q3 cohort",
      TRUE ~ "Average Q3 cohort"
    )
  ) %>%
  select(
    cohort,
    year,
    spawning_season,
    recruitment_season,
    ssb_q2 = ssb,
    rec_q3 = rec_obs,
    rec_pred_om_q3 = rec_pred_om,
    rec_resid_om_q3 = rec_resid_om,
    rec_mult_om_q3 = rec_mult_om,
    log_rec_q3_anom,
    cohort_class
  )

write_csv(
  cohort_strength_q3,
  file.path(tab_dir, "35_Q3_cohort_strength_classification.csv")
)

cohort_summary_q3 <- cohort_strength_q3 %>%
  group_by(cohort_class) %>%
  summarise(
    n = n(),
    mean_rec_q3 = mean(rec_q3, na.rm = TRUE),
    geomean_rec_q3 = exp(mean(log(rec_q3), na.rm = TRUE)),
    mean_resid_q3 = mean(rec_resid_om_q3, na.rm = TRUE),
    median_resid_q3 = median(rec_resid_om_q3, na.rm = TRUE),
    mean_ssb_q2 = mean(ssb_q2, na.rm = TRUE),
    median_ssb_q2 = median(ssb_q2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cohort_summary_q3,
  file.path(tab_dir, "36_Q3_cohort_strength_summary.csv")
)

p_cohort_strength_q3 <- ggplot(
  cohort_strength_q3,
  aes(year, rec_q3)
) +
  geom_line(alpha = 0.4) +
  geom_point(aes(shape = cohort_class), size = 2.3) +
  theme_bw() +
  labs(
    title = "Historical recruitment cohorts in Q3",
    subtitle = "Strong and weak cohorts based on standardized log recruitment in Q3",
    x = "Cohort year",
    y = "Recruitment in Q3",
    shape = "Cohort class"
  )

save_plot(
  p_cohort_strength_q3,
  "21_Q3_cohort_strength_classification.png"
)

p_cohort_resid_q3 <- ggplot(
  cohort_strength_q3,
  aes(year, rec_resid_om_q3)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = cohort_class), size = 2.3) +
  geom_line(alpha = 0.4) +
  theme_bw() +
  labs(
    title = "Q3 recruitment residuals by cohort",
    subtitle = "Strong and weak Q3 cohorts relative to the fitted BH relationship",
    x = "Cohort year",
    y = "Q3 recruitment residual",
    shape = "Cohort class"
  )

save_plot(
  p_cohort_resid_q3,
  "22_Q3_recruitment_residuals_by_cohort_strength.png"
)

# ------------------------------------------------------------
# 8o. Are strong/weak Q3 cohorts associated with Q2 SSB?
# ------------------------------------------------------------

cohort_ssb_q2_summary <- cohort_strength_q3 %>%
  group_by(cohort_class) %>%
  summarise(
    n = n(),
    mean_ssb_q2 = mean(ssb_q2, na.rm = TRUE),
    median_ssb_q2 = median(ssb_q2, na.rm = TRUE),
    min_ssb_q2 = min(ssb_q2, na.rm = TRUE),
    max_ssb_q2 = max(ssb_q2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cohort_ssb_q2_summary,
  file.path(tab_dir, "37_Q3_cohort_strength_vs_Q2_SSB_summary.csv")
)

cohort_ssb_q2_lm <- lm(log_rec_q3_anom ~ ssb_q2, data = cohort_strength_q3)

cohort_ssb_q2_lm_table <- broom::tidy(cohort_ssb_q2_lm) %>%
  mutate(
    model = "Standardized log recruitment Q3 ~ SSB Q2",
    interpretation = case_when(
      term == "ssb_q2" & p.value < 0.05 & estimate > 0 ~
        "Positive association between Q2 SSB and Q3 cohort strength",
      term == "ssb_q2" & p.value < 0.05 & estimate < 0 ~
        "Negative association between Q2 SSB and Q3 cohort strength",
      term == "ssb_q2" & p.value >= 0.05 ~
        "No clear linear association between Q2 SSB and Q3 cohort strength",
      TRUE ~ NA_character_
    )
  )

write_csv(
  cohort_ssb_q2_lm_table,
  file.path(tab_dir, "38_Q3_cohort_strength_vs_Q2_SSB_lm.csv")
)

p_cohort_ssb_q2 <- ggplot(
  cohort_strength_q3,
  aes(ssb_q2, log_rec_q3_anom)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(aes(shape = cohort_class), size = 2.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_bw() +
  labs(
    title = "Q3 cohort strength vs Q2 SSB",
    subtitle = "Quarterly diagnostic: SSB_Q2(year t) -> Recruitment_Q3(year t)",
    x = "SSB in Q2",
    y = "Standardized log recruitment in Q3",
    shape = "Cohort class"
  )

save_plot(
  p_cohort_ssb_q2,
  "23_Q3_cohort_strength_vs_Q2_SSB.png"
)

# ------------------------------------------------------------
# 8p. Updated quarterly conditioning decision table
# ------------------------------------------------------------

conditioning_decisions_quarterly <- tibble(
  question = c(
    "Should a stock-recruitment relationship be used?",
    "Which SR relationship should be used in the base OM?",
    "Should Ricker overcompensation be used?",
    "Should recruitment deviations be lognormal?",
    "Should recruitment deviations be autocorrelated?",
    "Should productivity regimes be included?",
    "Should biomass-dependent recruitment variance be included?",
    "Should environmental covariates be included?",
    "Should maternal or spawner-age effects be included?",
    "Should sporadic recruitment events be included?",
    "Should SR parameter uncertainty be explored?"
  ),
  quarterly_diagnostic = c(
    "Q2 SSB linked to Q3 recruitment",
    "BH vs segmented vs Ricker comparison using Q2 SSB and Q3 recruitment",
    "Ricker comparison and recruitment at high Q2 SSB",
    "Q3 recruitment residual distribution, QQ plot, Shapiro-Wilk and tail diagnostics",
    "ACF, PACF, Ljung-Box, AR(1), runs test on Q3 recruitment residuals",
    "Changepoint, breakpoints and productivity diagnostics using Q3 recruitment",
    "Absolute Q3 recruitment residuals vs Q2 SSB, biomass-class tests",
    "No tested mechanistic environmental driver",
    "Q3 recruitment residuals vs Q2 spawner age structure",
    "Extreme Q3 residual years and Q3 cohort strength classification",
    "Bootstrap BH parameters using Q2 SSB and Q3 recruitment"
  ),
  base_OM_decision = c(
    "Yes",
    "Beverton-Holt",
    "No",
    "Yes",
    "No",
    "No",
    "No",
    "No",
    "No, unless diagnostics show a clear relationship",
    "No as explicit process; retained through stochastic recruitment variability",
    "Not in base OM; use sensitivity scenarios"
  )
)

write_csv(
  conditioning_decisions_quarterly,
  file.path(tab_dir, "39_conditioning_decision_table_quarterly.csv")
)
