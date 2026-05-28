# ============================================================
# Recruitment conditioning diagnostics for OM
# Anchovy 9a South
# OM vs SS3 SR relationship + residual diagnostics
# ============================================================

rm(list = ls())

library(FLBEIA)
library(ggplotFL)
library(r4ss)
library(icesTAF)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(here)

set.seed(123)

# ============================================================
# 1. Settings
# ============================================================

wd <- here()
setwd(wd)

first.yr <- 1989
last.obs.yr <- 2024
proj.yr <- 2025
proj.nyr <- 10

hist.yrs <- first.yr:last.obs.yr
proj.yrs <- proj.yr:(proj.yr + proj.nyr - 1)
last.yr <- max(proj.yrs)

ni <- 50
ns <- 4

ss.rec <- 3
ss.ssb <- 2

Blim <- 4721
Bpa  <- 6561

stk_ane9aS_rds <- "boot/data/stk_ane9aS.rds"
ss3_rds        <- "boot/data/ss3_ane9aS.rds"

out.dir  <- here("outputs/recruitment_diagnostics")
data.dir <- here("data/processed")

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Load historical stock object
# ============================================================

ane.stock <- readRDS(file.path(stk_ane9aS_rds))
ane.stock <- propagate(ane.stock, ni)

name(ane.stock) <- "ANE"
desc(ane.stock) <- "WGHANSA2025 assessment output"
ane.stock@harvest@units <- "f"

# ============================================================
# 3. Biological timing
# ============================================================

fspwn <- harvest.spwn(ane.stock)
fspwn[] <- 0
fspwn[,,, "1", ] <- 1
harvest.spwn(ane.stock) <- fspwn

mspwn <- m.spwn(ane.stock)
mspwn[] <- 0
mspwn[,,, "1", ] <- 1
m.spwn(ane.stock) <- mspwn

mat(ane.stock)[ac(0), ] <- 0
mat(ane.stock)[,,, 1] <- 0
mat(ane.stock)[,,, 3] <- 0
mat(ane.stock)[,,, 4] <- 0

range(ane.stock, "minfbar") <- 3
range(ane.stock, "maxfbar") <- 3

# ============================================================
# 4. Create FLBiol
# ============================================================

ane <- FLBiol(
  n = stock.n(ane.stock),
  wt = stock.wt(ane.stock),
  m = m(ane.stock),
  spwn = m.spwn(ane.stock),
  mat = mat(ane.stock),
  fec = predictModel(
    FLQuants(fec = ane.stock@mat * 0 + 1),
    model = ~fec
  ),
  name = "ANE",
  desc = ane.stock@desc,
  range = ane.stock@range
)

units(ane)$m <- units(fec(ane)) <- units(mat(ane)) <- ""

ane <- window(ane, start = first.yr, end = last.yr)

biols <- FLBiols(ANE = ane)

# ============================================================
# 5. Fill projection biological inputs
# ============================================================

mean.yrs <- ac((last.obs.yr - 2):last.obs.yr)

m(biols$ANE)[, ac(proj.yrs), ]    <- yearMeans(m(biols$ANE[, mean.yrs]))
wt(biols$ANE)[, ac(proj.yrs), ]   <- yearMeans(wt(biols$ANE[, mean.yrs]))
mat(biols$ANE)[, ac(proj.yrs), ]  <- yearMeans(mat(biols$ANE[, mean.yrs]))
fec(biols$ANE)[, ac(proj.yrs), ]  <- yearMeans(fec(biols$ANE[, mean.yrs]))
spwn(biols$ANE)[, ac(proj.yrs), ] <- yearMeans(spwn(biols$ANE[, mean.yrs]))

# ============================================================
# 6. Fit Beverton-Holt stock-recruitment model for OM
# ============================================================

rec_h <- rec(ane)[,,, ss.rec, drop = FALSE]
ssb_h <- ssb(ane)[,,, ss.ssb, drop = FALSE]

# Use one iteration only for SR fitting
rec_fit <- iter(rec_h, 1)
ssb_fit <- iter(ssb_h, 1)

rec_fit <- window(rec_fit, start = first.yr, end = last.obs.yr)
ssb_fit <- window(ssb_fit, start = first.yr, end = last.obs.yr)

sr_check <- tibble(
  year = as.numeric(dimnames(rec_fit)$year),
  ssb = as.numeric(ssb_fit),
  rec = as.numeric(rec_fit)
)

stopifnot(
  nrow(sr_check) == length(hist.yrs),
  all(is.finite(sr_check$ssb)),
  all(is.finite(sr_check$rec)),
  all(sr_check$ssb > 0),
  all(sr_check$rec > 0)
)

mod2_bh <- FLSR(
  rec = rec_fit,
  ssb = ssb_fit,
  model = bevholt
)

fit2_bh <- fmle(mod2_bh)

hist_resid <- as.numeric(residuals(fit2_bh))
hist_resid <- hist_resid[is.finite(hist_resid)]

sigmaR_om <- sqrt(var(hist_resid, na.rm = TRUE))

om_pars <- as.numeric(params(fit2_bh))
names(om_pars) <- dimnames(params(fit2_bh))$params

a_om <- om_pars["a"]
b_om <- om_pars["b"]

if (!is.finite(a_om) | !is.finite(b_om)) {
  stop("OM parameters a or b are not finite. Check params(fit2_bh).")
}

# ============================================================
# 7. Load SS3 and extract SR parameters
# ============================================================

ss3rep <- readRDS(here(ss3_rds))

pars <- ss3rep$parameters
dq   <- ss3rep$derived_quants

find_value <- function(df, patterns, label_col = "Label", value_col = "Value") {
  pattern <- paste(patterns, collapse = "|")
  
  out <- df %>%
    filter(str_detect(.data[[label_col]], regex(pattern, ignore_case = TRUE)))
  
  if (nrow(out) == 0) return(NA_real_)
  
  out[[value_col]][1]
}

lnR0_ss3 <- find_value(
  pars,
  c("SR_LN\\(R0\\)", "SR_Ln\\(R0\\)", "SR_log\\(R0\\)", "SR_LN_R0", "SR_lnR0")
)

R0_ss3 <- exp(lnR0_ss3)

h_ss3 <- find_value( pars, c("SR_BH_steep", "SR_steep", "steep"))

sigmaR_ss3 <- find_value( pars, c("SR_sigmaR", "SR_sigma", "sigmaR", "SigmaR", "rec_sigma", "ln\\(R\\)_sigma"))

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
  parameter = c(
    "a_OM", "b_OM",
    "R0_SS3", "h_SS3", "SSB0_SS3",
    "sigmaR_used", "sigmaR_OM_resid"
  ),
  value = c(
    a_om, b_om,
    R0_ss3, h_ss3, SSB0_ss3,
    sigmaR_ss3, sigmaR_om
  ),
  source = c(
    "FLSR", "FLSR",
    "SS3", "SS3", "SS3",
    sigmaR_source, "OM residuals"
  )
)

print(sr_parameter_table)

# ============================================================
# 8. Define SR functions and curves
# ============================================================

ss3_bh <- function(SSB, R0, h, SSB0) {
  (4 * h * R0 * SSB) /
    (SSB0 * (1 - h) + SSB * (5 * h - 1))
}

om_bh <- function(SSB, a, b) {
  a * SSB / (b + SSB)
}

make_lognormal_bands <- function(rec_med, sigmaR) {
  tibble(
    rec_lower_95 = rec_med * exp(-0.5 * sigmaR^2 - 2.00 * sigmaR),
    rec_upper_95 = rec_med * exp(-0.5 * sigmaR^2 + 2.00 * sigmaR),
    rec_lower_80 = rec_med * exp(-0.5 * sigmaR^2 - 1.28 * sigmaR),
    rec_upper_80 = rec_med * exp(-0.5 * sigmaR^2 + 1.28 * sigmaR)
  )
}

sr_obs <- sr_check %>%
  mutate(
    rec_pred_om = om_bh(ssb, a_om, b_om),
    rec_pred_ss3 = ss3_bh(ssb, R0_ss3, h_ss3, SSB0_ss3),
    resid_om = log(rec) - log(rec_pred_om),
    resid_ss3 = log(rec) - log(rec_pred_ss3)
  )

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

# ============================================================
# 9. OM residual diagnostics requested by reviewer
# ============================================================

rec_diag_om <- sr_obs %>%
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
    sigmaR_lower_95 = -0.5 * sigmaR_ss3^2 - 2 * sigmaR_ss3,
    sigmaR_upper_95 = -0.5 * sigmaR_ss3^2 + 2 * sigmaR_ss3,
    below_sigmaR_95 = rec_resid_om < sigmaR_lower_95,
    above_sigmaR_95 = rec_resid_om > sigmaR_upper_95,
    inside_sigmaR_95 = !below_sigmaR_95 & !above_sigmaR_95
  )

om_resid_summary <- rec_diag_om %>%
  summarise(
    n_years = n(),
    mean_resid = mean(rec_resid_om, na.rm = TRUE),
    sd_resid = sd(rec_resid_om, na.rm = TRUE),
    sigmaR_used = sigmaR_ss3,
    sigmaR_source = sigmaR_source,
    q05 = quantile(rec_resid_om, 0.05, na.rm = TRUE),
    q10 = quantile(rec_resid_om, 0.10, na.rm = TRUE),
    median = median(rec_resid_om, na.rm = TRUE),
    q90 = quantile(rec_resid_om, 0.90, na.rm = TRUE),
    q95 = quantile(rec_resid_om, 0.95, na.rm = TRUE),
    min_resid = min(rec_resid_om, na.rm = TRUE),
    max_resid = max(rec_resid_om, na.rm = TRUE),
    prop_inside_sigmaR_95 = mean(inside_sigmaR_95, na.rm = TRUE),
    prop_below_sigmaR_95 = mean(below_sigmaR_95, na.rm = TRUE),
    prop_above_sigmaR_95 = mean(above_sigmaR_95, na.rm = TRUE)
  )

om_resid_tests <- tibble(
  test = c(
    "Ljung-Box lag 1",
    "Ljung-Box lag 3",
    "Ljung-Box lag 5",
    "Shapiro-Wilk normality"
  ),
  statistic = c(
    unname(Box.test(rec_diag_om$rec_resid_om, lag = 1, type = "Ljung-Box")$statistic),
    unname(Box.test(rec_diag_om$rec_resid_om, lag = 3, type = "Ljung-Box")$statistic),
    unname(Box.test(rec_diag_om$rec_resid_om, lag = 5, type = "Ljung-Box")$statistic),
    unname(shapiro.test(rec_diag_om$rec_resid_om)$statistic)
  ),
  p_value = c(
    Box.test(rec_diag_om$rec_resid_om, lag = 1, type = "Ljung-Box")$p.value,
    Box.test(rec_diag_om$rec_resid_om, lag = 3, type = "Ljung-Box")$p.value,
    Box.test(rec_diag_om$rec_resid_om, lag = 5, type = "Ljung-Box")$p.value,
    shapiro.test(rec_diag_om$rec_resid_om)$p.value
  )
)

extreme_residual_years <- rec_diag_om %>%
  arrange(rec_resid_om) %>%
  select(
    year, ssb, rec_obs, rec_pred_om,
    rec_resid_om, rec_mult_om,
    resid_percentile, residual_class,
    inside_sigmaR_95
  )

print(om_resid_summary)
print(om_resid_tests)
print(extreme_residual_years)

# ============================================================
# 10. Quantitative OM vs SS3 comparison
# ============================================================

compare_at_obs <- sr_obs %>%
  mutate(
    ratio_om_ss3 = rec_pred_om / rec_pred_ss3,
    diff_percent = 100 * (rec_pred_om - rec_pred_ss3) / rec_pred_ss3,
    inside_95_sigmaR = rec >= rec_pred_ss3 * exp(-0.5 * sigmaR_ss3^2 - 2 * sigmaR_ss3) &
      rec <= rec_pred_ss3 * exp(-0.5 * sigmaR_ss3^2 + 2 * sigmaR_ss3),
    inside_80_sigmaR = rec >= rec_pred_ss3 * exp(-0.5 * sigmaR_ss3^2 - 1.28 * sigmaR_ss3) &
      rec <= rec_pred_ss3 * exp(-0.5 * sigmaR_ss3^2 + 1.28 * sigmaR_ss3)
  )

summary_compare <- compare_at_obs %>%
  summarise(
    mean_ratio = mean(ratio_om_ss3, na.rm = TRUE),
    min_ratio = min(ratio_om_ss3, na.rm = TRUE),
    max_ratio = max(ratio_om_ss3, na.rm = TRUE),
    mean_diff_percent = mean(diff_percent, na.rm = TRUE),
    mean_abs_percent_diff = mean(abs(diff_percent), na.rm = TRUE),
    prop_inside_80_sigmaR = mean(inside_80_sigmaR, na.rm = TRUE),
    prop_inside_95_sigmaR = mean(inside_95_sigmaR, na.rm = TRUE)
  )

print(summary_compare)

# ============================================================
# 11. Plots
# ============================================================

p_sr_compare_sigmaR <- ggplot() +
  geom_ribbon(
    data = ss3_curve,
    aes(x = ssb, ymin = rec_lower_95, ymax = rec_upper_95),
    alpha = 0.12
  ) +
  geom_ribbon(
    data = ss3_curve,
    aes(x = ssb, ymin = rec_lower_80, ymax = rec_upper_80),
    alpha = 0.22
  ) +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(
    data = sr_obs,
    aes(x = ssb, y = rec),
    alpha = 0.65,
    size = 1.8
  ) +
  geom_line(
    data = ss3_curve,
    aes(x = ssb, y = rec_med, colour = "SS3 BH"),
    linewidth = 1.2
  ) +
  geom_line(
    data = om_curve,
    aes(x = ssb, y = rec_med, colour = "OM FLSR BH"),
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    title = "Stock-recruitment relationship: OM vs SS3",
    subtitle = paste0(
      "Shaded bands: sigmaR envelope (", sigmaR_source,
      "); sigmaR = ", round(sigmaR_ss3, 3),
      ". Vertical lines: Blim and Bpa."
    ),
    x = "SSB",
    y = "Recruitment",
    colour = "Curve"
  )

p_om_resid_ts <- ggplot(rec_diag_om, aes(year, rec_resid_om)) +
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

p_om_resid_hist <- ggplot(rec_diag_om, aes(rec_resid_om)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 12,
    fill = "grey80",
    colour = "black"
  ) +
  geom_density(linewidth = 1, adjust = 1.3) +
  stat_function(
    fun = dnorm,
    args = list(
      mean = -0.5 * sigmaR_ss3^2,
      sd = sigmaR_ss3
    ),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    title = "Distribution of OM recruitment residuals",
    subtitle = "Dashed curve: normal distribution using sigmaR",
    x = "Recruitment residual",
    y = "Density"
  )

p_om_qq <- ggplot(rec_diag_om, aes(sample = rec_resid_om)) +
  stat_qq() +
  stat_qq_line() +
  theme_bw() +
  labs(
    title = "QQ plot of OM recruitment residuals",
    x = "Theoretical quantiles",
    y = "Observed residual quantiles"
  )

p_om_resid_ssb <- ggplot(rec_diag_om, aes(ssb, rec_resid_om)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(aes(yintercept = sigmaR_lower_95), linetype = "dotted") +
  geom_hline(aes(yintercept = sigmaR_upper_95), linetype = "dotted") +
  geom_vline(xintercept = Blim, linetype = "dotted") +
  geom_vline(xintercept = Bpa, linetype = "dashed") +
  geom_point(aes(shape = residual_class), size = 2) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_bw() +
  labs(
    title = "OM recruitment residuals vs SSB",
    subtitle = "Dotted horizontal lines: sigmaR 95%; vertical lines: Blim and Bpa",
    x = "SSB",
    y = "Recruitment residual"
  )

print(p_sr_compare_sigmaR)
print(p_om_resid_ts)
print(p_om_resid_hist)
print(p_om_qq)
print(p_om_resid_ssb)

ggsave(file.path(out.dir, "08_SR_curve_OM_vs_SS3_sigmaR.png"),
       p_sr_compare_sigmaR, width = 8, height = 5)

ggsave(file.path(out.dir, "09_OM_residuals_time_series_sigmaR.png"),
       p_om_resid_ts, width = 8, height = 5)

ggsave(file.path(out.dir, "10_OM_residual_distribution_vs_sigmaR.png"),
       p_om_resid_hist, width = 8, height = 5)

ggsave(file.path(out.dir, "11_OM_residual_QQplot.png"),
       p_om_qq, width = 8, height = 5)

ggsave(file.path(out.dir, "12_OM_residuals_vs_SSB_sigmaR.png"),
       p_om_resid_ssb, width = 8, height = 5)

png(file.path(out.dir, "13_OM_residual_ACF.png"), width = 900, height = 700)
acf(rec_diag_om$rec_resid_om, main = "ACF of OM recruitment residuals")
dev.off()

png(file.path(out.dir, "14_OM_residual_PACF.png"), width = 900, height = 700)
pacf(rec_diag_om$rec_resid_om, main = "PACF of OM recruitment residuals")
dev.off()

# ============================================================
# 12. Save outputs
# ============================================================

save(
  fit2_bh,
  mod2_bh,
  sr_obs,
  sr_parameter_table,
  om_curve,
  ss3_curve,
  rec_diag_om,
  om_resid_summary,
  om_resid_tests,
  extreme_residual_years,
  compare_at_obs,
  summary_compare,
  R0_ss3,
  h_ss3,
  SSB0_ss3,
  sigmaR_ss3,
  sigmaR_om,
  sigmaR_source,
  a_om,
  b_om,
  file = file.path(data.dir, "recruitment_conditioning_diagnostics_OM_vs_SS3.RData")
)

message("Saved: ", file.path(data.dir, "recruitment_conditioning_diagnostics_OM_vs_SS3.RData"))