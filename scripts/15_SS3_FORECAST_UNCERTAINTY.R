# ============================================================
# 15_SS3_FORECAST_UNCERTAINTY.R
#
# Characterisation of within-assessment uncertainty propagated
# to the SS3 short-term forecast for the probabilistic
# Fcap + Besc management procedure.
#
# IMPORTANT:
# This uncertainty is distinct from retrospective assessment error.
# Retrospective assessment error is parameterised separately in
# 14_build_assessment_error.R.
#
# OBJECTIVES:
# 1. Characterise historical variation in SS3 within-assessment uncertainty.
# 2. Evaluate propagation of terminal recruitment uncertainty to
#    first-year forecast SSB uncertainty.
# 3. Characterise contemporary forecast uncertainty using the
#    100 bootstrap SS3 models used to condition the OMs.
# 4. Examine the internal SS3 covariance structure supporting the
#    R_terminal -> SSB_forecast relationship.
# 5. Define a small set of candidate implementations for the
#    probabilistic Fcap + Besc management procedure.
#
# Historical assessments are used primarily as a diagnostic of
# uncertainty propagation through time. They are NOT assumed to
# represent a stationary distribution of future assessment precision.
#
# Bootstrap assessments represent alternative conditionings under
# the contemporary assessment configuration and are therefore used
# to characterise plausible contemporary forecast precision.
# ============================================================

rm(list=ls())

library(r4ss)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(here)

# ============================================================
# PART I. DEFINITIONS AND ANALYSIS SETTINGS
# ============================================================

# ============================================================
# 1. DIRECTORIES
# ============================================================

retro_dir <- here("data","retro_moving")
rds_dir <- here("data","Rdata")
boot_dir <- here("data","bootstrap")
runs_boot_dir <- file.path(boot_dir,"bootstrap_runs")

forecast_hist_dir <- here("data","forecast_uncertainty")
forecast_boot_dir <- file.path(forecast_hist_dir,"bootstrap_OMs")

output_dir <- here("outputs","forecast_uncertainty")
table_dir <- file.path(output_dir,"tables")
figure_dir <- file.path(output_dir,"figures")

dir.create(forecast_hist_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(forecast_boot_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

# ============================================================
# 2. ANALYSIS AND FORECAST SETTINGS
# ============================================================
# All assessments are projected under the same diagnostic fishing
# scenario to isolate variation in SS3 forecast uncertainty from
# variation in historical fishing assumptions.
#
# F = 1 is a diagnostic scenario only. It is NOT a management
# procedure and is not used as an Fcap value.
# ============================================================

base_year <- 2024L
first_year <- 2004L
Besc <- 6561

propf <- c(0.138,0.417,0.334,0.111)
Ftest <- 1
Fseason <- Ftest*propf

stopifnot(abs(sum(propf)-1)<1e-10)
stopifnot(abs(sum(Fseason)-Ftest)<1e-10)

# ============================================================
# PART II. HISTORICAL SS3 ASSESSMENTS
# ============================================================

# ============================================================
# 3. IDENTIFY HISTORICAL ASSESSMENTS
# ============================================================
# Moving retrospective models are treated here as 21 historical
# fitted SS3 assessments with terminal years 2004-2024.
#
# Retrospective revisions are NOT analysed here. For each fitted
# assessment we use the model-based StdDev reported by SS3.
# ============================================================

peels <- 0:-(base_year-first_year)
retro_dirs <- file.path(retro_dir,paste0("retro",peels))
terminal_years <- base_year+peels

retro_check <- data.frame(terminal_year=terminal_years,
                          dir=retro_dirs,
                          report=file.exists(file.path(retro_dirs,"Report.sso")),
                          covar=file.exists(file.path(retro_dirs,"covar.sso")),
                          forecast=file.exists(file.path(retro_dirs,"forecast.ss")),
                          exe=file.exists(file.path(retro_dirs,"ss3_linux")))

stopifnot(nrow(retro_check)==21)
stopifnot(all(retro_check$report))
stopifnot(all(retro_check$covar))
stopifnot(all(retro_check$forecast))
stopifnot(all(retro_check$exe))

retro_check

# ============================================================
# 4. PREPARE STANDARDISED HISTORICAL FORECASTS
# ============================================================

forecast_dirs <- file.path(forecast_hist_dir,paste0("assessment_",terminal_years))
files_ss3 <- c("control.SS","data.SS","forecast.ss","ss3_linux","starter.ss","wtatage.ss")

for(i in seq_along(retro_dirs)) {
  dir.create(forecast_dirs[i],recursive=TRUE,showWarnings=FALSE)
  file.copy(file.path(retro_dirs[i],files_ss3),file.path(forecast_dirs[i],files_ss3),overwrite=TRUE)
  Sys.chmod(file.path(forecast_dirs[i],"ss3_linux"),mode="0755")
}

for(i in seq_along(forecast_dirs)) {
  terminal_year <- terminal_years[i]
  forecast_years <- (terminal_year+1):(terminal_year+2)
  fore <- r4ss::SS_readforecast(file.path(forecast_dirs[i],"forecast.ss"),verbose=FALSE)
  fore$Nforecastyrs <- 2
  fore$InputBasis <- 99
  fore$fcast_rec_option <- 0
  fore$fcast_rec_val <- 1
  fore$ForeCatch <- data.frame(Year=rep(forecast_years,each=4),Seas=rep(1:4,2),Fleet=rep(1:4,2),F=rep(Fseason,2))
  r4ss::SS_writeforecast(fore,dir=forecast_dirs[i],file="forecast.ss",overwrite=TRUE,verbose=FALSE)
}

# ============================================================
# 5. VALIDATE HISTORICAL FORECAST CONFIGURATION
# ============================================================

forecast_config_check <- map_dfr(seq_along(forecast_dirs),function(i) {
  fore <- r4ss::SS_readforecast(file.path(forecast_dirs[i],"forecast.ss"),verbose=FALSE)
  data.frame(terminal_year=terminal_years[i],Nforecastyrs=fore$Nforecastyrs,
             InputBasis=fore$InputBasis,
             rec_option=fore$fcast_rec_option,
             first_forecast_year=min(fore$ForeCatch$year),
             last_forecast_year=max(fore$ForeCatch$year),
             n_rows=nrow(fore$ForeCatch),
             sum_F_year1=sum(fore$ForeCatch$catch_or_F[fore$ForeCatch$year==terminal_years[i]+1]),
             sum_F_year2=sum(fore$ForeCatch$catch_or_F[fore$ForeCatch$year==terminal_years[i]+2]))
})

stopifnot(all(forecast_config_check$Nforecastyrs==2))
stopifnot(all(forecast_config_check$InputBasis==99))
stopifnot(all(forecast_config_check$rec_option==0))
stopifnot(all(forecast_config_check$first_forecast_year==forecast_config_check$terminal_year+1))
stopifnot(all(forecast_config_check$last_forecast_year==forecast_config_check$terminal_year+2))
stopifnot(all(forecast_config_check$n_rows==8))
stopifnot(all(abs(forecast_config_check$sum_F_year1-Ftest)<1e-10))
stopifnot(all(abs(forecast_config_check$sum_F_year2-Ftest)<1e-10))

forecast_config_check

# ============================================================
# 6. RUN OR VERIFY HISTORICAL FORECASTS
# ============================================================
# Existing complete forecasts are retained. SS3 is run only when
# Report.sso, covar.sso or ss3.cor is missing.
# ============================================================

forecast_run_status <- map_dfr(seq_along(forecast_dirs),function(i) {
  files_complete <- all(file.exists(file.path(forecast_dirs[i],c("Report.sso","covar.sso","ss3.cor"))))
  if(!files_complete) {
    oldwd <- getwd()
    setwd(forecast_dirs[i])
    run_output <- system2("./ss3_linux",stdout=TRUE,stderr=TRUE)
    setwd(oldwd)
    completed <- any(grepl("Run has completed",run_output))
  } else {
    completed <- TRUE
  }
  data.frame(terminal_year=terminal_years[i],
             completed=completed,
             report=file.exists(file.path(forecast_dirs[i],"Report.sso")),
             covar=file.exists(file.path(forecast_dirs[i],"covar.sso")),
             cor=file.exists(file.path(forecast_dirs[i],"ss3.cor")))
})

stopifnot(all(forecast_run_status$completed))
stopifnot(all(forecast_run_status$report))
stopifnot(all(forecast_run_status$covar))
stopifnot(all(forecast_run_status$cor))

forecast_run_status

# ============================================================
# 7. VALIDATE 2024 FORECAST AGAINST THE OPERATIONAL BENCHMARK
# ============================================================
# The 2024 assessment is used to confirm that the standardised
# forecast reproduces the expected SS3 derived quantities and
# their model-based standard deviations.
# ============================================================

base_forecast_dir <- forecast_dirs[terminal_years==base_year]
ss3_base_forecast <- r4ss::SS_output(dir=base_forecast_dir,covar=TRUE,verbose=FALSE)
dq_base <- ss3_base_forecast$derived_quants

base_validation <- dq_base %>% 
  filter(Label %in% c("SSB_2024","SSB_2025","SSB_2026","Recr_2024","Recr_2025","Recr_2026")) %>% 
  select(Label,Value,StdDev) %>% mutate(CV=StdDev/Value)

base_validation

base_SSB <- base_validation %>% filter(Label=="SSB_2025")
base_R <- base_validation %>% filter(Label=="Recr_2024")

base_SSB_CV <- base_SSB$CV
base_R_CV <- base_R$CV

base_SSB_P05 <- base_SSB$Value+qnorm(0.05)*base_SSB$StdDev
base_p_below_Besc <- pnorm(Besc,mean=base_SSB$Value,sd=base_SSB$StdDev)

base_forecast_summary <- data.frame(assessment_year=2024,
                                    forecast_year=2025,
                                    SSB=base_SSB$Value,
                                    SSB_SD=base_SSB$StdDev,
                                    SSB_CV=base_SSB_CV,
                                    SSB_P05=base_SSB_P05,
                                    p_below_Besc=base_p_below_Besc,
                                    R_terminal=base_R$Value,
                                    R_terminal_SD=base_R$StdDev,
                                    R_terminal_CV=base_R_CV)

base_forecast_summary

stopifnot(abs(base_SSB$Value-13929.3)<0.1)
stopifnot(abs(base_SSB$StdDev-3336.83)<0.1)

# ============================================================
# 8. EXTRACT HISTORICAL WITHIN-ASSESSMENT UNCERTAINTY
# ============================================================
# R_terminal_CV is the model-based uncertainty in terminal
# recruitment within each fitted assessment.
#
# SSB_CV is the model-based uncertainty propagated by SS3 to
# spawning biomass in the first forecast year.
# ============================================================

forecast_uncertainty <- map_dfr(seq_along(forecast_dirs),function(i) {
  terminal_year <- terminal_years[i]
  forecast_year <- terminal_year+1
  out <- r4ss::SS_output(dir=forecast_dirs[i],covar=TRUE,verbose=FALSE)
  dq <- out$derived_quants
  ssb <- dq[dq$Label==paste0("SSB_",forecast_year),]
  rec <- dq[dq$Label==paste0("Recr_",terminal_year),]
  data.frame(assessment_year=terminal_year,forecast_year=forecast_year,SSB=ssb$Value,SSB_SD=ssb$StdDev,SSB_CV=ssb$StdDev/ssb$Value,R_terminal=rec$Value,R_terminal_SD=rec$StdDev,R_terminal_CV=rec$StdDev/rec$Value)
}) %>% arrange(assessment_year)

stopifnot(nrow(forecast_uncertainty)==21)
stopifnot(all(is.finite(as.matrix(forecast_uncertainty[,-c(1,2)]))))
stopifnot(all(forecast_uncertainty$SSB>0))
stopifnot(all(forecast_uncertainty$SSB_SD>0))
stopifnot(all(forecast_uncertainty$R_terminal>0))
stopifnot(all(forecast_uncertainty$R_terminal_SD>0))

# ============================================================
# 9. HISTORICAL DIAGNOSTIC OF UNCERTAINTY PROPAGATION
# ============================================================
# Historical assessments test whether the relationship between
# terminal recruitment uncertainty and first-year forecast SSB
# uncertainty is a consistent feature of the SS3 assessment.
#
# These assessments are NOT used as a stationary distribution
# of future assessment precision.
# ============================================================

forecast_uncertainty <- forecast_uncertainty %>% mutate(CV_ratio=SSB_CV/R_terminal_CV)

hist_cor_pearson <- cor(forecast_uncertainty$R_terminal_CV,forecast_uncertainty$SSB_CV,method="pearson")
hist_cor_spearman <- cor(forecast_uncertainty$R_terminal_CV,forecast_uncertainty$SSB_CV,method="spearman")

fit_hist_cv <- lm(SSB_CV~R_terminal_CV,data=forecast_uncertainty)
fit_hist_cv_0 <- lm(SSB_CV~0+R_terminal_CV,data=forecast_uncertainty)

hist_uncertainty_summary <- forecast_uncertainty %>% 
  summarise(n=dplyr::n(),
            R_CV_mean=mean(R_terminal_CV),
            R_CV_median=median(R_terminal_CV),
            SSB_CV_mean=mean(SSB_CV),
            SSB_CV_median=median(SSB_CV),
            CV_ratio_median=median(CV_ratio),
            pearson=hist_cor_pearson,
            spearman=hist_cor_spearman)

hist_uncertainty_summary

# ============================================================
# 10. TEMPORAL STRUCTURE OF HISTORICAL UNCERTAINTY
# ============================================================
# Temporal patterns are examined to determine whether historical
# CVs can reasonably be treated as exchangeable observations.
#
# Changes in the quantity and composition of assessment information
# through time mean that temporal trends are interpreted as changes
# in assessment precision, not as a stationary stochastic process.
# ============================================================

hist_year_cor_pearson <- cor(forecast_uncertainty$assessment_year,forecast_uncertainty$SSB_CV,method="pearson")
hist_year_cor_spearman <- cor(forecast_uncertainty$assessment_year,forecast_uncertainty$SSB_CV,method="spearman")

acf_SSB_CV <- acf(forecast_uncertainty$SSB_CV,plot=FALSE)
acf_R_CV <- acf(forecast_uncertainty$R_terminal_CV,plot=FALSE)

hist_lag1_SSB <- as.numeric(acf_SSB_CV$acf[2])
hist_lag1_R <- as.numeric(acf_R_CV$acf[2])

hist_cor_Rlevel_RCV <- cor(forecast_uncertainty$R_terminal,forecast_uncertainty$R_terminal_CV,method="pearson")
hist_cor_SSBlevel_SSBCV <- cor(forecast_uncertainty$SSB,forecast_uncertainty$SSB_CV,method="pearson")

forecast_uncertainty <- forecast_uncertainty %>% 
  mutate(period=cut(assessment_year,
                    breaks=c(2003,2012,2018,2024),
                    labels=c("2004-2012","2013-2018","2019-2024")))

historical_period_summary <- forecast_uncertainty %>% group_by(period) %>% 
  summarise(n=dplyr::n(),
            SSB_CV_mean=mean(SSB_CV),
            SSB_CV_median=median(SSB_CV),
            SSB_CV_min=min(SSB_CV),
            SSB_CV_max=max(SSB_CV),
            R_CV_mean=mean(R_terminal_CV),
            R_CV_median=median(R_terminal_CV),.groups="drop")

historical_period_summary

# ============================================================
# PART III. CONTEMPORARY BOOTSTRAP CONDITIONINGS
# ============================================================

# ============================================================
# 11. IDENTIFY THE 100 BOOTSTRAP SS3 MODELS
# ============================================================
# The 100 conditioned OMs originate from 100 validated bootstrap
# refits of the 2024 SS3 assessment.
#
# Each bootstrap refit retains its own Hessian-derived uncertainty.
# The bootstrap-to-OM mapping is retained for traceability, but this
# does not imply that the same bootstrap-specific CV must later be
# assigned to the corresponding OM in the assessment emulator.
# ============================================================

iteration_map <- read.csv(file.path(rds_dir,"iteration_map.csv"))
boot_models_100 <- file.path(runs_boot_dir,iteration_map$bootstrap)

bootstrap_model_check <- data.frame(iter=iteration_map$iter,
                                    bootstrap=iteration_map$bootstrap,
                                    dir=boot_models_100,
                                    report=file.exists(file.path(boot_models_100,"Report.sso")),
                                    covar=file.exists(file.path(boot_models_100,"covar.sso")),
                                    cor=file.exists(file.path(boot_models_100,"ss3.cor")),
                                    forecast=file.exists(file.path(boot_models_100,"forecast.ss")),
                                    exe=file.exists(file.path(boot_models_100,"ss3_linux")))

stopifnot(nrow(bootstrap_model_check)==100)
stopifnot(all(bootstrap_model_check$report))
stopifnot(all(bootstrap_model_check$covar))
stopifnot(all(bootstrap_model_check$cor))
stopifnot(all(bootstrap_model_check$forecast))
stopifnot(all(bootstrap_model_check$exe))

# ============================================================
# 12. EXTRACT TERMINAL RECRUITMENT UNCERTAINTY
# ============================================================

om_R_uncertainty <- map_dfr(seq_along(boot_models_100),function(i) {
  out <- r4ss::SS_output(dir=boot_models_100[i],covar=TRUE,verbose=FALSE)
  dq <- out$derived_quants
  rec <- dq[dq$Label=="Recr_2024",c("Value","StdDev")]
  data.frame(iter=iteration_map$iter[i],
             bootstrap=iteration_map$bootstrap[i],
             R_2024=rec$Value,
             R_2024_SD=rec$StdDev,
             R_2024_CV=rec$StdDev/rec$Value)
})

stopifnot(nrow(om_R_uncertainty)==100)
stopifnot(all(is.finite(om_R_uncertainty$R_2024_CV)))
stopifnot(all(om_R_uncertainty$R_2024_CV>0))

# ============================================================
# 13. PREPARE STANDARDISED BOOTSTRAP FORECASTS
# ============================================================
# The same diagnostic forecast used for historical assessments is
# applied to every bootstrap model: F = 1 with identical seasonal
# allocation and two forecast years.
#
# Forecasts are run in copies so the original bootstrap models used
# for OM conditioning remain unchanged.
# ============================================================

om_forecast_dirs <- file.path(forecast_boot_dir,iteration_map$bootstrap)

for(i in seq_along(boot_models_100)) {
  dir.create(om_forecast_dirs[i],recursive=TRUE,showWarnings=FALSE)
  missing_forecast <- !all(file.exists(file.path(om_forecast_dirs[i],c("Report.sso","covar.sso","ss3.cor"))))
  if(missing_forecast) file.copy(file.path(boot_models_100[i],files_ss3),file.path(om_forecast_dirs[i],files_ss3),overwrite=TRUE)
  Sys.chmod(file.path(om_forecast_dirs[i],"ss3_linux"),mode="0755")
}

for(i in seq_along(om_forecast_dirs)) {
  output_complete <- all(file.exists(file.path(om_forecast_dirs[i],c("Report.sso","covar.sso","ss3.cor"))))
  if(!output_complete) {
    fore <- r4ss::SS_readforecast(file.path(om_forecast_dirs[i],"forecast.ss"),verbose=FALSE)
    fore$Nforecastyrs <- 2
    fore$InputBasis <- 99
    fore$fcast_rec_option <- 0
    fore$fcast_rec_val <- 1
    fore$ForeCatch <- data.frame(Year=rep(2025:2026,each=4),Seas=rep(1:4,2),Fleet=rep(1:4,2),F=rep(Fseason,2))
    r4ss::SS_writeforecast(fore,dir=om_forecast_dirs[i],file="forecast.ss",overwrite=TRUE,verbose=FALSE)
  }
}

# ============================================================
# 14. VALIDATE BOOTSTRAP FORECAST CONFIGURATION
# ============================================================

bootstrap_forecast_config <- map_dfr(seq_along(om_forecast_dirs),function(i) {
  fore <- r4ss::SS_readforecast(file.path(om_forecast_dirs[i],"forecast.ss"),verbose=FALSE)
  data.frame(iter=iteration_map$iter[i],
             bootstrap=iteration_map$bootstrap[i],
             Nforecastyrs=fore$Nforecastyrs,
             InputBasis=fore$InputBasis,
             rec_option=fore$fcast_rec_option,
             n_rows=nrow(fore$ForeCatch),
             sum_F_2025=sum(fore$ForeCatch$catch_or_F[fore$ForeCatch$year==2025]),
             sum_F_2026=sum(fore$ForeCatch$catch_or_F[fore$ForeCatch$year==2026]))
})

stopifnot(all(bootstrap_forecast_config$Nforecastyrs==2))
stopifnot(all(bootstrap_forecast_config$InputBasis==99))
stopifnot(all(bootstrap_forecast_config$rec_option==0))
stopifnot(all(bootstrap_forecast_config$n_rows==8))
stopifnot(all(abs(bootstrap_forecast_config$sum_F_2025-Ftest)<1e-10))
stopifnot(all(abs(bootstrap_forecast_config$sum_F_2026-Ftest)<1e-10))

# ============================================================
# 15. RUN OR VERIFY BOOTSTRAP FORECASTS
# ============================================================
# Existing complete forecasts are reused. This prevents unnecessary
# repetition of the 100 SS3 runs.
# ============================================================

bootstrap_run_status <- map_dfr(seq_along(om_forecast_dirs),function(i) {
  files_complete <- all(file.exists(file.path(om_forecast_dirs[i],c("Report.sso","covar.sso","ss3.cor"))))
  if(!files_complete) {
    oldwd <- getwd()
    setwd(om_forecast_dirs[i])
    run_output <- system2("./ss3_linux",stdout=TRUE,stderr=TRUE)
    setwd(oldwd)
    completed <- any(grepl("Run has completed",run_output))
  } else {
    completed <- TRUE
  }
  data.frame(iter=iteration_map$iter[i],
             bootstrap=iteration_map$bootstrap[i],
             completed=completed,report=file.exists(file.path(om_forecast_dirs[i],"Report.sso")),
             covar=file.exists(file.path(om_forecast_dirs[i],"covar.sso")),
             cor=file.exists(file.path(om_forecast_dirs[i],"ss3.cor")))
})

stopifnot(all(bootstrap_run_status$completed))
stopifnot(all(bootstrap_run_status$report))
stopifnot(all(bootstrap_run_status$covar))
stopifnot(all(bootstrap_run_status$cor))

# ============================================================
# 16. EXTRACT CONTEMPORARY FORECAST SSB UNCERTAINTY
# ============================================================

om_SSB_uncertainty <- map_dfr(seq_along(om_forecast_dirs),function(i) {
  out <- r4ss::SS_output(dir=om_forecast_dirs[i],covar=TRUE,verbose=FALSE)
  dq <- out$derived_quants
  ssb <- dq[dq$Label=="SSB_2025",c("Value","StdDev")]
  data.frame(iter=iteration_map$iter[i],
             bootstrap=iteration_map$bootstrap[i],
             SSB_2025=ssb$Value,SSB_2025_SD=ssb$StdDev,
             SSB_2025_CV=ssb$StdDev/ssb$Value)
})

stopifnot(nrow(om_SSB_uncertainty)==100)
stopifnot(all(om_SSB_uncertainty$iter==iteration_map$iter))
stopifnot(all(om_SSB_uncertainty$bootstrap==iteration_map$bootstrap))
stopifnot(all(is.finite(om_SSB_uncertainty$SSB_2025_CV)))
stopifnot(all(om_SSB_uncertainty$SSB_2025_CV>0))

# ============================================================
# 17. LINK TERMINAL R AND FIRST-YEAR FORECAST SSB UNCERTAINTY
# ============================================================

om_forecast_uncertainty <- left_join(om_R_uncertainty,
                                     om_SSB_uncertainty,by=c("iter","bootstrap")) %>% 
  mutate(CV_ratio=SSB_2025_CV/R_2024_CV)

stopifnot(nrow(om_forecast_uncertainty)==100)

boot_cor_pearson <- cor(om_forecast_uncertainty$R_2024_CV,
                        om_forecast_uncertainty$SSB_2025_CV,method="pearson")
boot_cor_spearman <- cor(om_forecast_uncertainty$R_2024_CV,
                         om_forecast_uncertainty$SSB_2025_CV,method="spearman")

fit_boot_cv <- lm(SSB_2025_CV~R_2024_CV,data=om_forecast_uncertainty)
fit_boot_cv_0 <- lm(SSB_2025_CV~0+R_2024_CV,data=om_forecast_uncertainty)

bootstrap_uncertainty_summary <- om_forecast_uncertainty %>% 
  summarise(n=dplyr::n(),
            R_CV_mean=mean(R_2024_CV),
            R_CV_median=median(R_2024_CV),
            SSB_CV_mean=mean(SSB_2025_CV),
            SSB_CV_median=median(SSB_2025_CV),
            SSB_CV_P05=quantile(SSB_2025_CV,0.05),
            SSB_CV_P95=quantile(SSB_2025_CV,0.95),
            CV_ratio_median=median(CV_ratio),
            pearson=boot_cor_pearson,
            spearman=boot_cor_spearman)

bootstrap_uncertainty_summary

# ============================================================
# PART IV. INTERNAL SS3 UNCERTAINTY STRUCTURE
# ============================================================

# ============================================================
# 18. DIRECT WITHIN-ASSESSMENT R2024 -> SSB2025 RELATIONSHIP
# ============================================================
# SS3 CoVar contains correlations among parameters and derived
# quantities from the model covariance structure.
#
# This provides a direct within-assessment diagnostic that is
# distinct from correlations calculated across historical or
# bootstrap assessments.
# ============================================================

cor_R2024_SSB2025 <- ss3_base_forecast$CoVar %>% 
  filter((label.i=="Recr_2024" & label.j=="SSB_2025") | (label.i=="SSB_2025" & label.j=="Recr_2024"))

cor_dev_R2024 <- ss3_base_forecast$CoVar %>% 
  filter((label.i=="Main_RecrDev_2024" & label.j=="Recr_2024") | (label.i=="Recr_2024" & label.j=="Main_RecrDev_2024"))

cor_dev_SSB2025 <- ss3_base_forecast$CoVar %>% 
  filter((label.i=="Main_RecrDev_2024" & label.j=="SSB_2025") | (label.i=="SSB_2025" & label.j=="Main_RecrDev_2024"))

internal_uncertainty_correlations <- bind_rows(cor_R2024_SSB2025,cor_dev_R2024,cor_dev_SSB2025) %>% 
  select(Par..i,Par..j,label.i,label.j,corr)

internal_uncertainty_correlations

stopifnot(nrow(cor_R2024_SSB2025)==1)

# ============================================================
# PART V. EVIDENCE SYNTHESIS
# ============================================================

# ============================================================
# 19. HISTORICAL VS CONTEMPORARY UNCERTAINTY
# ============================================================
# Historical assessments demonstrate consistency in the propagation
# of terminal recruitment uncertainty to first-year forecast SSB.
#
# Contemporary bootstrap conditionings characterise the magnitude
# and between-conditioning variation in forecast precision under
# the current assessment configuration.
# ============================================================

Table15_1 <- bind_rows(
  data.frame(Evidence="Historical assessments",
             n=nrow(forecast_uncertainty),
             R_CV=mean(forecast_uncertainty$R_terminal_CV),
             SSB_CV=mean(forecast_uncertainty$SSB_CV),
             Pearson_r=hist_cor_pearson,
             Spearman_r=hist_cor_spearman,
             Median_CV_ratio=median(forecast_uncertainty$CV_ratio)),
  data.frame(Evidence="Contemporary bootstrap",
             n=nrow(om_forecast_uncertainty),
             R_CV=median(om_forecast_uncertainty$R_2024_CV),
             SSB_CV=median(om_forecast_uncertainty$SSB_2025_CV),
             Pearson_r=boot_cor_pearson,
             Spearman_r=boot_cor_spearman,
             Median_CV_ratio=median(om_forecast_uncertainty$CV_ratio)),
  data.frame(Evidence="Base assessment 2024",
             n=1,
             R_CV=base_R_CV,
             SSB_CV=base_SSB_CV,
             Pearson_r=NA,
             Spearman_r=NA,
             Median_CV_ratio=base_SSB_CV/base_R_CV))

Table15_1

# ============================================================
# 20. BASE ASSESSMENT RELATIVE TO CONTEMPORARY BOOTSTRAPS
# ============================================================

base_SSB_percentile   <- mean(om_forecast_uncertainty$SSB_2025_CV<=base_SSB_CV)
base_R_percentile     <- mean(om_forecast_uncertainty$R_2024_CV<=base_R_CV)
base_ratio            <- base_SSB_CV/base_R_CV
base_ratio_percentile <- mean(om_forecast_uncertainty$CV_ratio<=base_ratio)

base_vs_bootstrap <- data.frame(metric=c("R_terminal_CV","SSB_forecast_CV","SSB_CV_to_R_CV_ratio"),
                                base=c(base_R_CV,base_SSB_CV,base_ratio),
                                bootstrap_percentile=c(base_R_percentile,base_SSB_percentile,base_ratio_percentile))

base_vs_bootstrap

# ============================================================
# 21. CANDIDATE IMPLEMENTATION SCENARIOS
# ============================================================
# S0 provides the deterministic benchmark.
#
# S1 assumes that the precision estimated for the current 2024
# assessment remains constant throughout the projection.
#
# S2 recognises uncertainty in contemporary assessment precision.
# A CV is sampled from the empirical distribution obtained from
# the 100 current bootstrap SS3 conditionings and held constant
# within a trajectory.
#
# S2 does NOT automatically pair OM_i with bootstrap CV_i because
# OM conditioning uncertainty and assessment-emulator precision
# represent conceptually distinct components of the MSE.
#
# No annual iid resampling is assumed because the bootstrap
# distribution does not estimate temporal variation in precision.
# ============================================================

uncertainty_scenarios <- data.frame(
  scenario=c("S0","S1","S2"),
  name=c("Deterministic","Current assessment","Contemporary uncertainty"),
  CV_source=c("None","Base 2024 SS3 forecast","100 contemporary bootstrap SS3 forecasts"),
  implementation=c("CV = 0",paste0("CV = ",round(base_SSB_CV,5)),"CV_i sampled from empirical bootstrap distribution"),
  temporal_assumption=c("Not applicable","Constant contemporary precision","CV_i fixed within trajectory"),
  rationale=c("Deterministic benchmark",
              "Current SS3 estimate of first-year forecast precision",
              "Plausible variation in forecast precision under contemporary assessment conditioning"))

Table15_2 <- uncertainty_scenarios
Table15_2

# ============================================================
# PART VI. REPORT TABLES AND FIGURES
# ============================================================

# ============================================================
# 22. SAVE REPORT TABLES
# ============================================================

write.csv(Table15_1,file.path(table_dir,"Table15_1_uncertainty_evidence.csv"),row.names=FALSE)
write.csv(Table15_2,file.path(table_dir,"Table15_2_candidate_scenarios.csv"),row.names=FALSE)
write.csv(historical_period_summary,file.path(table_dir,"Table15_S1_historical_periods.csv"),row.names=FALSE)
write.csv(base_vs_bootstrap,file.path(table_dir,"Table15_S2_base_vs_bootstrap.csv"),row.names=FALSE)
write.csv(internal_uncertainty_correlations,file.path(table_dir,"Table15_S3_internal_SS3_correlations.csv"),row.names=FALSE)

# ============================================================
# 23. FIGURE 15.1 - HISTORICAL FORECAST UNCERTAINTY
# ============================================================
# Purpose:
# Show that assessment precision changed substantially through
# time, supporting the decision not to treat the 21 historical
# assessments as a stationary distribution of future precision.
# ============================================================

fig15_1_data <- forecast_uncertainty %>% 
  select(assessment_year,R_terminal_CV,SSB_CV) %>% 
  pivot_longer(cols=c(R_terminal_CV,SSB_CV),names_to="quantity",values_to="CV") %>% 
  mutate(quantity=recode(quantity,R_terminal_CV="Terminal recruitment",SSB_CV="First-year forecast SSB"))

fig15_1 <- ggplot(fig15_1_data,aes(x=assessment_year,y=CV,linetype=quantity,shape=quantity)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(2004,2024,4)) +
  labs(x="Assessment year",y="Coefficient of variation",linetype=NULL,shape=NULL) +
  theme_bw() +
  theme(legend.position="top",panel.grid.minor=element_blank())

ggsave(file.path(figure_dir,"Fig15_1_historical_forecast_uncertainty.png"),fig15_1,width=8,height=5,dpi=300)
ggsave(file.path(figure_dir,"Fig15_1_historical_forecast_uncertainty.pdf"),fig15_1,width=8,height=5)

# ============================================================
# 24. FIGURE 15.2 - R -> FORECAST SSB UNCERTAINTY PROPAGATION
# ============================================================
# Purpose:
# Demonstrate that the strong relationship between uncertainty in
# terminal recruitment and first-year forecast SSB occurs both
# across historical assessments and across contemporary bootstrap
# conditionings.
#
# Regression lines are diagnostic only and are NOT used directly
# to parameterise the shortcut MSE.
# ============================================================

fig15_2_hist <- forecast_uncertainty %>% transmute(R_CV=R_terminal_CV,SSB_CV=SSB_CV,source="Historical assessments")
fig15_2_boot <- om_forecast_uncertainty %>% transmute(R_CV=R_2024_CV,SSB_CV=SSB_2025_CV,source="Contemporary bootstrap")
fig15_2_base <- data.frame(R_CV=base_R_CV,SSB_CV=base_SSB_CV,source="Base assessment 2024")
fig15_2_data <- bind_rows(fig15_2_hist,fig15_2_boot)

fig15_2 <- ggplot(fig15_2_data,aes(x=R_CV,y=SSB_CV,shape=source)) +
  geom_point(size=2,alpha=0.7) +
  geom_smooth(aes(linetype=source),method="lm",se=FALSE,linewidth=0.7) +
  geom_point(data=fig15_2_base,aes(x=R_CV,y=SSB_CV),inherit.aes=FALSE,shape=8,size=4) +
  labs(x="Terminal recruitment CV",y="First-year forecast SSB CV",shape=NULL,linetype=NULL) +
  theme_bw() +
  theme(legend.position="top",panel.grid.minor=element_blank())

ggsave(file.path(figure_dir,"Fig15_2_R_SSB_uncertainty_propagation.png"),fig15_2,width=7,height=6,dpi=300)
ggsave(file.path(figure_dir,"Fig15_2_R_SSB_uncertainty_propagation.pdf"),fig15_2,width=7,height=6)

# ============================================================
# 25. FIGURE 15.3 - CONTEMPORARY FORECAST PRECISION
# ============================================================
# Purpose:
# Characterise the empirical distribution of first-year forecast
# SSB CV across the 100 contemporary bootstrap conditionings and
# show the position of the base 2024 assessment.
# ============================================================

fig15_3 <- ggplot(om_forecast_uncertainty,aes(x=SSB_2025_CV)) +
  geom_histogram(bins=20,fill="grey80",colour="black") +
  geom_vline(xintercept=base_SSB_CV,linetype="dashed",linewidth=0.8) +
  geom_vline(xintercept=median(om_forecast_uncertainty$SSB_2025_CV),linetype="dotted",linewidth=0.8) +
  annotate("text",x=base_SSB_CV,y=Inf,label=paste0("Base = ",round(base_SSB_CV,3)),angle=90,vjust=-0.5,hjust=1.1,size=3.5) +
  annotate("text",x=median(om_forecast_uncertainty$SSB_2025_CV),y=Inf,label=paste0("Median = ",round(median(om_forecast_uncertainty$SSB_2025_CV),3)),angle=90,vjust=-0.5,hjust=1.1,size=3.5) +
  labs(x="First-year forecast SSB CV",y="Number of bootstrap assessments") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())

ggsave(file.path(figure_dir,"Fig15_3_contemporary_SSB_CV_distribution.png"),fig15_3,width=7,height=5,dpi=300)
ggsave(file.path(figure_dir,"Fig15_3_contemporary_SSB_CV_distribution.pdf"),fig15_3,width=7,height=5)

# ============================================================
# 26. SAVE CHARACTERISATION OBJECT
# ============================================================
# This object contains all information required for subsequent
# implementation of forecast uncertainty in the management
# procedure. SS3 forecasts therefore do not need to be rerun.
# ============================================================

forecast_uncertainty_characterisation <- list(
  settings=list(base_year=base_year,
                first_year=first_year,
                Besc=Besc,
                Ftest=Ftest,
                propf=propf),
  historical=forecast_uncertainty,
  historical_periods=historical_period_summary,
  contemporary=om_forecast_uncertainty,
  base=base_forecast_summary,
  internal_correlations=internal_uncertainty_correlations,
  base_vs_bootstrap=base_vs_bootstrap,
  scenarios=uncertainty_scenarios,
  diagnostics=list(hist_cor_pearson=hist_cor_pearson,
                   hist_cor_spearman=hist_cor_spearman,
                   boot_cor_pearson=boot_cor_pearson,
                   boot_cor_spearman=boot_cor_spearman,
                   hist_year_cor_pearson=hist_year_cor_pearson,
                   hist_year_cor_spearman=hist_year_cor_spearman,
                   hist_lag1_SSB=hist_lag1_SSB,
                   hist_lag1_R=hist_lag1_R,
                   hist_cor_Rlevel_RCV=hist_cor_Rlevel_RCV,
                   hist_cor_SSBlevel_SSBCV=hist_cor_SSBlevel_SSBCV)
)

saveRDS(forecast_uncertainty_characterisation,file.path(rds_dir,"forecast_uncertainty_characterisation.rds"))

# ============================================================
# 27. FINAL VALIDATION
# ============================================================

stopifnot(nrow(forecast_uncertainty_characterisation$historical)==21)
stopifnot(nrow(forecast_uncertainty_characterisation$contemporary)==100)
stopifnot(nrow(forecast_uncertainty_characterisation$scenarios)==3)
stopifnot(file.exists(file.path(rds_dir,"forecast_uncertainty_characterisation.rds")))

cat("\nSS3 forecast uncertainty characterisation completed.\n")
cat("Historical assessments:",nrow(forecast_uncertainty),"\n")
cat("Contemporary bootstrap assessments:",nrow(om_forecast_uncertainty),"\n")
cat("Base SSB forecast CV:",round(base_SSB_CV,4),"\n")
cat("Bootstrap median SSB forecast CV:",round(median(om_forecast_uncertainty$SSB_2025_CV),4),"\n")
cat("Historical R-CV / SSB-CV Pearson r:",round(hist_cor_pearson,4),"\n")
cat("Bootstrap R-CV / SSB-CV Pearson r:",round(boot_cor_pearson,4),"\n")
cat("Within-assessment R2024 / SSB2025 correlation:",round(cor_R2024_SSB2025$corr,4),"\n")

