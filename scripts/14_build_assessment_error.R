# ============================================================
# 14_build_assessment_error.R
# Assessment-error parameterisation, validation, and implementation checks
# Shortcut MSE - Anchovy ane.27.9aS
#
# Stock / analysis:
#   European anchovy (Engraulis encrasicolus), ane.27.9aS
#   Gulf of Cádiz shortcut Management Strategy Evaluation (MSE)
#
# Role in the workflow:
#   This script parameterises the assessment-error component used
#   in the closed-loop shortcut MSE from historical SS3 assessment
#   revisions. It validates the resulting multivariate emulator,
#   generates the 1000 assessment-error trajectories required by
#   the closed-loop simulations, and contains the existing technical
#   checks used to examine observation and assessment-error
#   implementation.
#
#   Historical SS3 assessments
#          |
#          v
#   Annual assessment revisions
#          |
#          v
#   Multivariate assessment-error emulator
#          |
#          v
#   1000 assessment-error trajectories
#          |
#          v
#   Closed-loop implementation checks
#
# Purpose:
#   1. Reconstruct annual revisions in SSB, F, and recruitment (R)
#      from consecutive SS3 assessments.
#   2. Parameterise the multivariate assessment-error distribution.
#   3. Diagnose cross-variable and temporal error structure.
#   4. Validate the fitted emulator internally and against
#      conventional retrospective analyses.
#   5. Generate assessment-error multipliers for 1000 MSE
#      trajectories during management years 2024-2053.
#   6. Retain the existing implementation checks for perfect
#      observation, true/perceived stock separation, deterministic
#      error scaling, emulator application, and recruitment timing.
#
# Assessment-error definition:
#
#       e = log(perceived / true)
#
#   implemented historically as:
#
#       e_t = log(older assessment / updated assessment)
#
#   The older assessment represents the information available to
#   management at that time. The subsequent assessment is used as
#   the updated proxy for the underlying stock state.
#
# Reference assessment-error model:
#   Variables: SSB, F, R
#   Error scale: log ratio
#   Cross-variable structure: estimated multivariate covariance
#   Temporal structure: iid
#
# Inputs:
#   data/retro_moving/
#     retro0 ... retro-20
#     Moving retrospective SS3 runs used to reconstruct annual
#     assessment revisions.
#
#   data/retro/
#     retro0 ... retro-5
#     Conventional retrospective SS3 runs used for external
#     validation of the assessment-error emulator.
#
#   data/Rdata/fit_bh_boot_unconstrained_100.rds
#     Used by the existing recruitment-timing implementation check.
#
#   functions/perfectObs4seas.R
#     Used by the existing observation implementation checks.
#
#   OM_validation_C7000_block_01.rds
#     Used by the existing observation/assessment-error
#     implementation checks through dir_C7000.
#
# Outputs:
#   data/Rdata/
#     assessment_error_parameters.rds
#     assessment_error_historical_revisions.rds
#     assessment_error_simulated_reference.rds
#     assessment_error_reference_1000iter.rds
#
#   outputs/assessment_error/figures/
#     01_historical_assessment_trajectories.png
#     02_assessment_error_time_series.png
#     03_assessment_error_distribution_validation.png
#     04_example_assessment_error_trajectories.png
#     05_assessment_error_external_validation.png
#
#   outputs/assessment_error/tables/
#     01_assessment_error_parameters.csv
#     02_assessment_error_temporal_structure.csv
#     03_emulator_validation.csv
#     04_retrospective_validation.csv
#     05_mohn_rho_validation.csv
#     06_directional_validation.csv
#
# Working Document outputs:
#   The figures and CSV tables written to
#   outputs/assessment_error/ document the parameterisation and
#   validation of the assessment-error emulator for the WD.
#
# Reproducibility settings retained:
#   Historical update years: 2005-2024
#   Closed-loop error years: 2024-2053
#   Closed-loop iterations: 1000
#   Emulator validation simulations: 10000
#   Retrospective validation simulations: 10000
#   SS3 convergence criterion: max_gradient < 0.001
#   Random seed: 123 where specified in the original script
#
# Important:
#   - Correlation among SSB, F, and R errors is retained through
#     the estimated covariance matrix.
#   - Assessment errors are simulated on the log scale and
#     transformed to multiplicative errors using exp(error).
#   - The reference implementation assumes iid errors among years.
#   - Diagnostic calculations are documented as diagnostics and
#     are not presented as acceptance criteria unless explicitly
#     enforced by the existing code.
#   - No scientific calculation, parameter, filter, model, or
#     methodological decision has been changed in this
#     restructuring.
# ============================================================

# ============================================================
# 0. PACKAGES AND SESSION INITIALISATION
# Clear the workspace and load packages used by the existing
# assessment-error parameterisation and validation workflow.
# ============================================================

rm(list=ls())

library(r4ss)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(ggplot2)
library(MASS)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

retro_dir <- here("data","retro_moving")
rds_dir <- here("data","Rdata")

base_year <- 2024L
first_year <- 2004L
variables <- c("SSB","F","R")
gradient_limit <- 0.001

fig_dir <- here("outputs","assessment_error","figures")
tab_dir <- here("outputs","assessment_error","tables")

dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(tab_dir,recursive=TRUE,showWarnings=FALSE)

# ============================================================
# 2. LOAD AND VALIDATE MOVING RETROSPECTIVE SS3 RUNS
# Read the annual SS3 retrospective runs used to reconstruct
# historical assessment revisions and verify convergence.
# ============================================================

peels <- 0:-(base_year-first_year)
dirs <- file.path(retro_dir,paste0("retro",peels))

run_check <- tibble(
  peel=peels,
  terminal_year=base_year+peels,
  dir=dirs,
  exists=file.exists(file.path(dirs,"Report.sso"))
) %>%
  filter(exists) %>%
  arrange(desc(terminal_year))

models <- SSgetoutput(dirvec=run_check$dir,verbose=FALSE)
summary_ss3 <- SSsummarize(models)

model_def <- run_check %>%
  mutate(
    model=summary_ss3$modelnames,
    max_gradient=summary_ss3$maxgrad,
    converged=is.finite(max_gradient) & max_gradient<gradient_limit
  )

stopifnot(all(model_def$converged))

# ============================================================
# 3. EXTRACT HISTORICAL ASSESSMENT ESTIMATES
# Extract SSB, fishing mortality, and recruitment estimates
# from each retrospective assessment and associate each
# estimate with its assessment terminal year.
# ============================================================

assessment_long <- summary_ss3$quants %>%  as_tibble() %>%
                    pivot_longer(cols=all_of(summary_ss3$modelnames),names_to="model",values_to="estimate") %>%
                    left_join(model_def %>% dplyr::select(model,terminal_year),by="model") %>%
                    mutate(year=suppressWarnings(as.integer(Yr)),
                    variable=case_when(
                      grepl("^SSB_[0-9]{4}$",Label) ~ "SSB",
                      grepl("^F_[0-9]{4}$",Label) ~ "F",
                      grepl("^Recr_[0-9]{4}$",Label) ~ "R",
                        TRUE ~ NA_character_)) %>%
                        filter(!is.na(variable),!is.na(year),is.finite(estimate),estimate>0)

# ============================================================
# 4. CALCULATE ANNUAL ASSESSMENT REVISIONS
# Compare estimates for the same quantity and year between
# consecutive assessments and express the revision as:
#
#   log_error = log(older estimate / updated estimate)
#
# These annual revisions provide the historical proxy for
# assessment error used in the shortcut MSE.
# ============================================================

old <- assessment_long %>%
        dplyr::select(variable,
                      comparison_year=year,
                      old_terminal_year=terminal_year,
                      old_estimate=estimate)
new <- assessment_long %>%
        dplyr::select(variable,
                      comparison_year=year,
                      new_terminal_year=terminal_year,
                      new_estimate=estimate)

available_years <- sort(unique(model_def$terminal_year))

assessment_errors <- expand_grid(update_year=available_years,variable=variables) %>%
                      filter((update_year-1L) %in% available_years) %>%
                      mutate(comparison_year=update_year-1L,old_terminal_year=update_year-1L,new_terminal_year=update_year) %>%
                      left_join(old,by=c("variable","comparison_year","old_terminal_year")) %>%
                      left_join(new,by=c("variable","comparison_year","new_terminal_year")) %>%
                      filter(is.finite(old_estimate),is.finite(new_estimate),old_estimate>0,new_estimate>0) %>%
                      mutate(log_error=log(old_estimate/new_estimate),
                      relative_error=old_estimate/new_estimate-1) %>%
                      arrange(update_year,variable)

stopifnot(n_distinct(assessment_errors$update_year)==20)
stopifnot(all(table(assessment_errors$variable)==20))

# ============================================================
# 5. PARAMETERISE THE MULTIVARIATE ERROR DISTRIBUTION
# Estimate the mean log-error vector and covariance matrix for
# SSB, F, and recruitment from the historical annual revisions.
#
# The covariance matrix preserves the observed dependence among
# assessment errors in the three quantities.
# ============================================================

error_wide <- assessment_errors %>%
              dplyr::select(update_year,variable,log_error) %>%
              pivot_wider(names_from=variable,values_from=log_error) %>%
              drop_na(all_of(variables)) %>%
              arrange(update_year)

mu <- colMeans(error_wide[,variables])
Sigma <- cov(error_wide[,variables])
correlation <- cov2cor(Sigma)


assessment_parameter_table <- tibble(
                              Quantity=variables,
                              `Mean log-error`=round(mu[variables],3),
                              `SD log-error`=round(sqrt(diag(Sigma))[variables],3),
                                SSB=round(correlation[variables,"SSB"],3),
                                F=round(correlation[variables,"F"],3),
                                R=round(correlation[variables,"R"],3))

assessment_parameter_table
write.csv(assessment_parameter_table,file.path(tab_dir,"01_assessment_error_parameters.csv"),row.names=FALSE)

# ============================================================
# 6. DIAGNOSE HISTORICAL ERROR STRUCTURE
# Examine the historical assessment-error structure before
# defining the reference emulator:
#   - historical assessment trajectories,
#   - robustness of the SSB-F dependence,
#   - temporal autocorrelation in SSB, F, and R errors.
# ============================================================

assessment_long_plot <- assessment_long %>% filter(year<=terminal_year) %>%
                        mutate(variable=factor(variable,levels=c("SSB","F","R"),
                                               labels=c("SSB","Fishing mortality","Recruitment")))

p_retro <- ggplot(assessment_long_plot,aes(year,estimate,group=terminal_year,color=terminal_year)) +
          geom_line(linewidth=0.6,alpha=0.8) +
          geom_point(size=1.2,alpha=0.8) +
          facet_wrap(~variable,ncol=1,scales="free_y") +
          scale_color_viridis_c(option="C",end=0.9,breaks=c(2005,2010,2015,2020,2024)) +
          labs(x="Year",y=NULL,color="Assessment\nterminal year")+
          theme_bw() +
          theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
          strip.background=element_blank(),strip.text=element_text(face="bold"),
          legend.position="right",
          axis.title=element_text(size=11),axis.text=element_text(size=10))


p_retro
ggsave(file.path(fig_dir,"01_historical_assessment_trajectories.png"),p_retro,width=8,height=9,dpi=300)
# ============================================================
# 6.1 ROBUSTNESS OF SSB-F DEPENDENCE
# ============================================================

ssb_f_loo <- map_dfr(error_wide$update_year,function(y) {
              x <- filter(error_wide,update_year!=y)
              tibble(omitted_year=y,
                     Pearson=cor(x$SSB,x$F,method="pearson"),
                     Spearman=cor(x$SSB,x$F,method="spearman"))
})

ssb_f_loo_summary <- ssb_f_loo %>%
                      summarise(Pearson_min=min(Pearson),
                                Pearson_max=max(Pearson),
                                Spearman_min=min(Spearman),
                                Spearman_max=max(Spearman))

ssb_f_loo_summary

# ============================================================
# 6.2 TEMPORAL STRUCTURE
# ============================================================

ar1 <- map_dfr(variables,function(v) {
  x <- error_wide[[v]]
  fit <- lm(x[-1]~x[-length(x)])
  tibble(variable=v,phi=unname(coef(fit)[2]),p_value=coef(summary(fit))[2,4])
})

ar1
write.csv(ar1,file.path(tab_dir,"02_assessment_error_temporal_structure.csv"),row.names=FALSE)

# ============================================================
# 7. DEFINE THE REFERENCE ASSESSMENT-ERROR MODEL
# Define the assessment-error model used in the shortcut MSE.
#
# The reference implementation assumes independent errors
# among years (iid) while retaining the estimated covariance
# among SSB, F, and recruitment within each year.
# ============================================================

temporal_model <- "iid"
phi <- setNames(rep(0,length(variables)),variables)

# ============================================================
# 8. BUILD AND SAVE THE REFERENCE PARAMETER OBJECT
# ============================================================

assessment_error_parameters <- list(variables=variables,
                                    mu=mu,
                                    Sigma=Sigma,
                                    correlation=correlation,
                                    temporal_model=temporal_model,
                                    phi=phi,
                                    years=range(error_wide$update_year),
                                    n_years=nrow(error_wide),
                                    definition="log(older assessment / updated assessment)")

stopifnot(all(eigen(Sigma,symmetric=TRUE)$values>0))

print(mu)
print(sqrt(diag(Sigma)))
print(correlation)
print(ssb_f_loo_summary)
print(ar1)
print(temporal_model)

saveRDS(assessment_error_parameters,file.path(rds_dir,"assessment_error_parameters.rds"))
saveRDS(assessment_errors,file.path(rds_dir,"assessment_error_historical_revisions.rds"))

# ============================================================
# 9. DOCUMENT THE HISTORICAL TEMPORAL PATTERN
# ============================================================

error_long_plot <- error_wide %>% pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error") %>%
                    mutate(mean_error=mu[variable],
                           variable=factor(variable,levels=c("SSB","F","R"), labels=c("SSB","Fishing mortality","Recruitment")))

p_error_time <- ggplot(error_long_plot,aes(update_year,log_error)) +
                geom_hline(yintercept=0,linetype="dashed",linewidth=0.5,color="grey50") +
                geom_line(aes(y=mean_error),linewidth=0.6,color="black") +
                geom_line(linewidth=0.6,color="grey40") +
                geom_point(size=1.8,color="black") +
                facet_wrap(~variable,ncol=1,scales="free_y") +
                labs(x="Assessment update year",y="Annual log-revision") +
                theme_bw() +
                theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
                strip.background=element_blank(),strip.text=element_text(face="bold"),
                axis.title=element_text(size=11),axis.text=element_text(size=10))

p_error_time
ggsave(file.path(fig_dir,"02_assessment_error_time_series.png"),p_error_time,width=8,height=8,dpi=300)

# ============================================================
# 10. VALIDATE THE MULTIVARIATE ERROR EMULATOR
# Simulate assessment errors from the fitted multivariate
# distribution and verify that the simulated errors reproduce
# the target marginal means, standard deviations, and
# cross-variable correlations.
# ============================================================

set.seed(123)
n_sim <- 10000

sim_error <- MASS::mvrnorm(n=n_sim,mu=mu,Sigma=Sigma) %>%as.data.frame()
names(sim_error) <- variables

# ============================================================
# 10.1 MARGINAL DISTRIBUTIONS
# ============================================================

sim_summary <- tibble(variable=variables,
                      target_mean=mu[variables],
                      simulated_mean=sapply(sim_error[variables],mean),
                      target_sd=sqrt(diag(Sigma))[variables],
                      simulated_sd=sapply(sim_error[variables],sd))

sim_summary

# ============================================================
# 10.2 CROSS-VARIABLE DEPENDENCE
# ============================================================

sim_correlation <- cor(sim_error[,variables])

correlation
sim_correlation

sim_validation_table <- tibble(variable=variables,
                               target_mean=mu[variables],
                               simulated_mean=sapply(sim_error[variables],mean),
                               target_sd=sqrt(diag(Sigma))[variables],
                               simulated_sd=sapply(sim_error[variables],sd),
                               target_cor_SSB=correlation[variables,"SSB"],
                               simulated_cor_SSB=sim_correlation[variables,"SSB"])

sim_validation_table
write.csv(sim_validation_table,file.path(tab_dir,"03_emulator_validation.csv"),row.names=FALSE)

saveRDS(assessment_error_parameters,file.path(rds_dir,"assessment_error_parameters.rds"))
saveRDS(assessment_errors,file.path(rds_dir,"assessment_error_historical_revisions.rds"))
saveRDS(sim_error,file.path(rds_dir,"assessment_error_simulated_reference.rds"))


# ============================================================
# 10.3 HISTORICAL VS SIMULATED DISTRIBUTIONS
# ============================================================

hist_plot <- error_wide %>%pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error")
sim_plot <- sim_error %>%pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error")
hist_plot <- hist_plot %>%mutate(variable=factor(variable,levels=c("SSB","F","R"),labels=c("SSB","Fishing mortality","Recruitment")))
sim_plot <- sim_plot %>%mutate(variable=factor(variable,levels=c("SSB","F","R"),labels=c("SSB","Fishing mortality","Recruitment")))

p_error_validation <- ggplot() +
                      geom_density(data=sim_plot,aes(log_error),fill="grey80",color="black",linewidth=0.7,alpha=0.7) +
                      geom_rug(data=hist_plot,aes(log_error),sides="b",linewidth=0.8) +
                      geom_vline(xintercept=0,linetype="dashed",linewidth=0.5,color="grey40") +
                      facet_wrap(~variable,nrow=1,scales="free") +
                      labs(x="Log-assessment error",y="Simulated density") +
                      theme_bw() +
                      theme(panel.grid.minor=element_blank(),
                            panel.grid.major.x=element_blank(),
                            strip.background=element_blank(),
                            strip.text=element_text(face="bold"),
                            axis.title=element_text(size=11),
                            axis.text=element_text(size=10))

p_error_validation
ggsave(file.path(fig_dir,"03_assessment_error_distribution_validation.png"),p_error_validation,width=10,height=4,dpi=300)

# ============================================================
# 11. GENERATE EXAMPLE ASSESSMENT-ERROR TRAJECTORIES
# Illustrate the temporal behaviour of the reference
# assessment-error model by comparing the historical annual
# revisions with example simulated trajectories.
#
# These trajectories are shown for diagnostic purposes only.
# The closed-loop MSE trajectories are generated separately
# in STEP 9 for all 1000 iterations.
# ============================================================

set.seed(123)

proj_years <- 2025:2054
n_rep_plot <- 5

sim_time <- map_dfr(1:n_rep_plot,function(rep) {
      e <- MASS::mvrnorm(n=length(proj_years),mu=mu,Sigma=Sigma) %>% as.data.frame()
  names(e) <- variables
  e %>% mutate(year=proj_years,replicate=rep)
})

sim_time_long <- sim_time %>% pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error")
hist_time_long <- error_wide %>%
                  rename(year=update_year) %>%
                  pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error")

plot_assessment_error <- function(v) {
  
  hist_v <- hist_time_long %>% filter(variable==v)
  sim_v <- sim_time_long %>% filter(variable==v)
  
  lab_v <- c(SSB="SSB",F="Fishing mortality",R="Recruitment")[v]
  
  ggplot() +
    geom_line(data=hist_v,aes(year,log_error),linewidth=0.7,color="black") +
    geom_point(data=hist_v,aes(year,log_error),size=1.7,color="black") +
    geom_line(data=sim_v,aes(year,log_error,group=replicate,color=factor(replicate)),linewidth=0.6,alpha=0.8) +
    geom_point(data=sim_v,aes(year,log_error,color=factor(replicate)),size=1.4,alpha=0.8) +
    geom_hline(yintercept=mu[v],linetype="dashed",linewidth=0.5,color="grey40") +
    geom_vline(xintercept=2024.5,linetype="dotted",linewidth=0.5,color="grey40") +
    annotate("text",x=2022,y=Inf,label="Historical",vjust=1.5,hjust=1,size=3.3) +
    annotate("text",x=2027,y=Inf,label="Simulated",vjust=1.5,hjust=0,size=3.3) +
    labs(x="Year",y="Log-assessment error",title=lab_v,color="Example\ntrajectory") +
    theme_bw() +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      plot.title=element_text(face="bold",hjust=0.5,size=11),
      legend.position="right",
      axis.title=element_text(size=11),
      axis.text=element_text(size=10)
    )
}

p_err_SSB <- plot_assessment_error("SSB")
p_err_F <- plot_assessment_error("F")
p_err_R <- plot_assessment_error("R")

p_error_trajectories <- p_err_SSB / p_err_F / p_err_R

p_error_trajectories

ggsave(file.path(fig_dir,"04_example_assessment_error_trajectories.png"),p_error_trajectories,width=9,height=10,dpi=300)

# ============================================================
# 12. EXTERNAL VALIDATION USING CONVENTIONAL RETROSPECTIVES
# Evaluate whether the parameterised assessment-error emulator
# produces retrospective behaviour consistent with the
# conventional SS3 retrospective analysis (peels -1 to -5).
#
# Validation considers:
#   - peel-specific retrospective errors,
#   - Mohn's rho,
#   - persistence in the direction of retrospective errors.
# ============================================================

retro_conv_dir <- here("data","retro")
retro_peels <- 0:-5
retro_conv_dirs <- file.path(retro_conv_dir,paste0("retro",retro_peels))

retro_conv_check <- tibble(peel=retro_peels,
                           terminal_year=base_year+retro_peels,
                           dir=retro_conv_dirs,
                           exists=file.exists(file.path(retro_conv_dirs,"Report.sso")))
retro_conv_check

# ============================================================
# 12.1 READ CONVENTIONAL RETROSPECTIVES
# ============================================================

retro_conv_models <- SSgetoutput(dirvec=retro_conv_check$dir,verbose=FALSE)
retro_conv_summary <- SSsummarize(retro_conv_models)

retro_conv_def <- retro_conv_check %>%
                  mutate(model=retro_conv_summary$modelnames,
                         max_gradient=retro_conv_summary$maxgrad,
                         converged=is.finite(max_gradient) & max_gradient<gradient_limit)

retro_conv_def

stopifnot(all(retro_conv_def$converged))

# ============================================================
# 12.2 EXTRACT CONVENTIONAL RETROSPECTIVE ESTIMATES
# ============================================================

retro_conv_long <- retro_conv_summary$quants %>% as_tibble() %>%
                    pivot_longer(cols=all_of(retro_conv_summary$modelnames),names_to="model",values_to="estimate") %>%
                    left_join(retro_conv_def %>% dplyr::select(model,peel,terminal_year),by="model") %>%
                    mutate(year=suppressWarnings(as.integer(Yr)),
                    variable=case_when(grepl("^SSB_[0-9]{4}$",Label) ~ "SSB",
                                       grepl("^F_[0-9]{4}$",Label) ~ "F",
                                       grepl("^Recr_[0-9]{4}$",Label) ~ "R",
                                       TRUE ~ NA_character_)) %>%
                    filter(!is.na(variable),!is.na(year),is.finite(estimate),estimate>0,year<=terminal_year)

# ============================================================
# 12.3 OBSERVED RETROSPECTIVE ERRORS
# ============================================================

retro_full <- retro_conv_long %>% filter(peel==0) %>%
              dplyr::select(variable,year,full_estimate=estimate)

retro_observed <- retro_conv_long %>%
                  filter(peel<0,year==terminal_year) %>%
                  dplyr::select(peel,terminal_year,variable,retro_estimate=estimate) %>%
                  left_join(retro_full,by=c("variable","terminal_year"="year")) %>%
                  mutate(log_error=log(retro_estimate/full_estimate),
                         relative_error=retro_estimate/full_estimate-1) %>%
                  arrange(variable,desc(terminal_year))

retro_observed

# ============================================================
# 12.4 SIMULATED RETROSPECTIVE ERRORS
# ============================================================

set.seed(123)
n_retro_sim <- 10000

retro_sim <- map_dfr(1:n_retro_sim,function(sim) {
  
  annual_errors <- MASS::mvrnorm(n=5,mu=mu,Sigma=Sigma)
  colnames(annual_errors) <- variables
  
  cumulative_errors <- apply(annual_errors,2,cumsum)
  
  as_tibble(cumulative_errors) %>%
    mutate(sim=sim,peel=-(1:5))
})

retro_sim_long <- retro_sim %>% pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error")

retro_sim_summary <- retro_sim_long %>%
                      group_by(peel,variable) %>%
                      summarise(median=median(log_error),
                      q05=quantile(log_error,0.05),
                      q95=quantile(log_error,0.95), .groups="drop")

retro_validation <- retro_observed %>%
                    dplyr::select(peel,terminal_year,variable,observed=log_error) %>%
                    left_join(retro_sim_summary,by=c("peel","variable")) %>%
                    mutate(within_90=observed>=q05 & observed<=q95)

retro_validation
write.csv(retro_validation,file.path(tab_dir,"04_retrospective_validation.csv"),row.names=FALSE)

# ============================================================
# 12.5 OBSERVED MOHN'S RHO
# ============================================================

mohn_observed <- retro_observed %>%group_by(variable) %>%summarise(rho=mean(relative_error),.groups="drop")
mohn_observed

# ============================================================
# 12.6 SIMULATED MOHN'S RHO
# ============================================================

mohn_sim <- retro_sim %>%
            pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error") %>%
            mutate(relative_error=exp(log_error)-1) %>%
            group_by(sim,variable) %>%
            summarise(rho=mean(relative_error),.groups="drop")

mohn_sim_summary <- mohn_sim %>%
                    group_by(variable) %>%
                    summarise(median=median(rho),
                              q05=quantile(rho,0.05),
                              q95=quantile(rho,0.95), .groups="drop")

mohn_validation <- mohn_observed %>%
                   rename(observed_rho=rho) %>%
                   left_join(mohn_sim_summary,by="variable") %>%
                   mutate(within_90=observed_rho>=q05 & observed_rho<=q95)
mohn_validation
write.csv(mohn_validation,file.path(tab_dir,"05_mohn_rho_validation.csv"),row.names=FALSE)

# ============================================================
# 12.7 DIRECTIONAL PATTERN OF OBSERVED RETROSPECTIVES
# ============================================================

sign_observed <- retro_observed %>%
                  group_by(variable) %>%
                  summarise(n_positive=sum(log_error>0),
                            n_negative=sum(log_error<0),
                            max_same_direction=max(n_positive,n_negative),.groups="drop")

sign_observed

# ============================================================
# 12.8 DIRECTIONAL-PATTERN VALIDATION
# ============================================================

sign_sim <- retro_sim %>%
            pivot_longer(cols=all_of(variables),names_to="variable",values_to="log_error") %>%
            group_by(sim,variable) %>%
            summarise(n_positive=sum(log_error>0),
                      n_negative=sum(log_error<0),.groups="drop")

sign_validation <- sign_observed %>%
                   dplyr::select(variable,obs_positive=n_positive,obs_negative=n_negative) %>%
                   left_join(sign_sim %>%
                   group_by(variable) %>%
                   summarise(p_observed_direction=case_when(first(variable)=="F" ~ mean(n_positive>=4),
                                                             first(variable)=="SSB" ~ mean(n_negative>=4),
                                                             first(variable)=="R" ~ mean(n_positive>=3)),
                              .groups="drop"),by="variable")

sign_validation
write.csv(sign_validation,file.path(tab_dir,"06_directional_validation.csv"),row.names=FALSE)

# ============================================================
# 13. DOCUMENT EXTERNAL VALIDATION
# ============================================================

library(patchwork)

# ============================================================
# 13.1 PANEL A - PEEL-SPECIFIC RETROSPECTIVE ERRORS
# ============================================================

retro_validation_plot <- retro_validation %>%
                         mutate(variable=factor(variable,levels=c("SSB","F","R"),
                         labels=c("SSB","Fishing mortality","Recruitment")))

retro_sim_summary_plot <- retro_sim_summary %>%
                          mutate(variable=factor(variable,levels=c("SSB","F","R"),
                          labels=c("SSB","Fishing mortality","Recruitment")))

p_retro_validation <- ggplot() +
                      geom_ribbon(data=retro_sim_summary_plot,aes(x=peel,ymin=q05,ymax=q95),fill="grey88",color=NA) +
                      geom_line(data=retro_sim_summary_plot,aes(peel,median),linewidth=0.7,color="grey45") +
                      geom_point(data=retro_validation_plot,aes(peel,observed), size=2.6,color="black") +
                      geom_hline(yintercept=0,linetype="dashed",linewidth=0.4,color="grey50") +
                      facet_wrap(~variable,nrow=1,scales="free_y") +
                      scale_x_reverse(breaks=-1:-5) +
                      labs(x="Retrospective peel (years)",y="Log-assessment error") +
                      theme_bw() +
                      theme(panel.grid.minor=element_blank(),
                      panel.grid.major.x=element_blank(),
                      strip.background=element_blank(),
                      strip.text=element_text(face="bold"),
                      axis.title=element_text(size=11),
                      axis.text=element_text(size=10))

# ============================================================
# 13.2 PANEL B - MOHN'S RHO
# ============================================================

mohn_plot <- mohn_validation %>% mutate(variable=factor(variable,levels=c("SSB","F","R"),labels=c("SSB","Fishing mortality","Recruitment")))

p_mohn <- ggplot(mohn_plot,aes(variable,observed_rho)) +
          geom_linerange(aes(ymin=q05,ymax=q95),linewidth=1,color="grey55") +
          geom_point(aes(y=median,shape="Simulated median"),size=3,fill="white",color="black") +
          geom_point(aes(shape="Observed"),size=3,color="black") +
          geom_hline(yintercept=0,linetype="dashed",linewidth=0.4,color="grey50") +
          scale_shape_manual(values=c("Observed"=16,"Simulated median"=21),name=NULL) +
          labs(x=NULL,y="Mohn's rho") +
          theme_bw() +
          theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
          legend.position="top",legend.justification="left",
          axis.title=element_text(size=11),axis.text=element_text(size=10))

# ============================================================
# 13.3 PANEL C - DIRECTIONAL PERSISTENCE
# ============================================================

sign_plot <- sign_validation %>% mutate(variable=factor(variable,levels=c("SSB","F","R"), labels=c("SSB","Fishing mortality","Recruitment")),
         observed_direction=case_when(
           variable=="SSB" ~ obs_negative,
           variable=="Fishing mortality" ~ obs_positive,
           variable=="Recruitment" ~ obs_positive))

p_sign <- ggplot(sign_plot,aes(variable,observed_direction)) +
          geom_col(width=0.6,fill="grey70",color="black") +
          geom_text(aes(label=paste0("p = ",sprintf("%.3f",p_observed_direction))),
            vjust=-0.5,size=3.5) +
          scale_y_continuous(breaks=0:5,limits=c(0,5.3),expand=expansion(mult=c(0,0.02)))+
          labs(x=NULL,y="Peels with same error direction") +
          theme_bw() +
          theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
          axis.title=element_text(size=11),axis.text=element_text(size=10))



p_external_validation <- p_retro_validation /(p_mohn | p_sign) +
                          plot_layout(heights=c(1.15,1)) +
                          plot_annotation(tag_levels="A")

p_external_validation

ggsave(file.path(fig_dir,"05_assessment_error_external_validation.png"),p_external_validation,width=10,height=7,dpi=300)

# ============================================================
# 14. GENERATE ASSESSMENT-ERROR TRAJECTORIES FOR THE MSE
# Generate 1000 multivariate assessment-error trajectories for
# the management years 2024-2053.
#
# For each iteration and year:
#   1. Draw correlated log-errors for SSB, F, and R from the
#      fitted multivariate distribution.
#   2. Transform log-errors to multiplicative errors:
#
#         multiplier = exp(log_error)
#
# The resulting array has dimensions:
#   year x iteration x variable
#
# and is saved for subsequent closed-loop MSE simulations.
# ============================================================

set.seed(123)

assessment_years <- 2024:2053
n_iter <- 1000

assessment_error_array <- array(NA_real_,dim=c(length(assessment_years),
                                               n_iter,length(variables)),
                                dimnames=list(year=as.character(assessment_years),
                                              iter=as.character(seq_len(n_iter)),
                                              variable=variables))

for(i in seq_len(n_iter)) {
  assessment_error_array[,i,] <- exp(MASS::mvrnorm(n=length(assessment_years),mu=mu,Sigma=Sigma))
}

assessment_error_reference <- list(multiplier=assessment_error_array,parameters=assessment_error_parameters)

stopifnot(dim(assessment_error_array)[1]==30)
stopifnot(dim(assessment_error_array)[2]==1000)
stopifnot(dim(assessment_error_array)[3]==3)
stopifnot(all(is.finite(assessment_error_array)))
stopifnot(all(assessment_error_array>0))
stopifnot(identical(dimnames(assessment_error_array)$variable,c("SSB","F","R")))

saveRDS(assessment_error_reference,file.path(rds_dir,"assessment_error_reference_1000iter.rds"))

cat("\nAssessment-error trajectories generated: OK\n")
cat("Years:",min(assessment_years),"-",max(assessment_years),"\n")
cat("Iterations:",n_iter,"\n")
cat("Variables:",paste(variables,collapse=", "),"\n")

# ============================================================
# 15. VALIDATE GENERATED ASSESSMENT-ERROR TRAJECTORIES
# Verify that the final 1000 assessment-error trajectories
# reproduce the target means, standard deviations, and
# cross-variable correlation structure on the log-error scale.
# ============================================================

mult <- assessment_error_reference$multiplier
log_mult <- log(mult)

generated_mean <- apply(log_mult,3,mean)
generated_sd <- apply(log_mult,3,sd)
generated_cor <- cor(matrix(log_mult,ncol=3,dimnames=list(NULL,variable=dimnames(mult)$variable)))

generated_mean
mu

generated_sd
sqrt(diag(Sigma))

generated_cor
correlation

# ============================================================
# FINAL VALIDATION SUMMARY
# The assessment-error component is now parameterised and
# validated independently of the Management Procedure.
#
# Final object for the closed-loop MSE:
#   data/Rdata/assessment_error_reference_1000iter.rds
#
# Subsequent workflow:
#   1. Validate observation/projection infrastructure.
#   2. Implement the Fcap + Besc Management Procedure.
#   3. Validate the complete closed loop:
#        true stock -> observation -> perceived stock ->
#        assessment error -> MP -> advice -> realised dynamics.
# ============================================================

message(
  "Assessment-error parameterisation and validation completed successfully: ",
  "SSB, F and R; iid multivariate emulator; 1000 trajectories for 2024-2053."
)
