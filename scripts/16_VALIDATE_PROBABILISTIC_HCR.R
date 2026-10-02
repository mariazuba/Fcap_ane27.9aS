# ============================================================
# 16_VALIDATE_PROBABILISTIC_HCR.R
#
# Validation of probabilistic Fcap + Bescapement HCR
# Anchovy 9a South shortcut MSE
#
# Objective:
# Validate the probabilistic escapement constraint before its
# implementation in the FLBEIA closed-loop management procedure.
#
# The probabilistic HCR requires:
#
#   P(SSBforecast < Besc) <= alpha
#
# which under the normal approximation is equivalent to:
#
#   Qalpha(SSBforecast) >= Besc
#
# with:
#
#   SD(SSBforecast) = CV_SSB * mean(SSBforecast)
#
# Validation sequence:
# 1. Load forecast-uncertainty characterisation
# 2. Define probabilistic escapement functions
# 3. Validate probability calculation against SS3
# 4. Derive minimum mean SSB required for escapement
# 5. Define and validate the one-year projection engine
# 6. Load SS3 benchmark inputs
# 7. Validate F response and characterise informative Fcap range
# 8. Apply forecast-uncertainty scenarios
# 9. Validate the probabilistic Fcap + Besc HCR
# 10. Validate deterministic equivalence
# 11. Validate current-assessment uncertainty
# 12. Evaluate contemporary bootstrap uncertainty
# 13. Summarise scenario responses
# 14. Validate HCR decision branches across stock states
# 15. Produce methodological tables and figures
# 16. Save validation outputs
# 17. Print final validation summary
# ============================================================

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(r4ss)

# ============================================================ 
# 1. DIRECTORIES AND SETTINGS 
# ============================================================ 
# Defines input/output directories and the reference HCR settings 
# used throughout the validation: Fcap, Besc and the maximum 
# acceptable probability of falling below Besc (alpha = 0.05). 
# ============================================================
rds_dir <- here("data","Rdata")
out_dir <- here("outputs","MP")
fig_dir <- file.path(out_dir,"figures")
table_dir <- file.path(out_dir,"tables")

dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)

Fcap <- 1.5
Besc <- 6561
alpha <- 0.05

# ============================================================ 
# 1A. LOAD FORECAST-UNCERTAINTY CHARACTERISATION 
# ============================================================ 
# Loads the forecast-uncertainty object generated in Script 15. 
# The base SS3 forecast CV represents the current-assessment 
# uncertainty, while the 100 bootstrap CVs describe the empirical 
# distribution of contemporary forecast uncertainty. 
# ============================================================

uncertainty_file <- file.path(rds_dir,"forecast_uncertainty_characterisation.rds")

stopifnot(file.exists(uncertainty_file))

uncertainty_objects <- readRDS(uncertainty_file)

stopifnot(all(c("base","contemporary","scenarios") %in% names(uncertainty_objects)))

om_forecast_uncertainty <- uncertainty_objects$contemporary
base_SSB_CV <- uncertainty_objects$base$SSB_CV
bootstrap_SSB_CV <- om_forecast_uncertainty$SSB_2025_CV

stopifnot(length(base_SSB_CV)==1)
stopifnot(is.finite(base_SSB_CV))
stopifnot(base_SSB_CV>0)
stopifnot(nrow(om_forecast_uncertainty)==100)
stopifnot(all(is.finite(bootstrap_SSB_CV)))
stopifnot(all(bootstrap_SSB_CV>0))

cat("Base 2024 SSB forecast CV =",round(base_SSB_CV,4),"\n")
cat("Bootstrap median SSB forecast CV =",round(median(bootstrap_SSB_CV),4),"\n")

# ============================================================ 
# 1B. DEFINE UNCERTAINTY SCENARIOS 
# ============================================================ 
# Defines the three uncertainty scenarios used to test the HCR: 
# S0: deterministic forecast (CV = 0); 
# S1: current-assessment forecast uncertainty; 
# S2: contemporary uncertainty represented by the empirical 
# bootstrap distribution of forecast SSB CVs. 
# ============================================================

scenario_summary <- data.frame(scenario=c("S0","S1","S2"),
                               name=c("Deterministic","Current assessment","Contemporary uncertainty"),
                               CV_reference=c(0,base_SSB_CV,median(bootstrap_SSB_CV)),
                               implementation=c("CV = 0",
                                                "Base 2024 SS3 forecast CV",
                                                "CV sampled from contemporary bootstrap distribution"))

scenario_summary

# ============================================================ 
# 2. PROBABILISTIC ESCAPEMENT FUNCTIONS 
# ============================================================ 
# Defines the probabilistic escapement criterion. 
# 
# For a forecast mean SSB and its associated CV, the functions 
# calculate the forecast SD, the lower 5% SSB quantile and the 
# probability that forecast SSB falls below Besc. 
# 
# The HCR is considered acceptable when: 
# 
# P(SSBforecast < Besc) <= alpha 
# 
# which is equivalent to: 
# 
# Qalpha(SSBforecast) >= Besc 
# ============================================================
get_SSB_risk <- function(SSB,CV,Besc=6561,alpha=0.05) {
  if(CV==0) return(data.frame(SSB=SSB,CV=CV,SSB_SD=0,SSB_P05=SSB,p_below=as.numeric(SSB<Besc),acceptable=SSB>=Besc))
  SSB_SD <- CV*SSB
  SSB_P05 <- qnorm(alpha,mean=SSB,sd=SSB_SD)
  p_below <- pnorm(Besc,mean=SSB,sd=SSB_SD)
  data.frame(SSB=SSB,CV=CV,SSB_SD=SSB_SD,SSB_P05=SSB_P05,p_below=p_below,acceptable=p_below<=alpha)
}

get_SSB_required <- function(CV,Besc=6561,alpha=0.05) {
  if(CV==0) return(Besc)
  denominator <- 1+qnorm(alpha)*CV
  if(denominator<=0) return(Inf)
  Besc/denominator
}

# ============================================================ 
# 3. VALIDATE PROBABILITY CALCULATION AGAINST SS3 
# ============================================================ 
# Checks that the probabilistic emulator reproduces the probability 
# obtained directly from the SS3 first-year forecast mean and SD. 
# 
# This validates the translation of SS3 forecast uncertainty into 
# the probability calculation subsequently used by the shortcut HCR. 
# ============================================================

cat("\n============================================================\n")
cat("VALIDATION 1: SS3 PROBABILITY BENCHMARK\n")
cat("============================================================\n")

SSB_SS3 <- 13929.3
SSB_SD_SS3 <- 3336.83
SSB_CV_SS3 <- SSB_SD_SS3/SSB_SS3

SS3_probability_direct <- pnorm(Besc,mean=SSB_SS3,sd=SSB_SD_SS3)

SS3_probability_emulator <- get_SSB_risk(SSB=SSB_SS3,CV=SSB_CV_SS3,Besc=Besc,alpha=alpha)

Table16_1 <- data.frame(source=c("SS3 direct","Probabilistic emulator"),SSB_mean=c(SSB_SS3,SS3_probability_emulator$SSB),SSB_SD=c(SSB_SD_SS3,SS3_probability_emulator$SSB_SD),SSB_CV=c(SSB_CV_SS3,SS3_probability_emulator$CV),SSB_P05=c(qnorm(alpha,mean=SSB_SS3,sd=SSB_SD_SS3),SS3_probability_emulator$SSB_P05),p_below_Besc=c(SS3_probability_direct,SS3_probability_emulator$p_below))

Table16_1

stopifnot(abs(SS3_probability_direct-SS3_probability_emulator$p_below)<1e-12)
stopifnot(abs(SSB_SD_SS3-SS3_probability_emulator$SSB_SD)<1e-6)
stopifnot((SS3_probability_emulator$p_below<=alpha)==(SS3_probability_emulator$SSB_P05>=Besc))

cat("PASS: probability calculation reproduces SS3 benchmark\n")
cat("P(SSB < Besc) =",round(SS3_probability_emulator$p_below,6),"\n")

# ============================================================ 
# 4. MINIMUM MEAN SSB REQUIRED FOR ESCAPEMENT 
# ============================================================ 
# Calculates the minimum forecast mean SSB required to satisfy the 
# 5% escapement-risk criterion for different levels of uncertainty. 
# 
# This demonstrates an important property of the probabilistic HCR: 
# increasing forecast uncertainty requires a progressively larger 
# expected SSB to maintain the same probability of exceeding Besc. 
# ============================================================

cat("\n============================================================\n")
cat("VALIDATION 2: REQUIRED MEAN SSB\n")
cat("============================================================\n")

CV_values <- c(S0=0,S1=base_SSB_CV,S2_P05=as.numeric(quantile(bootstrap_SSB_CV,0.05)),S2_P50=as.numeric(quantile(bootstrap_SSB_CV,0.50)),S2_P95=as.numeric(quantile(bootstrap_SSB_CV,0.95)))

Table16_2 <- data.frame(case=names(CV_values),CV=as.numeric(CV_values),SSB_required=sapply(CV_values,get_SSB_required,Besc=Besc,alpha=alpha))

Table16_2$buffer_above_Besc <- Table16_2$SSB_required-Besc
Table16_2$relative_buffer <- Table16_2$SSB_required/Besc
Table16_2$p_at_required <- mapply(function(SSB,CV) get_SSB_risk(SSB=SSB,CV=CV,Besc=Besc,alpha=alpha)$p_below,Table16_2$SSB_required,Table16_2$CV)

Table16_2

Table16_2_ordered <- Table16_2 %>% arrange(CV)

stopifnot(Table16_2$SSB_required[Table16_2$case=="S0"]==Besc)
stopifnot(all(diff(Table16_2_ordered$SSB_required)>=0))
stopifnot(abs(Table16_2$p_at_required[Table16_2$case=="S0"]-0)<1e-12)
stopifnot(all(abs(Table16_2$p_at_required[Table16_2$case!="S0"]-alpha)<1e-10))

cat("PASS: required mean SSB increases with forecast uncertainty\n")
cat("PASS: SSB_required corresponds exactly to the 5% risk boundary\n")

# ============================================================ 
# 5. SS3-LIKE ONE-YEAR PROJECTION ENGINE 
# ============================================================ 
# Defines the one-year seasonal projection used by the shortcut HCR 
# to evaluate candidate fishing mortality values. 
# 
# Population abundance is projected sequentially through four 
# quarters. Fishing mortality is distributed among quarters using 
# the observed seasonal F proportions, recruitment enters in Q3, 
# and spawning-stock biomass is evaluated at the start of Q2. 
# 
# Natural mortality is supplied as an annual instantaneous rate and 
# converted to quarterly scale (M * 0.25). Fishing mortality is 
# already allocated among quarters through propf and is therefore 
# not multiplied again by dt. 
# ============================================================

project_ss3_one_year <- function(N_Q1,M,Frel,W_spawn,W_catch,Fbase,propf,a,b,dt=0.25,M_is_seasonal=FALSE) {
  ages <- names(N_Q1)
  seasons <- 1:4
  Fseason <- Fbase*propf
  N <- matrix(0,nrow=length(ages),ncol=4,dimnames=list(age=ages,season=seasons))
  Fage <- matrix(0,nrow=length(ages),ncol=4,dimnames=list(age=ages,season=seasons))
  Z <- matrix(0,nrow=length(ages),ncol=4,dimnames=list(age=ages,season=seasons))
  catch_num <- matrix(0,nrow=length(ages),ncol=4,dimnames=list(age=ages,season=seasons))
  catch_bio <- numeric(4)
  N[,1] <- N_Q1
  
  for(s in seasons) {
    Ms <- if(M_is_seasonal) M[,s] else M*dt
    Fs <- Fseason[s]*Frel[,s]
    Zs <- Ms+Fs
    Fage[,s] <- Fs
    Z[,s] <- Zs
    catch_num[,s] <- N[,s]*(Fs/Zs)*(1-exp(-Zs))
    catch_bio[s] <- sum(catch_num[,s]*W_catch[,s])
    if(s==1) SSB_Q2 <- sum(N[,1]*exp(-Zs)*W_spawn)
    if(s<4) N[,s+1] <- N[,s]*exp(-Zs)
    if(s==2) {
      R_Q3 <- a*SSB_Q2/(b+SSB_Q2)
      N["0",3] <- R_Q3
    }
  }
  
  N_endQ4 <- N[,4]*exp(-Z[,4])
  N_next_Q1 <- c(`0`=0,`1`=unname(N_endQ4["0"]),`2`=unname(N_endQ4["1"]),`3`=unname(N_endQ4["2"]+N_endQ4["3"]))
  
  list(Fbase=Fbase,Fseason=Fseason,N=N,Fage=Fage,Z=Z,SSB_Q2=SSB_Q2,R_Q3=R_Q3,catch_num=catch_num,catch_bio=catch_bio,TAC=sum(catch_bio),N_endQ4=N_endQ4,N_next_Q1=N_next_Q1)
}

# ============================================================ 
# 6. LOAD SS3 BENCHMARK INPUTS 
# ============================================================ 
# Extracts the population numbers-at-age, seasonal selectivity, 
# catch weights and spawning weights required to reproduce the 
# first forecast year from the reference SS3 assessment. 
# 
# These inputs provide an independent benchmark against which the 
# shortcut one-year projection is validated. 
# ============================================================
forecast_dir <- here("data","ss3_forecast_validation")

ss3_fore <- r4ss::SS_output(dir=forecast_dir,covar=FALSE,verbose=FALSE)

ss3_N <- ss3_fore$natage %>% filter(Yr==2025,`Beg/Mid`=="B")
ss3_sel <- ss3_fore$ageselex %>% filter(Yr==2025,Factor=="Asel")

N_Q1 <- unlist(ss3_N %>% filter(Seas==1) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)
names(N_Q1) <- 0:3

sel_Q1 <- unlist(ss3_sel %>% filter(Fleet==1) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)
sel_Q2 <- unlist(ss3_sel %>% filter(Fleet==2) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)
sel_Q3 <- unlist(ss3_sel %>% filter(Fleet==3) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)
sel_Q4 <- unlist(ss3_sel %>% filter(Fleet==4) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)

sel <- cbind(Q1=sel_Q1,Q2=sel_Q2,Q3=sel_Q3,Q4=sel_Q4)
rownames(sel) <- 0:3

W_catch <- matrix(NA_real_,nrow=4,ncol=4,dimnames=list(age=0:3,season=1:4))

for(s in 1:4) W_catch[,s] <- unlist(ss3_fore$wtatage %>% filter(year==2025,seas==s,fleet==s) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)

W_spawn <- unlist(ss3_fore$wtatage %>% filter(year==2025,seas==2,fleet==-2) %>% dplyr::select(`0`,`1`,`2`,`3`),use.names=FALSE)
names(W_spawn) <- 0:3

M <- c(`0`=2.97,`1`=1.60,`2`=2.48,`3`=2.48)
propf <- c(0.138,0.417,0.334,0.111)

a_SS3 <- 4*0.8*8292080/(5*0.8-1)
b_SS3 <- 17217.7*(1-0.8)/(5*0.8-1)

# ============================================================ 
# 7. VALIDATE F RESPONSE 
# ============================================================ 
# Evaluates the one-year projection over the reference F range. 
# 
# The basic expected response is checked before applying uncertainty: 
# increasing F should decrease forecast SSB and increase catch over 
# the F range relevant to the candidate management procedures. 
# ============================================================

cat("\n============================================================\n")
cat("VALIDATION 3: F RESPONSE\n")
cat("============================================================\n")

F_grid <- seq(0,Fcap,by=0.01)

project_F <- function(F) {
  p <- project_ss3_one_year(N_Q1=N_Q1,M=M,Frel=sel,W_spawn=W_spawn,W_catch=W_catch,Fbase=F,propf=propf,a=a_SS3,b=b_SS3,dt=0.25,M_is_seasonal=FALSE)
  data.frame(F=F,SSB=p$SSB_Q2,TAC=p$TAC)
}

F_response <- bind_rows(lapply(F_grid,project_F))

stopifnot(all(diff(F_response$SSB)<=0))
stopifnot(all(diff(F_response$TAC)>=0))

cat("PASS: SSB decreases monotonically with F\n")
cat("PASS: catch increases monotonically with F over tested range\n")

# ============================================================ 
# 7A. VALIDATE PROJECTION AGAINST SS3 FORECAST 
# ============================================================ 
# Compares the shortcut one-year projection with an independent SS3 
# first-year forecast benchmark at F = 1.112. 
# 
# Agreement in predicted catch provides a direct validation of the 
# seasonal mortality scaling and catch calculation used by the 
# shortcut projection engine. 
# ============================================================

proj_F1112 <- project_ss3_one_year(N_Q1=N_Q1,M=M,Frel=sel,W_spawn=W_spawn,W_catch=W_catch,Fbase=1.112,propf=propf,a=a_SS3,b=b_SS3,dt=0.25,M_is_seasonal=FALSE)

SS3_projection_validation <- data.frame(F=1.112,SS3_catch=7309,projected_catch=proj_F1112$TAC,difference=proj_F1112$TAC-7309,relative_difference=(proj_F1112$TAC-7309)/7309*100)

SS3_projection_validation

projection_tolerance_pct <- 0.1

stopifnot(abs(SS3_projection_validation$relative_difference) <= projection_tolerance_pct)

cat("PASS: one-year projection reproduces SS3 first-year catch benchmark", "within", projection_tolerance_pct, "% tolerance\n")


# ============================================================ 
# 7B. CHARACTERISE WITHIN-YEAR RESPONSE AT HIGH F 
# ============================================================ 
# Examines how increasing annual F affects seasonal fishing 
# mortality, vulnerable biomass, catch and exploitation. 
# 
# This diagnostic identifies within-year depletion and seasonal 
# saturation effects that may occur at high F and helps determine 
# whether very high Fcap values provide useful additional contrast 
# when defining candidate management procedures. 
# ============================================================

diagnose_high_F <- function(Fbase) {
  p <- project_ss3_one_year(N_Q1=N_Q1,M=M,Frel=sel,W_spawn=W_spawn,W_catch=W_catch,Fbase=Fbase,propf=propf,a=a_SS3,b=b_SS3,dt=0.25,M_is_seasonal=FALSE)
  bind_rows(lapply(1:4,function(s) {
    biomass_start <- sum(p$N[,s]*W_catch[,s])
    vulnerable_biomass <- sum(p$N[,s]*W_catch[,s]*sel[,s])
    data.frame(Fbase=Fbase,season=paste0("Q",s),Fseason=p$Fseason[s],biomass_start=biomass_start,vulnerable_biomass=vulnerable_biomass,catch=p$catch_bio[s],exploitation=ifelse(vulnerable_biomass>0,p$catch_bio[s]/vulnerable_biomass,NA_real_))
  }))
}

F_test <- c(0.5,1.0,1.112,1.5,2.0,2.139,2.5,3.0)

high_F_diag <- bind_rows(lapply(F_test,diagnose_high_F))

high_F_diag

# ============================================================ 
# 7C. CHARACTERISE INFORMATIVE FCAP RANGE 
# ============================================================ 
# Explores a broad F range before defining the candidate Fcap grid. 
# 
# For each F, the diagnostic records forecast SSB, seasonal and 
# annual catch, marginal catch gain and SSB loss. This allows the 
# informative Fcap range to be identified from the stock-fishery 
# response rather than selected arbitrarily. 
# ============================================================

F_grid_range <- seq(0.25,3.0,by=0.25)

project_F_range <- function(F) {
  p <- project_ss3_one_year(N_Q1=N_Q1,M=M,Frel=sel,W_spawn=W_spawn,W_catch=W_catch,Fbase=F,propf=propf,a=a_SS3,b=b_SS3,dt=0.25,M_is_seasonal=FALSE)
  data.frame(F=F,SSB=p$SSB_Q2,Catch_Q1=p$catch_bio[1],Catch_Q2=p$catch_bio[2],Catch_Q3=p$catch_bio[3],Catch_Q4=p$catch_bio[4],Catch_total=p$TAC)
}

Fcap_range_diag <- bind_rows(lapply(F_grid_range,project_F_range))

Fcap_range_diag <- Fcap_range_diag %>% mutate(delta_F=F-lag(F),delta_catch=Catch_total-lag(Catch_total),marginal_catch=delta_catch/delta_F,catch_gain_pct=100*delta_catch/lag(Catch_total),SSB_loss=lag(SSB)-SSB,SSB_loss_pct=100*SSB_loss/lag(SSB))

Fcap_range_diag

# ============================================================ 
# 8. APPLY UNCERTAINTY SCENARIOS TO F RESPONSE 
# ============================================================ 
# Applies S0, S1 and the median S2 forecast uncertainty to the same 
# biological F-response curve. 
# 
# The underlying population projection is unchanged; uncertainty 
# modifies only the perceived risk associated with each forecast 
# SSB and therefore the management decision. 
# ============================================================
evaluate_scenario <- function(F_response,CV,scenario) {
  risk <- bind_rows(lapply(F_response$SSB,get_SSB_risk,CV=CV,Besc=Besc,alpha=alpha))
  bind_cols(F_response,risk %>% select(CV,SSB_SD,SSB_P05,p_below,acceptable)) %>% mutate(scenario=scenario)
}

risk_S0 <- evaluate_scenario(F_response,0,"S0")
risk_S1 <- evaluate_scenario(F_response,base_SSB_CV,"S1")
risk_S2_P50 <- evaluate_scenario(F_response,median(bootstrap_SSB_CV),"S2_P50")

risk_response <- bind_rows(risk_S0,risk_S1,risk_S2_P50)

# ============================================================ 
# 8A. CHECK F RESPONSE UNDER UNCERTAINTY 
# ============================================================ 
# Verifies that the probabilistic quantities behave consistently 
# with increasing fishing pressure. 
# 
# As F increases, forecast SSB and its lower 5% quantile should 
# decrease, while the probability of falling below Besc should 
# increase. 
# ============================================================
Table16_Fresponse <- risk_response %>% group_by(scenario) %>% 
                      summarise(CV=first(CV),
                                SSB_F0=first(SSB),
                                SSB_Fcap=last(SSB),
                                P05_F0=first(SSB_P05),
                                P05_Fcap=last(SSB_P05),
                                pbelow_F0=first(p_below),
                                pbelow_Fcap=last(p_below),
                                Fcap_acceptable=last(acceptable),.groups="drop")

Table16_Fresponse

stopifnot(all(diff(risk_S0$SSB)<=0))
stopifnot(all(diff(risk_S1$SSB_P05)<=0))
stopifnot(all(diff(risk_S2_P50$SSB_P05)<=0))
stopifnot(all(diff(risk_S1$p_below)>=0))
stopifnot(all(diff(risk_S2_P50$p_below)>=0))

cat("PASS: SSB decreases with F\n")
cat("PASS: lower 5% SSB quantile decreases with F\n")
cat("PASS: escapement risk increases with F\n")

# ============================================================ 
# 8B. IDENTIFY APPROXIMATE 5% RISK BOUNDARY 
# ============================================================ 
# Identifies the largest F on the evaluation grid that satisfies 
# the escapement criterion for each uncertainty scenario. 
# 
# This provides a first numerical illustration of how increasing 
# forecast uncertainty can reduce the fishing mortality compatible 
# with the same Besc and risk tolerance. 
# ============================================================
risk_boundary <- risk_response %>% group_by(scenario) %>% 
                filter(acceptable) %>% slice_max(F,n=1,with_ties=FALSE) %>% 
                select(scenario,CV,F,SSB,SSB_P05,p_below,TAC) %>% ungroup()

SSB_required_check <- risk_boundary %>% 
                      mutate(SSB_required=sapply(CV,get_SSB_required,
                                                 Besc=Besc,
                                                 alpha=alpha),
                             difference=SSB-SSB_required)

risk_boundary

SSB_required_check

# ============================================================ 
# 8C. CHECK FCAP ACROSS CONTEMPORARY BOOTSTRAP CVS 
# ============================================================ 
# Evaluates whether Fcap satisfies the escapement criterion for each 
# of the 100 contemporary bootstrap forecast CVs. 
# 
# Analytical CV thresholds are also calculated for the transition 
# from Fcap to reduced F and from reduced F to complete closure. 
# ============================================================

proj_Fcap <- project_ss3_one_year(N_Q1=N_Q1,
                                  M=M,
                                  Frel=sel,
                                  W_spawn=W_spawn,
                                  W_catch=W_catch,
                                  Fbase=Fcap,
                                  propf=propf,
                                  a=a_SS3,
                                  b=b_SS3,
                                  dt=0.25,
                                  M_is_seasonal=FALSE)

risk_Fcap_bootstrap <- bind_rows(lapply(seq_along(bootstrap_SSB_CV),function(i) {
  
                        risk <- get_SSB_risk(SSB=proj_Fcap$SSB_Q2,
                                             CV=bootstrap_SSB_CV[i],
                                             Besc=Besc,
                                             alpha=alpha)
  data.frame(draw=i,
             CV=bootstrap_SSB_CV[i],
             SSB=proj_Fcap$SSB_Q2,
             SSB_P05=risk$SSB_P05,
             p_below=risk$p_below,
             acceptable=risk$acceptable)
}))

CV_max_Fcap <- (Besc/proj_Fcap$SSB_Q2-1)/qnorm(alpha)

proj_F0 <- project_ss3_one_year(N_Q1=N_Q1,
                                M=M,
                                Frel=sel,
                                W_spawn=W_spawn,
                                W_catch=W_catch,
                                Fbase=0,
                                propf=propf,
                                a=a_SS3,
                                b=b_SS3,
                                dt=0.25,
                                M_is_seasonal=FALSE)

CV_max_F0 <- (Besc/proj_F0$SSB_Q2-1)/qnorm(alpha)

stopifnot(all(risk_Fcap_bootstrap$acceptable==(bootstrap_SSB_CV<=CV_max_Fcap)))

cat("PASS: numerical and analytical Fcap acceptance are identical\n")

# ============================================================ 
# 8D. ANALYTICAL CLASSIFICATION OF HCR BRANCHES 
# ============================================================ 
# Uses the analytical CV thresholds to classify each bootstrap CV 
# into the expected HCR decision branch: 
# 
# Fcap : Fcap already satisfies the risk criterion 
# Reduced F : a value between 0 and Fcap is required 
# Closure : the criterion cannot be satisfied even at F = 0 
# 
# These classifications provide an independent benchmark for the 
# numerical HCR implementation. 
# ============================================================

bootstrap_branch_expected <- data.frame(draw=seq_along(bootstrap_SSB_CV),CV=bootstrap_SSB_CV)

bootstrap_branch_expected$branch_expected <- ifelse(bootstrap_branch_expected$CV<=CV_max_Fcap,"Fcap",
                                                    ifelse(bootstrap_branch_expected$CV<=CV_max_F0,"Reduced F","Closure"))

table(bootstrap_branch_expected$branch_expected)

# ============================================================ 
# 9. PROBABILISTIC FCAP + BESC HCR 
# ============================================================ 
# Implements the probabilistic management rule. 
# 
# The HCR first evaluates Fcap. If Fcap satisfies the escapement 
# criterion it is retained. Otherwise, if the criterion can still 
# be satisfied at F = 0, root finding identifies the maximum F for 
# which Q05(SSBforecast) = Besc. If the criterion cannot be met even 
# at F = 0, the fishery is closed. 
# ============================================================

get_F_HCR_prob <- function(N_Q1,M,Frel,W_spawn,W_catch,Fcap,Besc,propf,a,b,CV,alpha=0.05,dt=0.25,M_is_seasonal=FALSE) {
  run_F <- function(F) project_ss3_one_year(N_Q1=N_Q1,
                                            M=M,
                                            Frel=Frel,
                                            W_spawn=W_spawn,
                                            W_catch=W_catch,
                                            Fbase=F,
                                            propf=propf,
                                            a=a,
                                            b=b,
                                            dt=dt,
                                            M_is_seasonal=M_is_seasonal)
  
  get_risk <- function(proj) get_SSB_risk(SSB=proj$SSB_Q2,
                                          CV=CV,
                                          Besc=Besc,
                                          alpha=alpha)
  p0 <- run_F(0)
  pcap <- run_F(Fcap)
  risk0 <- get_risk(p0)
  riskcap <- get_risk(pcap)
  
  if(riskcap$acceptable) {
    Fstar <- Fcap
    proj <- pcap
    branch <- "Fcap"
  } else if(risk0$acceptable) {
    objective <- function(F) get_risk(run_F(F))$SSB_P05-Besc
    Fstar <- uniroot(objective,c(0,Fcap),tol=1e-10)$root
    proj <- run_F(Fstar)
    branch <- "Reduced F"
  } else {
    Fstar <- 0
    proj <- p0
    branch <- "Closure"
  }
  
  risk <- get_risk(proj)
  
  list(F=Fstar,
       TAC=proj$TAC,
       SSB=proj$SSB_Q2,
       SSB_SD=risk$SSB_SD,
       SSB_P05=risk$SSB_P05,
       p_below=risk$p_below,
       CV=CV,
       branch=branch,
       projection=proj)
}

# ============================================================ 
# 9A. RUN HCR ACROSS 100 BOOTSTRAP CVS 
# ============================================================ 
# Applies the probabilistic HCR to the complete empirical 
# distribution of contemporary forecast uncertainty. 
# 
# Biological conditions are held constant so that differences in 
# F advice arise only from differences in forecast uncertainty. 
# ============================================================

HCR_bootstrap <- bind_rows(lapply(seq_along(bootstrap_SSB_CV),function(i) {
  h <- get_F_HCR_prob(N_Q1=N_Q1,
                      M=M,
                      Frel=sel,
                      W_spawn=W_spawn,
                      W_catch=W_catch,
                      Fcap=Fcap,
                      Besc=Besc,
                      propf=propf,
                      a=a_SS3,
                      b=b_SS3,
                      CV=bootstrap_SSB_CV[i],
                      alpha=alpha,
                      dt=0.25,
                      M_is_seasonal=FALSE)
  data.frame(draw=i,
             CV=bootstrap_SSB_CV[i],
             F=h$F,
             TAC=h$TAC,
             SSB=h$SSB,
             SSB_SD=h$SSB_SD,
             SSB_P05=h$SSB_P05,
             p_below=h$p_below,
             branch=h$branch)
}))

# ============================================================ 
# 9B. VALIDATE HCR BRANCHES AGAINST ANALYTICAL BOUNDARIES 
# ============================================================ 
# Compares the numerical branch selected by the HCR with the branch 
# predicted independently from the analytical CV thresholds. 
# 
# Exact agreement verifies that the decision logic has been 
# implemented consistently. 
# ============================================================

HCR_validation <- HCR_bootstrap %>% left_join(bootstrap_branch_expected,by=c("draw","CV"))

stopifnot(all(HCR_validation$branch==HCR_validation$branch_expected))

cat("PASS: HCR branches reproduce analytical CV boundaries exactly\n")

# ============================================================ 
# 9C. VALIDATE REDUCED-F SOLUTIONS 
# ============================================================ 
# Checks the intermediate branch of the HCR. 
# 
# Every reduced-F solution must lie strictly between 0 and Fcap and 
# satisfy the risk boundary numerically: # # Q05(SSBforecast) = Besc 
# P(SSBforecast < Besc) = 0.05 
# ============================================================
HCR_reduced <- HCR_validation %>% filter(branch=="Reduced F") %>% arrange(CV)

stopifnot(nrow(HCR_reduced)>0)
stopifnot(all(HCR_reduced$F>0))
stopifnot(all(HCR_reduced$F<Fcap))
stopifnot(all(abs(HCR_reduced$SSB_P05-Besc)<1e-6))
stopifnot(all(abs(HCR_reduced$p_below-alpha)<1e-8))

cat("PASS: reduced-F solutions lie strictly between 0 and Fcap\n")
cat("PASS: reduced-F solutions satisfy P05 = Besc\n")
cat("PASS: reduced-F solutions satisfy P(SSB < Besc) = 0.05\n")

# ============================================================ 
# 9D. VALIDATE FCAP AND CLOSURE SOLUTIONS 
# ============================================================ 
# Checks the two limiting HCR branches. 
# 
# Fcap must only be retained when the risk criterion is satisfied, 
# while closure must occur only when the criterion cannot be 
# satisfied even in the absence of fishing. 
# ============================================================

HCR_Fcap <- HCR_validation %>% filter(branch=="Fcap")
HCR_closure <- HCR_validation %>% filter(branch=="Closure")

stopifnot(nrow(HCR_Fcap)>0)
stopifnot(nrow(HCR_closure)>0)
stopifnot(all(HCR_Fcap$F==Fcap))
stopifnot(all(HCR_Fcap$p_below<=alpha))
stopifnot(all(HCR_closure$F==0))
stopifnot(all(HCR_closure$TAC==0))
stopifnot(all(HCR_closure$p_below>alpha))

cat("PASS: Fcap branch satisfies the risk criterion at Fcap\n")
cat("PASS: closure occurs only when the criterion cannot be satisfied at F = 0\n")

# ============================================================ 
# 10. VALIDATE DETERMINISTIC EQUIVALENCE 
# ============================================================ 
# Tests the limiting case with no forecast uncertainty. 
# 
# When CV = 0, the probabilistic HCR must reproduce exactly the 
# original deterministic Fcap + Besc rule. This confirms that the 
# probabilistic formulation is an extension of, rather than a 
# different implementation from, the deterministic rule. 
# ============================================================

cat("\n============================================================\n")
cat("VALIDATION 4: DETERMINISTIC EQUIVALENCE\n")
cat("============================================================\n")

get_F_HCR_det <- function(N_Q1,M,Frel,W_spawn,W_catch,Fcap,Besc,propf,a,b,dt=0.25,M_is_seasonal=FALSE) {
  run_F <- function(F) project_ss3_one_year(N_Q1=N_Q1,
                                            M=M,
                                            Frel=Frel,
                                            W_spawn=W_spawn,
                                            W_catch=W_catch,
                                            Fbase=F,
                                            propf=propf,
                                            a=a,
                                            b=b,
                                            dt=dt,
                                            M_is_seasonal=M_is_seasonal)
  p0 <- run_F(0)
  pcap <- run_F(Fcap)
  
  if(pcap$SSB_Q2>=Besc) {
    Fstar <- Fcap
    proj <- pcap
  } else if(p0$SSB_Q2>=Besc) {
    Fstar <- uniroot(function(F) run_F(F)$SSB_Q2-Besc,c(0,Fcap),tol=1e-10)$root
    proj <- run_F(Fstar)
  } else {
    Fstar <- 0
    proj <- p0
  }
  
  list(F=Fstar,TAC=proj$TAC,SSB=proj$SSB_Q2)
}

HCR_det <- get_F_HCR_det(N_Q1=N_Q1,
                         M=M,
                         Frel=sel,
                         W_spawn=W_spawn,
                         W_catch=W_catch,
                         Fcap=Fcap,
                         Besc=Besc,
                         propf=propf,
                         a=a_SS3,
                         b=b_SS3)

HCR_S0 <- get_F_HCR_prob(N_Q1=N_Q1,
                         M=M,
                         Frel=sel,
                         W_spawn=W_spawn,
                         W_catch=W_catch,
                         Fcap=Fcap,
                         Besc=Besc,
                         propf=propf,
                         a=a_SS3,
                         b=b_SS3,
                         CV=0)

stopifnot(abs(HCR_det$F-HCR_S0$F)<1e-8)
stopifnot(abs(HCR_det$TAC-HCR_S0$TAC)<1e-6)
stopifnot(abs(HCR_det$SSB-HCR_S0$SSB)<1e-6)

cat("PASS: CV = 0 reproduces deterministic HCR exactly\n")

# ============================================================ 
# 11. VALIDATE CURRENT-ASSESSMENT SCENARIO 
# ============================================================ 
# Applies the HCR using the forecast CV estimated from the current 
# SS3 assessment (S1). 
# 
# This represents the uncertainty level associated with the current 
# assessment and verifies that the selected advice satisfies the 
# specified 5% escapement-risk criterion. 
# ============================================================
cat("\n============================================================\n")
cat("VALIDATION 5: CURRENT-ASSESSMENT SCENARIO\n")
cat("============================================================\n")

HCR_S1 <- get_F_HCR_prob(N_Q1=N_Q1,
                         M=M,
                         Frel=sel,
                         W_spawn=W_spawn,
                         W_catch=W_catch,
                         Fcap=Fcap,
                         Besc=Besc,
                         propf=propf,
                         a=a_SS3,
                         b=b_SS3,
                         CV=base_SSB_CV)

stopifnot(HCR_S1$p_below<=alpha+1e-8)

cat("PASS: S1 satisfies probabilistic escapement constraint\n")

# ============================================================ 
# 12. EVALUATE CONTEMPORARY UNCERTAINTY SCENARIO 
# ============================================================ 
# Evaluates the HCR across the 100 contemporary bootstrap forecast 
# CVs (S2). 
# 
# This characterises how variation in plausible assessment precision 
# alone changes F advice, TAC and the HCR decision branch. 
# ============================================================

cat("\n============================================================\n")
cat("VALIDATION 6: CONTEMPORARY UNCERTAINTY\n")
cat("============================================================\n")

HCR_S2 <- HCR_bootstrap

stopifnot(all(HCR_S2$p_below<=alpha+1e-8 | HCR_S2$branch=="Closure"))

CV_S2_P50 <- median(bootstrap_SSB_CV)

HCR_S2_P50 <- get_F_HCR_prob(N_Q1=N_Q1,
                             M=M,
                             Frel=sel,
                             W_spawn=W_spawn,
                             W_catch=W_catch,
                             Fcap=Fcap,
                             Besc=Besc,
                             propf=propf,
                             a=a_SS3,
                             b=b_SS3,
                             CV=CV_S2_P50,
                             alpha=alpha)

# ============================================================ 
# 13. SUMMARY OF SCENARIO RESPONSES 
# ============================================================ 
# Summarises the management response under S0, S1 and S2-P50. 
# 
# S2-P50 is evaluated explicitly using the median bootstrap CV, 
# rather than constructed from marginal medians of different HCR 
# runs, ensuring that each row represents one internally consistent 
# management decision. 
# ============================================================
Table16_3 <- bind_rows(
  data.frame(scenario="S0",
             CV=0,
             F=HCR_S0$F,
             TAC=HCR_S0$TAC,
             SSB=HCR_S0$SSB,
             SSB_P05=HCR_S0$SSB_P05,
             p_below=HCR_S0$p_below,
             branch=HCR_S0$branch),
  data.frame(scenario="S1",
             CV=base_SSB_CV,
             F=HCR_S1$F,
             TAC=HCR_S1$TAC,
             SSB=HCR_S1$SSB,
             SSB_P05=HCR_S1$SSB_P05,
             p_below=HCR_S1$p_below,
             branch=HCR_S1$branch),
  data.frame(scenario="S2_P50",
             CV=CV_S2_P50,
             F=HCR_S2_P50$F,
             TAC=HCR_S2_P50$TAC,
             SSB=HCR_S2_P50$SSB,
             SSB_P05=HCR_S2_P50$SSB_P05,
             p_below=HCR_S2_P50$p_below,
             branch=HCR_S2_P50$branch))

Table16_3

# ============================================================ 
# 14. VALIDATE HCR BRANCH LOGIC ACROSS STOCK STATES 
# ============================================================ 
# Tests whether all three HCR branches can be generated across a 
# broad range of stock states. 
# 
# Rather than imposing arbitrary abundance multipliers, the script 
# scans a range of stock scalings, identifies the regions producing 
# Fcap, reduced-F and closure decisions, and then extracts one 
# representative example of each branch. 
# ============================================================
cat("\n============================================================\n")
cat("VALIDATION 7: HCR BRANCH LOGIC\n")
cat("============================================================\n")

test_branch <- function(scale_N,CV=base_SSB_CV) {
  N_test <- N_Q1*scale_N
  get_F_HCR_prob(N_Q1=N_test,
                 M=M,
                 Frel=sel,
                 W_spawn=W_spawn,
                 W_catch=W_catch,
                 Fcap=Fcap,
                 Besc=Besc,
                 propf=propf,
                 a=a_SS3,
                 b=b_SS3,
                 CV=CV,
                 alpha=alpha)
}

scale_grid <- seq(0.1,2.0,by=0.01)

branch_scan <- bind_rows(lapply(scale_grid,function(scale_N) {
  h <- test_branch(scale_N)
  data.frame(scale_N=scale_N,
             F=h$F,
             TAC=h$TAC,
             SSB=h$SSB,
             SSB_P05=h$SSB_P05,
             p_below=h$p_below,
             branch=h$branch)
}))

branch_counts <- branch_scan %>% count(branch)

branch_counts

stopifnot(all(c("Fcap","Reduced F","Closure") %in% branch_scan$branch))

branch_high <- branch_scan %>% filter(branch=="Fcap") %>% slice(1)
branch_mid <- branch_scan %>% filter(branch=="Reduced F") %>% slice(round(dplyr::n()/2))
branch_low <- branch_scan %>% filter(branch=="Closure") %>% slice_tail(n=1)

Table16_4 <- bind_rows(branch_high %>% 
                         mutate(stock_state="High"),branch_mid %>% 
                         mutate(stock_state="Intermediate"),branch_low %>% 
                         mutate(stock_state="Low"))

Table16_4 <- Table16_4 %>% select(stock_state,scale_N,F,TAC,SSB,SSB_P05,p_below,branch)

Table16_4

stopifnot(Table16_4$branch[Table16_4$stock_state=="High"]=="Fcap")
stopifnot(Table16_4$branch[Table16_4$stock_state=="Intermediate"]=="Reduced F")
stopifnot(Table16_4$branch[Table16_4$stock_state=="Low"]=="Closure")

cat("PASS: Fcap branch identified at high stock state\n")
cat("PASS: reduced-F branch identified at intermediate stock state\n")
cat("PASS: closure branch identified at low stock state\n")

# ============================================================ 
# 15. METHODOLOGICAL TABLES AND FIGURES 
# ============================================================ 
# Builds the principal diagnostics required to document the method 
# and support the subsequent selection of candidate MPs. 
# 
# Outputs describe: 
# - validation against SS3; 
# - effect of uncertainty on required SSB; 
# - stock and catch response to F; 
# - diminishing catch returns at high F; 
# - probabilistic escapement boundaries; 
# - HCR response to forecast uncertainty; 
# - seasonal catch response; 
# - frequency of each decision branch across bootstrap CVs. 
# 
# Branch frequencies in this validation are conditional responses 
# to the 100 forecast-CV values at a fixed stock state. They are not 
# closed-loop probabilities of Fcap, reduced F or closure. 
# ============================================================

Table16_5 <- HCR_S2 %>% count(branch,name="n") %>% mutate(percent=100*n/sum(n))

Table16_5

# ============================================================ 
# 15A. STOCK AND CATCH RESPONSE TO F 
# ============================================================ 
# Shows the biological and fishery response used to characterise 
# the informative Fcap range. SSB decreases with increasing F, 
# whereas annual catch increases with progressively smaller gains. 
# ============================================================

fig16_1a <- ggplot(Fcap_range_diag,aes(x=F,y=SSB)) +
  geom_line(linewidth=0.9) +
  geom_point(size=2) +
  labs(x="Annual fishing mortality",y="Spawning-stock biomass (t)") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank())

fig16_1b <- ggplot(Fcap_range_diag,aes(x=F,y=Catch_total)) +
  geom_line(linewidth=0.9) +
  geom_point(size=2) +
  labs(x="Annual fishing mortality",y="Annual catch (t)") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank())

# ============================================================ 
# 15B. DIMINISHING CATCH RETURNS 
# ============================================================ 
# Quantifies the reduction in incremental catch gains as F 
# increases. This diagnostic helps identify the region where 
# additional increases in F provide progressively less catch. 
# ============================================================

fig16_2 <- ggplot(Fcap_range_diag %>% filter(!is.na(catch_gain_pct)),aes(x=F,y=catch_gain_pct)) +
  geom_line(linewidth=0.9) +
  geom_point(size=2) +
  labs(x="Annual fishing mortality",y="Incremental catch gain (%)") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank())

# ============================================================ 
# 15C. PROBABILISTIC ESCAPEMENT CONSTRAINT 
# ============================================================ 
# Displays the lower 5% forecast SSB quantile across F under S0, 
# S1 and S2-P50. 
# 
# The intersection with Besc represents the maximum F compatible 
# with the 5% escapement-risk criterion under each uncertainty 
# scenario. 
# ============================================================
fig16_3 <- ggplot(risk_response,aes(x=F,y=SSB_P05,linetype=scenario)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept=Besc,linetype="dashed",linewidth=0.7) +
  labs(x="Annual fishing mortality",y="5th percentile of forecast SSB (t)",linetype="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="top",panel.grid.minor=element_blank())

# ============================================================ 
# 15D. HCR RESPONSE TO FORECAST UNCERTAINTY 
# ============================================================ 
# Shows how F advice changes across the empirical bootstrap CV 
# distribution while biological conditions remain fixed. 
# 
# The figure visualises the transition from Fcap to reduced F and, 
# at sufficiently high forecast uncertainty, to closure. 
# ============================================================

fig16_4 <- ggplot(HCR_S2,aes(x=CV,y=F,shape=branch)) +
  geom_point(size=2.3) +
  geom_hline(yintercept=Fcap,linetype="dashed",linewidth=0.6) +
  labs(x="Forecast SSB CV",y="Fishing mortality advice",shape="HCR branch") +
  theme_bw(base_size=12) +
  theme(legend.position="top",panel.grid.minor=element_blank())

# ============================================================ 
# 15E. SEASONAL CATCH RESPONSE 
# ============================================================ 
# Examines the contribution of each quarter to the annual catch 
# response as F increases. 
# 
# This diagnostic helps identify seasonal saturation or depletion 
# effects that are not apparent from annual catch alone. 
# ============================================================

seasonal_catch <- Fcap_range_diag %>% 
                  select(F,Catch_Q1,Catch_Q2,Catch_Q3,Catch_Q4) %>% 
                  pivot_longer(cols=starts_with("Catch_Q"),names_to="season",values_to="Catch")

seasonal_catch <- seasonal_catch %>% mutate(season=recode(season,Catch_Q1="Q1",
                                                          Catch_Q2="Q2",
                                                          Catch_Q3="Q3",
                                                          Catch_Q4="Q4"))

fig16_5 <- ggplot(seasonal_catch,aes(x=F,y=Catch,linetype=season)) +
  geom_line(linewidth=0.9) +
  geom_point(size=1.8) +
  labs(x="Annual fishing mortality",y="Seasonal catch (t)",linetype="Season") +
  theme_bw(base_size=12) +
  theme(legend.position="top",panel.grid.minor=element_blank())

# ============================================================ 
# 16. SAVE VALIDATION OUTPUTS 
# ============================================================ 
# Saves all numerical diagnostics and report-ready figures so that 
# methodological results can be reproduced without rerunning or 
# manually reconstructing intermediate calculations. 
# ============================================================
write.csv(Table16_1,file.path(table_dir,"Table16_1_SS3_probability_validation.csv"),row.names=FALSE)
write.csv(Table16_2,file.path(table_dir,"Table16_2_required_SSB_by_uncertainty.csv"),row.names=FALSE)
write.csv(Table16_3,file.path(table_dir,"Table16_3_HCR_scenario_response.csv"),row.names=FALSE)
write.csv(Table16_4,file.path(table_dir,"Table16_4_HCR_branch_validation.csv"),row.names=FALSE)
write.csv(Table16_5,file.path(table_dir,"Table16_5_HCR_branch_frequency.csv"),row.names=FALSE)
write.csv(SS3_projection_validation,file.path(table_dir,"Table16_SS3_projection_validation.csv"),row.names=FALSE)
write.csv(Table16_Fresponse,file.path(table_dir,"Table16_F_response_by_uncertainty.csv"),row.names=FALSE)
write.csv(Fcap_range_diag,file.path(table_dir,"Table16_Fcap_range_diagnostic.csv"),row.names=FALSE)
write.csv(high_F_diag,file.path(table_dir,"Table16_seasonal_F_diagnostic.csv"),row.names=FALSE)
write.csv(HCR_S2,file.path(table_dir,"Table16_S2_bootstrap_response.csv"),row.names=FALSE)
write.csv(branch_scan,file.path(table_dir,"Table16_stock_state_branch_scan.csv"),row.names=FALSE)

ggsave(file.path(fig_dir,"Fig16_1a_SSB_response_to_F.png"),fig16_1a,width=7,height=5,dpi=300)
ggsave(file.path(fig_dir,"Fig16_1b_catch_response_to_F.png"),fig16_1b,width=7,height=5,dpi=300)
ggsave(file.path(fig_dir,"Fig16_2_diminishing_catch_returns.png"),fig16_2,width=7,height=5,dpi=300)
ggsave(file.path(fig_dir,"Fig16_3_probabilistic_escapement_constraint.png"),fig16_3,width=8,height=6,dpi=300)
ggsave(file.path(fig_dir,"Fig16_4_HCR_response_to_forecast_uncertainty.png"),fig16_4,width=8,height=6,dpi=300)
ggsave(file.path(fig_dir,"Fig16_5_seasonal_catch_response.png"),fig16_5,width=8,height=6,dpi=300)

# ============================================================ 
# 17. FINAL VALIDATION SUMMARY 
# ============================================================ 
# Prints a concise record of the validation checks completed 
# successfully. The script should reach this section only after all 
# numerical consistency tests and stopifnot() checks have passed. 
# ============================================================
cat("\n============================================================\n")
cat("PROBABILISTIC HCR VALIDATION COMPLETED\n")
cat("============================================================\n")
cat("PASS: probability calculation reproduces SS3 benchmark\n")
cat("PASS: one-year projection reproduces SS3 first-year catch benchmark\n")
cat("PASS: P(SSB < Besc) <= 0.05 is equivalent to P05 >= Besc\n")
cat("PASS: required mean SSB increases with forecast uncertainty\n")
cat("PASS: SSB decreases and catch increases over tested F range\n")
cat("PASS: analytical and numerical HCR branch boundaries agree\n")
cat("PASS: CV = 0 reproduces deterministic HCR\n")
cat("PASS: S1 satisfies probabilistic escapement constraint\n")
cat("PASS: S2 evaluated across 100 contemporary bootstrap CVs\n")
cat("PASS: Fcap, reduced-F and closure branches validated\n")
cat("Tables saved to:",table_dir,"\n")
cat("Figures saved to:",fig_dir,"\n")
cat("============================================================\n")

