# ============================================================
# 17_TEST_PROBABILISTIC_HCR_FUNCTION.R
#
# Unit test of the deterministic/probabilistic Fcap-Besc HCR
# outside FLBEIA.
#
# Objectives:
#   1. Verify deterministic mode.
#   2. Verify probabilistic mode with CV = 0.
#   3. Verify S1 current-assessment uncertainty.
#   4. Verify S2 median contemporary uncertainty.
#   5. Verify Fcap / Reduced F / Closure branches.
#
# This script DOES NOT run FLBEIA.
# ============================================================
rm(list=ls())

library(FLCore)
library(FLBEIA)
library(dplyr)
library(here)

# ============================================================
# DIRECTORIES
# ============================================================

rds_dir <- here("data","Rdata")

# ============================================================
# LOAD REFERENCE OBJECTS
# ============================================================

fleets_ref <- readRDS(file.path(rds_dir,"fleets_conditioned_1000iter.rds"))
SRs_ref <- readRDS(file.path(rds_dir,"SRs_reference_conditioned_1000iter.rds"))

# ============================================================
# LOAD REFERENCE OM
# ============================================================

OM_test <- readRDS(here("data","mse","OM_validation_C0","OM_validation_C0_block_01.rds"))
# ============================================================
# SETTINGS
# ============================================================

Fcap <- 1.5
Besc <- 6561
alpha <- 0.05


forecast_uncertainty <- readRDS(file.path(rds_dir,"forecast_uncertainty_characterisation.rds"))

CV_S0 <- 0
CV_S1 <- forecast_uncertainty$base$SSB_CV
CV_S2 <- median(forecast_uncertainty$contemporary$SSB_2025_CV)

stopifnot(length(CV_S1)==1)
stopifnot(is.finite(CV_S1))
stopifnot(CV_S1>0)

stopifnot(length(CV_S2)==1)
stopifnot(is.finite(CV_S2))
stopifnot(CV_S2>0)

test_iter <- 1
ass.yr <- 2024
tac.yr <- 2025

# ============================================================
# BUILD REFERENCE INPUTS
# ============================================================

q_test <- OM_test$fleets$SEINE@metiers$ALL@catches$ANE@catch.q
q_test <- q_test[,as.character(tac.yr),,,,test_iter]

Frel <- sapply(1:4,function(s) as.numeric(q_test[,,,as.character(s),,])/max(as.numeric(q_test[,,,as.character(s),,])))
rownames(Frel) <- 0:3
colnames(Frel) <- paste0("Q",1:4)

propf <- as.numeric(OM_test$fleets.ctrl$seasonal.share$ANE[,as.character(tac.yr),,,,test_iter])

sr_params <- OM_test$SRs$ANE@params
a <- as.numeric(sr_params["a",as.character(tac.yr),"3",test_iter])
b <- as.numeric(sr_params["b",as.character(tac.yr),"3",test_iter])

M <- sapply(1:4,function(s) as.numeric(OM_test$biols$ANE@m[,as.character(tac.yr),,as.character(s),,test_iter]))
rownames(M) <- 0:3
colnames(M) <- paste0("Q",1:4)

mat_pred <- predict(OM_test$biols$ANE@mat)
wt_Q2 <- as.numeric(OM_test$biols$ANE@wt[,as.character(tac.yr),,"2",,test_iter])
mat_Q2 <- as.numeric(mat_pred[,as.character(tac.yr),,"2",,test_iter])
W_spawn <- wt_Q2*mat_Q2
names(W_spawn) <- 0:3

catch_test <- OM_test$fleets$SEINE@metiers$ALL@catches$ANE
W_catch <- sapply(1:4,function(s) as.numeric(catch_test@landings.wt[,as.character(tac.yr),,as.character(s),,test_iter]))
rownames(W_catch) <- 0:3
colnames(W_catch) <- paste0("Q",1:4)

N_Q1 <- as.numeric(OM_test$biols$ANE@n[,as.character(tac.yr),,"1",,test_iter])
names(N_Q1) <- 0:3

# ============================================================
# SS3-LIKE ONE-YEAR PROJECTION
# ============================================================

project_ss3_one_year <- function(
    N_Q1,
    M,
    Frel,
    W_spawn,
    W_catch,
    Fbase,
    propf,
    a,
    b,
    dt=0.25,
    M_is_seasonal=FALSE
) {
  
  ages <- names(N_Q1)
  seasons <- seq_len(ncol(Frel))
  Fseason <- Fbase*propf
  
  N <- matrix(0,nrow=length(ages),ncol=length(seasons),dimnames=list(age=ages,season=seasons))
  Fage <- matrix(0,nrow=length(ages),ncol=length(seasons),dimnames=list(age=ages,season=seasons))
  Z <- matrix(0,nrow=length(ages),ncol=length(seasons),dimnames=list(age=ages,season=seasons))
  catch_num <- matrix(0,nrow=length(ages),ncol=length(seasons),dimnames=list(age=ages,season=seasons))
  catch_bio <- numeric(length(seasons))
  
  get_Ms <- function(s) if (M_is_seasonal) M[,s] else M
  get_dt <- function() if (M_is_seasonal) 1 else dt
  
  N[,1] <- N_Q1
  dts <- get_dt()
  
  Ms <- get_Ms(1)
  Fage[,1] <- Fseason[1]*Frel[,1]
  Z[,1] <- Ms+Fage[,1]
  catch_num[,1] <- N[,1]*(Fage[,1]/Z[,1])*(1-exp(-Z[,1]*dts))
  catch_bio[1] <- sum(catch_num[,1]*W_catch[,1])
  N[,2] <- N[,1]*exp(-Z[,1]*dts)
  
  SSB_Q2 <- sum(N[,2]*W_spawn)
  R_Q3 <- a*SSB_Q2/(b+SSB_Q2)
  
  Ms <- get_Ms(2)
  Fage[,2] <- Fseason[2]*Frel[,2]
  Z[,2] <- Ms+Fage[,2]
  catch_num[,2] <- N[,2]*(Fage[,2]/Z[,2])*(1-exp(-Z[,2]*dts))
  catch_bio[2] <- sum(catch_num[,2]*W_catch[,2])
  N[,3] <- N[,2]*exp(-Z[,2]*dts)
  N["0",3] <- R_Q3
  
  Ms <- get_Ms(3)
  Fage[,3] <- Fseason[3]*Frel[,3]
  Z[,3] <- Ms+Fage[,3]
  catch_num[,3] <- N[,3]*(Fage[,3]/Z[,3])*(1-exp(-Z[,3]*dts))
  catch_bio[3] <- sum(catch_num[,3]*W_catch[,3])
  N[,4] <- N[,3]*exp(-Z[,3]*dts)
  
  Ms <- get_Ms(4)
  Fage[,4] <- Fseason[4]*Frel[,4]
  Z[,4] <- Ms+Fage[,4]
  catch_num[,4] <- N[,4]*(Fage[,4]/Z[,4])*(1-exp(-Z[,4]*dts))
  catch_bio[4] <- sum(catch_num[,4]*W_catch[,4])
  
  N_endQ4 <- N[,4]*exp(-Z[,4]*dts)
  
  N_next_Q1 <- c(
    `0`=0,
    `1`=unname(N_endQ4["0"]),
    `2`=unname(N_endQ4["1"]),
    `3`=unname(N_endQ4["2"]+N_endQ4["3"])
  )
  
  list(
    Fbase=Fbase,
    Fseason=Fseason,
    N=N,
    Fage=Fage,
    Z=Z,
    SSB_Q2=SSB_Q2,
    R_Q3=R_Q3,
    catch_num=catch_num,
    catch_bio=catch_bio,
    TAC=sum(catch_bio),
    N_endQ4=N_endQ4,
    N_next_Q1=N_next_Q1
  )
}

# ============================================================
# DETERMINISTIC / PROBABILISTIC HCR
# ============================================================

get_F_HCR <- function(
    N_Q1,
    M,
    Frel,
    W_spawn,
    W_catch,
    Fcap,
    Besc,
    propf,
    a,
    b,
    dt=0.25,
    M_is_seasonal=TRUE,
    use.probability=FALSE,
    CV=0,
    alpha=0.05
) {
  
  run_F <- function(F) {
    project_ss3_one_year(
      N_Q1=N_Q1,
      M=M,
      Frel=Frel,
      W_spawn=W_spawn,
      W_catch=W_catch,
      Fbase=F,
      propf=propf,
      a=a,
      b=b,
      dt=dt,
      M_is_seasonal=M_is_seasonal
    )
  }
  
  evaluate_F <- function(F) {
    
    proj <- run_F(F)
    SSB_mean <- proj$SSB_Q2
    
    if (!use.probability) {
      
      SSB_SD <- 0
      SSB_Qalpha <- SSB_mean
      p_below <- as.numeric(SSB_mean < Besc)
      acceptable <- SSB_mean >= Besc
      
    } else {
      
      SSB_SD <- CV*SSB_mean
      
      if (CV == 0) {
        
        SSB_Qalpha <- SSB_mean
        p_below <- as.numeric(SSB_mean < Besc)
        
      } else {
        
        SSB_Qalpha <- qnorm(alpha,mean=SSB_mean,sd=SSB_SD)
        p_below <- pnorm(Besc,mean=SSB_mean,sd=SSB_SD)
      }
      
      acceptable <- SSB_Qalpha >= Besc
    }
    
    list(
      projection=proj,
      SSB_mean=SSB_mean,
      SSB_SD=SSB_SD,
      SSB_Qalpha=SSB_Qalpha,
      p_below=p_below,
      acceptable=acceptable
    )
  }
  
  e0 <- evaluate_F(0)
  ecap <- evaluate_F(Fcap)
  
  if (ecap$acceptable) {
    
    Fstar <- Fcap
    eval_star <- ecap
    branch <- "Fcap"
    
  } else if (!e0$acceptable) {
    
    Fstar <- 0
    eval_star <- e0
    branch <- "closure"
    
  } else {
    
    obj <- function(F) evaluate_F(F)$SSB_Qalpha-Besc
    
    Fstar <- uniroot(obj,c(0,Fcap))$root
    eval_star <- evaluate_F(Fstar)
    branch <- "reduced"
  }
  
  proj <- eval_star$projection
  
  list(
    F=Fstar,
    TAC=proj$TAC,
    SSB_Q2=proj$SSB_Q2,
    R_Q3=proj$R_Q3,
    SSB_SD=eval_star$SSB_SD,
    SSB_Qalpha=eval_star$SSB_Qalpha,
    p_below=eval_star$p_below,
    acceptable=eval_star$acceptable,
    branch=branch,
    use_probability=use.probability,
    CV=CV,
    alpha=alpha,
    SSB_F0=e0$SSB_mean,
    SSB_Fcap=ecap$SSB_mean,
    Qalpha_F0=e0$SSB_Qalpha,
    Qalpha_Fcap=ecap$SSB_Qalpha,
    pbelow_F0=e0$p_below,
    pbelow_Fcap=ecap$p_below,
    projection=proj
  )
}


# ============================================================
# CHECK REFERENCE INPUTS
# ============================================================

stopifnot(length(N_Q1)==4)
stopifnot(all(dim(M)==c(4,4)))
stopifnot(all(dim(Frel)==c(4,4)))
stopifnot(length(W_spawn)==4)
stopifnot(all(dim(W_catch)==c(4,4)))
stopifnot(length(propf)==4)
stopifnot(length(a)==1)
stopifnot(length(b)==1)
stopifnot(abs(sum(propf)-1)<1e-6)

cat("\nREFERENCE INPUTS READY\n")
cat("Iteration:",test_iter,"\n")
cat("Forecast year:",tac.yr,"\n")
cat("a:",a,"\n")
cat("b:",b,"\n")
cat("Seasonal F proportions:",paste(round(propf,4),collapse=", "),"\n")

# ============================================================
# TEST 1 — ORIGINAL DETERMINISTIC HCR
# ============================================================

res_det <- get_F_HCR(
  N_Q1=N_Q1,
  M=M,
  Frel=Frel,
  W_spawn=W_spawn,
  W_catch=W_catch,
  Fcap=Fcap,
  Besc=Besc,
  propf=propf,
  a=a,
  b=b,
  M_is_seasonal=TRUE,
  use.probability=FALSE
)


# ============================================================
# TEST 2 — PROBABILISTIC SWITCH ON, CV = 0
# ============================================================

res_prob0 <- get_F_HCR(
  N_Q1=N_Q1,
  M=M,
  Frel=Frel,
  W_spawn=W_spawn,
  W_catch=W_catch,
  Fcap=Fcap,
  Besc=Besc,
  propf=propf,
  a=a,
  b=b,
  M_is_seasonal=TRUE,
  use.probability=TRUE,
  CV=CV_S0,
  alpha=alpha
)


# ============================================================
# TEST 3 — S1
# ============================================================

res_S1 <- get_F_HCR(
  N_Q1=N_Q1,
  M=M,
  Frel=Frel,
  W_spawn=W_spawn,
  W_catch=W_catch,
  Fcap=Fcap,
  Besc=Besc,
  propf=propf,
  a=a,
  b=b,
  M_is_seasonal=TRUE,
  use.probability=TRUE,
  CV=CV_S1,
  alpha=alpha
)


# ============================================================
# TEST 4 — S2 MEDIAN CV
# ============================================================

res_S2 <- get_F_HCR(
  N_Q1=N_Q1,
  M=M,
  Frel=Frel,
  W_spawn=W_spawn,
  W_catch=W_catch,
  Fcap=Fcap,
  Besc=Besc,
  propf=propf,
  a=a,
  b=b,
  M_is_seasonal=TRUE,
  use.probability=TRUE,
  CV=CV_S2,
  alpha=alpha
)

# ============================================================
# CHECK TESTS 1-4
# ============================================================

stopifnot(isTRUE(all.equal(res_det$F,res_prob0$F,tolerance=1e-10)))
stopifnot(isTRUE(all.equal(res_det$TAC,res_prob0$TAC,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(res_det$SSB_Q2,res_prob0$SSB_Q2,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(res_det$SSB_Qalpha,res_prob0$SSB_Qalpha,tolerance=1e-8)))
stopifnot(res_det$branch==res_prob0$branch)

stopifnot(res_S1$branch=="Fcap")
stopifnot(res_S1$F==Fcap)
stopifnot(res_S1$SSB_Qalpha>=Besc)
stopifnot(res_S1$p_below<=alpha)

stopifnot(res_S2$branch=="reduced")
stopifnot(res_S2$F>0)
stopifnot(res_S2$F<Fcap)
stopifnot(abs(res_S2$SSB_Qalpha-Besc)<1)
stopifnot(abs(res_S2$p_below-alpha)<1e-4)

cat("\nPASS: deterministic and probabilistic HCR tests completed\n")

# ============================================================
# TEST SUMMARY
# ============================================================

test_summary <- data.frame(
  Scenario=c("Deterministic","Probability CV=0","S1","S2 median"),
  CV=c(0,CV_S0,CV_S1,CV_S2),
  F=c(res_det$F,res_prob0$F,res_S1$F,res_S2$F),
  TAC=c(res_det$TAC,res_prob0$TAC,res_S1$TAC,res_S2$TAC),
  SSB=c(res_det$SSB_Q2,res_prob0$SSB_Q2,res_S1$SSB_Q2,res_S2$SSB_Q2),
  Q05=c(res_det$SSB_Qalpha,res_prob0$SSB_Qalpha,res_S1$SSB_Qalpha,res_S2$SSB_Qalpha),
  P_below=c(res_det$p_below,res_prob0$p_below,res_S1$p_below,res_S2$p_below),
  Branch=c(res_det$branch,res_prob0$branch,res_S1$branch,res_S2$branch)
)

print(test_summary)

# ============================================================
# TEST 5 — FCAP / REDUCED F / CLOSURE BRANCHES
# ============================================================

CV_closure_threshold <- (Besc/res_S2$SSB_F0-1)/qnorm(alpha)

CV_closure_test <- CV_closure_threshold+0.01

res_closure <- get_F_HCR(
  N_Q1=N_Q1,
  M=M,
  Frel=Frel,
  W_spawn=W_spawn,
  W_catch=W_catch,
  Fcap=Fcap,
  Besc=Besc,
  propf=propf,
  a=a,
  b=b,
  M_is_seasonal=TRUE,
  use.probability=TRUE,
  CV=CV_closure_test,
  alpha=alpha
)

cat("\nClosure CV threshold:",round(CV_closure_threshold,4),"\n")
cat("Closure test CV:",round(CV_closure_test,4),"\n")
cat("Q05 at F=0:",round(res_closure$Qalpha_F0,2),"\n")
cat("P below Besc at F=0:",round(res_closure$pbelow_F0,4),"\n")
cat("F:",round(res_closure$F,4),"\n")
cat("Branch:",res_closure$branch,"\n")

# ============================================================
# CHECK ALL THREE HCR BRANCHES
# ============================================================

stopifnot(res_S1$branch=="Fcap")
stopifnot(res_S2$branch=="reduced")
stopifnot(res_closure$branch=="closure")

stopifnot(res_S1$F==Fcap)
stopifnot(res_S2$F>0 & res_S2$F<Fcap)
stopifnot(res_closure$F==0)

stopifnot(res_S1$Qalpha_Fcap>=Besc)
stopifnot(abs(res_S2$SSB_Qalpha-Besc)<1)
stopifnot(res_closure$Qalpha_F0<Besc)

cat("\nPASS: Fcap / Reduced F / Closure branches correctly identified\n")

# ============================================================
# FINAL VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("PROBABILISTIC HCR FUNCTION VALIDATION COMPLETED\n")
cat("============================================================\n")