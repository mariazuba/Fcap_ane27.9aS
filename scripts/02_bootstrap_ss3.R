# ============================================================
# 02_bootstrap_ss3.R
# Bootstrap conditioning of the SS3 operating model
#
# Purpose:
#   Generate SS3 parametric bootstrap datasets, combine them
#   with bootstrap weight-at-age inputs, run the SS3 models,
#   validate convergence, and retain 100 valid historical
#   operating-model conditions for the MSE.
#
# INPUTS:
#   boot/data/ss3_model/
#     - Complete base SS3 model
#     - SS3 executable (ss3_linux)
#
#   data/bootstrap/wtatage_boot/
#     - wtatage_boot_001.ss ... wtatage_boot_152.ss
#
# REQUIRED SCRIPT:
#   scripts/03_read_ss3_to_flbeia.R
#
# OUTPUTS:
#   data/bootstrap/
#     - data_boot_*.ss
#     - bootstrap_runs/boot_*/
#
#   data/Rdata/
#     - stk_boot_100.rds
#     - iteration_map.csv
#     - bootstrap_plot_data.rds
#
# Notes:
#   - SS3 runs are controlled with RUN_SS3 and RUN_bootSS3.
#   - Bootstrap models are checked using maximum gradient
#     and Hessian diagnostics.
#   - 100 valid bootstrap models are retained.
#
# WARNING:
# This script prepares and overwrites SS3 input files even when
# RUN_SS3 and RUN_bootSS3 are FALSE.
#
# It also regenerates consolidated objects and iteration_map.csv.
# Do not run it merely to load existing conditioning objects.
#
# Existing SS3 bootstrap runs must be restored from the backup
# before processing their results.
#
# The current dataset contains 150 bootstrap data files and
# 150 matching weight-at-age files. The starter setting is 152.
# Preserve the existing bootstrap datasets to reproduce the
# conditioning used in the current MSE.
# ============================================================

library(icesTAF)
library(r4ss)
library(here)
library(tidyverse)
library(scales)
library(FLBEIA)

#Paso 1: crear carpeta donde ejecutaremos el bootstrap
run_dir <- here("boot", "data", "ss3_model")
boot_dir <-here("data", "bootstrap")
wt_boot_dir <- here("data","bootstrap","wtatage_boot")
runs_boot_dir <- file.path(boot_dir, "bootstrap_runs")

dir.create(runs_boot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(boot_dir, recursive = TRUE, showWarnings = FALSE)
r4ss::copy_SS_inputs(dir.old = run_dir, dir.new = boot_dir, copy_exe = TRUE, verbose = FALSE)


#Paso 2: Abrir y modificar archivo starter.ss

starter <- SS_readstarter(file.path(boot_dir, "starter.ss"),verbose = FALSE)
starter$N_bootstraps <- 152

SS_writestarter(starter,dir = boot_dir,file = "starter.ss",
                overwrite = TRUE,verbose = TRUE,warn = lifecycle::deprecated())

# Paso 3: GENERATE MODEL RUNS (RUN ONCE) ---------------------------------------------
# These blocks are disabled by default. Set RUN_SS3 <- TRUE to recreate all runs.

RUN_SS3 <- FALSE

if (RUN_SS3) {
  Sys.chmod(file.path(boot_dir, "ss3_linux"), mode = "0755")
  run(dir = boot_dir,exe = "ss3_linux",extras = "",
      skipfinished = TRUE,show_in_console = TRUE,verbose = TRUE)
}

# reviso que ss3 haya generado los 100 data_boot
# list.files(boot_dir)

# Creo una cartepa que guardarás las corridas de cada boot dentro de la carpeta bootstrap

# creo una lista de los archivos boot 
boot_files <- list.files(boot_dir,pattern = "^data_boot_[0-9]{3}\\.ss$",full.names = TRUE)
# creo el directorio para cada bootstrap
for(i in seq_along(boot_files)) {
  
  boot_name <- sprintf("boot_%03d",i)
  
  # directorio del bootstrap i
  boot_i_dir <- file.path(runs_boot_dir,boot_name)
  # Copiar inputs del modelo SS3 original
  r4ss::copy_SS_inputs(dir.old = run_dir,dir.new = boot_i_dir,copy_exe = TRUE,verbose = FALSE)
  # Reemplazar data.SS por el bootstrap correspondiente
  file.copy(from = boot_files[i],to = file.path(boot_i_dir, "data.SS"),overwrite = TRUE)
 
  # Bootstrap weight-at-age
  wt_file <- file.path(wt_boot_dir,paste0("wtatage_",boot_name,".ss"))
  file.copy(from = wt_file,to = file.path(boot_i_dir,"wtatage.ss"),overwrite = TRUE)
  # Modificar starter
  starter <- r4ss::SS_readstarter(file.path(boot_i_dir, "starter.ss"),verbose = FALSE)
  starter$datfile <- "data.SS"
  starter$N_bootstraps <- 0
  
  r4ss::SS_writestarter(starter, dir = boot_i_dir,overwrite = TRUE)
  
  cat("Prepared bootstrap", boot_name, "\n")
}

# correr  primer boot 
# Sys.chmod(file.path(runs_boot_dir, "boot_001", "ss3_linux"),mode = "0755")
#  r4ss::run(dir = file.path(runs_boot_dir, "boot_001"),exe = "ss3_linux",extras = "",
#            skipfinished = TRUE,show_in_console = TRUE,verbose = TRUE)

#'*============================================================================*
#'*============================================================================*
# Ejecuto todos los bootstrap
RUN_bootSS3 <- FALSE


if (RUN_bootSS3) {
  
  boot_models <- list.dirs(runs_boot_dir,recursive = FALSE,full.names = TRUE)
  # ordenar boot_001, boot_002, ..., boot_100
  boot_models <- boot_models[order(basename(boot_models))]
  
  for (i in seq_along(boot_models)) {
    model <- boot_models[i]
    cat(
      "\n----------------------------------\n",
      "Running ", basename(model),
      " (", i, "/", length(boot_models), ")\n",
      "----------------------------------\n",
      sep = "")
    
    Sys.chmod(file.path(model, "ss3_linux"),mode = "0755")
    r4ss::run(dir = model,exe = "ss3_linux",extras = "",
              skipfinished = TRUE,show_in_console = FALSE,verbose = TRUE)
  }
}
boot_models <- list.dirs(runs_boot_dir,recursive = FALSE,full.names = TRUE)
#################################################################################
# Ahora chekeo los 100 bootstrap 
# Leo los reportes, chekeo el gradiente maximo 
run_check <- purrr::map_dfr(
  boot_models,
  function(model) {
    console_file <- file.path(model,"console.output.txt")
    txt <- if (file.exists(console_file)) {readLines(console_file, warn = FALSE)
    } else {character()}
    tibble(
      bootstrap = basename(model),
      console_exists =file.exists(console_file),
      run_completed =any(grepl("Run has completed", txt, ignore.case = TRUE)),
      report_exists =file.exists(file.path(model, "Report.sso")),
      covar_exists =file.exists(file.path(model, "covar.sso")),
      cor_exists = file.exists(file.path(model, "ss3.cor")))
  }
)

run_check %>% summarise(
    n = dplyr::n(),
    completed = sum(run_completed),
    report = sum(report_exists),
    covar = sum(covar_exists),
    cor = sum(cor_exists))

run_check %>% filter(!run_completed |!report_exists) %>% print(n = Inf)
boot_check <- vector("list", length(boot_models))

for(i in seq_along(boot_models)) {
  model <- boot_models[i]
  cat("Checking",basename(model),i,"/",length(boot_models),"\n")
  out <- tryCatch(r4ss::SS_output(dir = model,covar = FALSE,verbose = FALSE),error = function(e) NULL)
  if(is.null(out)) {
    boot_check[[i]] <- tibble(bootstrap = basename(model),read_ok = FALSE,max_gradient = NA_real_)
  } else {
    boot_check[[i]] <- tibble(bootstrap = basename(model),read_ok = TRUE,max_gradient = out$maximum_gradient_component)
  }
  rm(out)
  # Más frecuente para cuidar memoria
  if(i %% 5 == 0) {gc()}
}
#'*============================================================================*
#'*============================================================================*
boot_check <- bind_rows(boot_check)
summary(boot_check$max_gradient)

boot_check %>%filter(is.na(max_gradient) |max_gradient > 0.001) %>%
  arrange(desc(max_gradient)) %>%print(n = Inf)

# Ahora checkeo la hessiana

hessian_check <- purrr::map_dfr(boot_models,
  function(model) {
    console_file <- file.path(model,"console.output.txt")
    txt <- if (file.exists(console_file)) {
      readLines(console_file, warn = FALSE)} else {character()}
    
    tibble(bootstrap = basename(model),
           hessian_started = any(grepl("Inverting Hessian",txt,ignore.case = TRUE)),
           se_calculated = any(grepl("Starting standard error calculations", txt, ignore.case = TRUE)),
           run_completed = any(grepl("Run has completed", txt, ignore.case = TRUE)),
           covar_exists = file.exists(file.path(model, "covar.sso")),
           cor_exists = file.exists(file.path(model, "ss3.cor")))
  })

hessian_check %>% summarise(n = dplyr::n(),
                            hessian_started = sum(hessian_started),
                            se_calculated = sum(se_calculated),
                            run_completed = sum(run_completed),
                            covar_exists = sum(covar_exists),
                            cor_exists = sum(cor_exists))

boot_diagnostics <- boot_check %>%left_join(hessian_check,by = "bootstrap")
boot_diagnostics <- boot_diagnostics %>%
                    mutate(gradient_ok = max_gradient <= 0.001,
                           hessian_ok = hessian_started & se_calculated & covar_exists,
                           model_ok = read_ok & gradient_ok & hessian_ok)

table(boot_diagnostics$model_ok)

boot_failed <- boot_diagnostics %>%
               filter(!model_ok) %>%
               select(bootstrap,max_gradient,gradient_ok,hessian_started,
                      se_calculated, covar_exists,cor_exists,run_completed,model_ok) %>%
               arrange(desc(max_gradient))

set.seed(123)

valid_boot <- boot_diagnostics %>% filter(model_ok) %>%slice_sample(n = 100) %>% select(bootstrap)
valid_boot <- valid_boot %>%arrange(bootstrap)

boot_files_100 <- boot_files[str_extract(basename(boot_files),"boot_[0-9]{3}") %in% valid_boot$bootstrap]
#ordenar por ID original
boot_files_100 <- boot_files_100[order(str_extract(basename(boot_files_100),"boot_[0-9]{3}"))]

length(boot_files_100)

# # Primer bootstrap válido seleccionado
# test_boot <- boot_models_100[1]
# 
# test_boot
# basename(test_boot)
# 
# #llamamos la función para convertir ss3 a flstock
# test_obj <- read_ss3_to_flbeia(wd = test_boot,output_name = paste0("ane9aS_",basename(test_boot)),rds_dir = rds_boot_dir)
boot_models_100 <- boot_models[basename(boot_models) %in% valid_boot$bootstrap]
rds_boot_dir <- file.path(getwd(),"data","Rdata")
dir.create(rds_boot_dir,recursive = TRUE,showWarnings = FALSE)

source(file.path(getwd(),"scripts","03_read_ss3_to_flbeia.R"))

stk_boot <- vector("list",length(boot_models_100))
names(stk_boot) <- basename(boot_models_100)

for(i in seq_along(boot_models_100)) {
  model_dir <- boot_models_100[i]
  boot_id <- basename(model_dir)
  cat("\nConverting ",boot_id," (",i,"/",length(boot_models_100),")\n",sep="")
  obj <- read_ss3_to_flbeia(wd=model_dir,output_name=paste0("ane9aS_",boot_id),rds_dir=rds_boot_dir)
  stk_boot[[boot_id]] <- obj$stk
  rm(obj)
  if(i %% 5 == 0) gc()
}

saveRDS(stk_boot,file.path(rds_boot_dir,"stk_boot_100.rds"))

iteration_map <- tibble(iter=seq_along(stk_boot),bootstrap=names(stk_boot))
write.csv(iteration_map,file.path(rds_boot_dir,"iteration_map.csv"),row.names=FALSE)

# ============================================================
# DATA FOR BOOTSTRAP PLOTS - 100 VALIDATED BOOTSTRAPS
# ============================================================

stopifnot(length(boot_files_100)==100)
stopifnot(length(boot_models_100)==100)
stopifnot(setequal(str_extract(basename(boot_files_100),"boot_[0-9]{3}"),basename(boot_models_100)))
# Leer output del modelo base
base_out <- r4ss::SS_output(dir=boot_dir,covar=FALSE,verbose=FALSE)

base_indices <- base_out$cpue %>% filter(Fleet %in% 5:8) %>% transmute(index_name=Fleet_name,year=Yr,obs=Obs,exp=Exp)
base_catch <- base_out$catch %>% filter(Fleet %in% 1:4,Yr>=1989) %>% transmute(seas=Seas,fleet=Fleet,year=Yr,obs=Obs,exp=Exp)
base_agecomp <- base_out$agedbase %>% filter(Fleet %in% 1:8,Used=="yes") %>% mutate(fleet_name=recode(as.character(Fleet),`1`="SEINE_Q1",`2`="SEINE_Q2",`3`="SEINE_Q3",`4`="SEINE_Q4",`5`="PELAGO",`6`="ECOCADIZ",`7`="BOCADEVA",`8`="ECORECLUTAS"),age=factor(Bin,levels=0:3,labels=c("0","1","2","3+"))) %>% select(year=Yr,fleet_name,age,observed=Obs,expected=Exp)
base_natage <- base_out$natage %>% filter(Yr>=1989,`Beg/Mid`=="B") %>% select(Yr,Seas,`0`:`3`) %>% pivot_longer(`0`:`3`,names_to="age",values_to="N")
base_ssb <- base_out$timeseries %>% filter(Yr>=1989,Seas==2,is.finite(SpawnBio)) %>% select(Yr,Seas,SSB=SpawnBio)
base_rec <- base_out$timeseries %>% filter(Yr>=1989,Seas==3,is.finite(Recruit_0)) %>% select(Yr,Seas,R=Recruit_0)
base_F <- base_out$timeseries %>% filter(Yr>=1989) %>% select(Yr,Seas,`F:_1`:`F:_4`) %>% pivot_longer(`F:_1`:`F:_4`,names_to="fleet",values_to="F") %>% mutate(fleet=recode(fleet,`F:_1`="SEINE_Q1",`F:_2`="SEINE_Q2",`F:_3`="SEINE_Q3",`F:_4`="SEINE_Q4")) %>% filter((fleet=="SEINE_Q1" & Seas==1) | (fleet=="SEINE_Q2" & Seas==2) | (fleet=="SEINE_Q3" & Seas==3) | (fleet=="SEINE_Q4" & Seas==4))
base_sel <- base_out$ageselex %>% filter(Factor=="Asel",Fleet %in% 1:8) %>% pivot_longer(`0`:`3`,names_to="age",values_to="Sel") %>% mutate(fleet_name=recode(as.character(Fleet),`1`="SEINE_Q1",`2`="SEINE_Q2",`3`="SEINE_Q3",`4`="SEINE_Q4",`5`="PELAGO",`6`="ECOCADIZ",`7`="BOCADEVA",`8`="ECORECLUTAS"),age=factor(age,levels=c("0","1","2","3")))
base_q <- base_out$estimated_non_dev_parameters %>% filter(grepl("^LnQ_base_",rownames(.))) %>% mutate(parameter=rownames(.),index_name=case_when(grepl("PELAGO",parameter)~"PELAGO",grepl("ECOCADIZ",parameter)~"ECOCADIZ",grepl("BOCADEVA",parameter)~"BOCADEVA",grepl("ECORECLUTAS",parameter)~"ECORECLUTAS"),q=exp(Value)) %>% select(index_name,lnq=Value,q)

# Datos bootstrap de entrda
boot_indices <- vector("list",length(boot_files_100))
boot_catch <- vector("list",length(boot_files_100))
boot_agecomp_fleet <- vector("list",length(boot_files_100))

for(i in seq_along(boot_files_100)) {
  dat <- r4ss::SS_readdat(file=boot_files_100[i],verbose=FALSE)
  boot_id <- str_extract(basename(boot_files_100[i]),"boot_[0-9]{3}")
  boot_indices[[i]] <- dat$CPUE %>% mutate(bootstrap=boot_id)
  boot_catch[[i]] <- dat$catch %>% filter(year>=1989,fleet==seas) %>% mutate(bootstrap=boot_id)
  boot_agecomp_fleet[[i]] <- dat$agecomp %>% filter(fleet %in% 1:8) %>% mutate(bootstrap=boot_id,fleet_name=recode(as.character(fleet),`1`="SEINE_Q1",`2`="SEINE_Q2",`3`="SEINE_Q3",`4`="SEINE_Q4",`5`="PELAGO",`6`="ECOCADIZ",`7`="BOCADEVA",`8`="ECORECLUTAS"))
  rm(dat)
  if(i %% 20==0) gc()
}

boot_indices <- bind_rows(boot_indices) %>% mutate(index_name=recode(as.character(index),`5`="PELAGO",`6`="ECOCADIZ",`7`="BOCADEVA",`8`="ECORECLUTAS"))
boot_catch <- bind_rows(boot_catch)
boot_agecomp_fleet <- bind_rows(boot_agecomp_fleet)
boot_agecomp_long <- boot_agecomp_fleet %>% pivot_longer(cols=a0:a3,names_to="age",values_to="N") %>% mutate(age=factor(age,levels=c("a0","a1","a2","a3"),labels=c("0","1","2","3+"))) %>% group_by(bootstrap,year,fleet,fleet_name) %>% mutate(prop=N/sum(N)) %>% ungroup()

# Outputs de los 100 SS3 validados
boot_natage <- vector("list",length(boot_models_100))
boot_ssb <- vector("list",length(boot_models_100))
boot_rec <- vector("list",length(boot_models_100))
boot_F <- vector("list",length(boot_models_100))
boot_sel <- vector("list",length(boot_models_100))
boot_q <- vector("list",length(boot_models_100))

#'*============================================================================*
#'*============================================================================*
for(i in seq_along(boot_models_100)) {

  out <- r4ss::SS_output(boot_models_100[i], covar = FALSE, verbose = FALSE)

  # N-at-age
  boot_natage[[i]] <- out$natage %>% filter(Yr >= 1989, `Beg/Mid` == "B") %>% select(Yr, Seas, `0`:`3`) %>%
                      pivot_longer(`0`:`3`, names_to = "age", values_to = "N") %>%
                       mutate(bootstrap=basename(boot_models_100[i]))
  # SSB
  boot_ssb[[i]] <- out$timeseries %>% filter(Yr>=1989,Seas==2,is.finite(SpawnBio)) %>% 
                    select(Yr,Seas,SSB=SpawnBio) %>% mutate(bootstrap=basename(boot_models_100[i]))
  
  # Recruitment
  boot_rec[[i]] <- out$timeseries %>% filter(Yr>=1989,Seas==3,is.finite(Recruit_0)) %>% 
                    select(Yr,Seas,R=Recruit_0) %>% mutate(bootstrap=basename(boot_models_100[i]))
  
  # Mortalidad por pesca
  boot_F[[i]] <- out$timeseries %>%filter(Yr >= 1989) %>%select(Yr, Seas, `F:_1`:`F:_4`) %>%
                 pivot_longer(`F:_1`:`F:_4`,names_to = "fleet",values_to = "F") %>%
                 mutate(bootstrap = basename(boot_models_100[i]),
                        fleet = recode(fleet,
                                       `F:_1` = "SEINE_Q1",
                                        `F:_2` = "SEINE_Q2",
                                         `F:_3` = "SEINE_Q3",
                                          `F:_4` = "SEINE_Q4")) %>%
                 filter((fleet == "SEINE_Q1" & Seas == 1) |(fleet == "SEINE_Q2" & Seas == 2) |
                        (fleet == "SEINE_Q3" & Seas == 3) |(fleet == "SEINE_Q4" & Seas == 4))

  # Selectividad
  boot_sel[[i]] <- out$ageselex %>%
                   filter(Factor == "Asel", Fleet %in% 1:8) %>%
                   pivot_longer(`0`:`3`, names_to = "age", values_to = "Sel") %>%
                   mutate(bootstrap = basename(boot_models_100[i]),
                   fleet_name = recode(as.character(Fleet),
                               `1`="SEINE_Q1", `2`="SEINE_Q2", `3`="SEINE_Q3", `4`="SEINE_Q4",
                               `5`="PELAGO", `6`="ECOCADIZ", `7`="BOCADEVA", `8`="ECORECLUTAS"),
                   age = factor(age, levels = c("0","1","2","3")))
  # Catchability
  pars        <- out$estimated_non_dev_parameters
  boot_q[[i]] <- pars %>%filter(grepl("^LnQ_base_", rownames(.))) %>%
                 mutate(parameter = rownames(.),
                        bootstrap = basename(boot_models_100[i]),
                        index_name = case_when(grepl("PELAGO", parameter) ~ "PELAGO",
                                  grepl("ECOCADIZ", parameter) ~ "ECOCADIZ",
                                  grepl("BOCADEVA", parameter) ~ "BOCADEVA",
                                  grepl("ECORECLUTAS", parameter) ~ "ECORECLUTAS"),
                        q = exp(Value)) %>%
                select(bootstrap, index_name, lnq = Value, q)

  rm(out, pars)
  if(i %% 10 == 0) gc()
}


boot_natage <- bind_rows(boot_natage)
boot_ssb <- bind_rows(boot_ssb)
boot_rec <- bind_rows(boot_rec)
boot_F <- bind_rows(boot_F)
boot_sel <- bind_rows(boot_sel)
boot_q <- bind_rows(boot_q)

# ============================================================
# SAVE DATA FOR BOOTSTRAP PLOTS
# ============================================================
boot_plot_data <- list(
  base_indices=base_indices,
  boot_indices=boot_indices,
  base_catch=base_catch,
  boot_catch=boot_catch,
  base_agecomp=base_agecomp,
  boot_agecomp=boot_agecomp_long,
  base_natage=base_natage,
  boot_natage=boot_natage,
  base_ssb=base_ssb,
  boot_ssb=boot_ssb,
  base_rec=base_rec,
  boot_rec=boot_rec,
  base_F=base_F,
  boot_F=boot_F,
  base_sel=base_sel,
  boot_sel=boot_sel,
  base_q=base_q,
  boot_q=boot_q
)

saveRDS(boot_plot_data,file.path(rds_boot_dir,"bootstrap_plot_data.rds"))
