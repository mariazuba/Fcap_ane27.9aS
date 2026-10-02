# Run from R/RStudio; each scenario runs in its own R process.
# Six runs: REF + five alternatives, all using global trajectory 1.
# MP 10: Fcap=1.50, Besc=6561. Full projection 2025-2054.
# No figures. Stop immediately if any scenario fails.
project_dir <- here::here()
runner <- file.path(project_dir,"scripts","27_RUN_MP_ROBUSTNESS.R")
stopifnot(file.exists(runner))
Rscript_bin <- file.path(R.home("bin"),"Rscript")
log_dir <- file.path(project_dir,"outputs","mse","robustness","validation","test","n1","logs")
dir.create(log_dir,recursive=TRUE,showWarnings=FALSE)
scenario_ids <- c("REF","R-OM","S-C1","S-C2","I0","I1")
old_workdir <- getwd()
setwd(project_dir)

tryCatch({
  for(scenario in scenario_ids) {
    cat("\nTesting",scenario,"- global trajectory 1 - MP 10\n")
    log_file <- file.path(log_dir,paste0(scenario,"_MP10_block01.log"))
    status <- system2(
      Rscript_bin,
      args=c(shQuote(runner),scenario,"10","1","test"),
      env="MSE_N_TEST=1",
      stdout=log_file,
      stderr=log_file
    )
    cat(paste(readLines(log_file,warn=FALSE),collapse="\n"),"\n")
    if(status!=0) stop("Scenario ",scenario," failed. Log: ",log_file)
  }
  cat("\nSIX SINGLE-TRAJECTORY RUNS COMPLETED.\n")
  cat("Execution/structural checks passed; comparative validation remains pending.\n")
}, finally=setwd(old_workdir))
