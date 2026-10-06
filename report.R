# ============================================================
# report.R
# Generate figures and render Working Documents.
# Uncomment the required figures or documents.
# Required data, tables and figures must already be available.
# ============================================================

# Figures
# source("scripts/02b_bootstrap_ss3_plots.R")
# source("scripts/23_PLOT_CANDIDATE_MPS_FINAL.R")

# Working Document paths
wd_dir <- here::here("report", "WD")
pdf_dir <- file.path(wd_dir, "pdf")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# Working Documents
 rmarkdown::render(file.path(wd_dir, "00WD_Introduction_MSE_ane27_9aS_ES.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "01WD_OM_conditioning_ES.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "02WD_recruitment_SR_ES.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "03WD_fishery_conditioning_seasonality_ES.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "04WD_reference_OM_open_loop_validation_ES.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "05WD_assessment_error_ES.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "06WD_Fcap_Besc.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "07WD_Candidate_MP_evaluation_reference_OM_ES_corrected.Rmd"), output_dir = pdf_dir)
 rmarkdown::render(file.path(wd_dir, "08WD_MP_robustness_and_implementation_ES.Rmd"), output_dir = pdf_dir)
 