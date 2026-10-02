# ============================================================
# 07_build_SR_reference.R
# Build bootstrap-specific stock-recruitment reference objects
#
# Stock / analysis:
#   European anchovy, ane.27.9aS (Gulf of Cádiz)
#
# Role in the workflow:
#   This script follows the bootstrap-specific Beverton-Holt fits
#   generated in 06_fit_SR_bootstrap.R. It links each conditioned OM
#   to its corresponding positive BH parameters and constructs the
#   FLSRsim stock-recruitment object used by subsequent OM expansion.
#
# Purpose:
#   - Link each conditioned OM to its corresponding BH parameters.
#   - Build the BH parameter array for the 100 conditioned OMs.
#   - Construct the FLSRsim stock-recruitment object.
#   - Preserve OM-specific sigmaR in the SR reference map.
#   - Diagnose variability among the 100 BH relationships.
#
# Inputs:
#   data/Rdata/biols_conditioned_100iter.rds
#   data/Rdata/iteration_map.csv
#   data/Rdata/bh_summary_positive.csv
#   data/Rdata/SR_historical_bootstrap.rds
#
# Outputs — computational objects:
#   data/Rdata/SRs_reference_conditioned_100iter.rds
#   data/Rdata/SR_ANE_reference_conditioned_100iter.rds
#   data/Rdata/SR_reference_iteration_map.csv
#
# Working Document outputs:
#   outputs/recruitment/figures/BH_reference_100OMs.png
#   outputs/recruitment/figures/BH_reference_100OMs.pdf
#   outputs/recruitment/tables/BH_reference_envelope.csv
#
#   The figure represents the 100 OM-specific Beverton-Holt
#   relationships and their median relationship over the historical
#   SSB range. The envelope table stores the 5th percentile, median
#   and 95th percentile predicted recruitment across OMs for the
#   common SSB grid..
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================

rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================
rds_boot_dir <- file.path(here(),"data","Rdata")
fig_dir <- file.path(here(),"outputs","recruitment","figures")
tab_dir <- file.path(here(),"outputs","recruitment","tables")

dir.create(fig_dir,recursive = TRUE,showWarnings = FALSE)
dir.create(tab_dir,recursive = TRUE,showWarnings = FALSE)

first.yr <- 1989
last.obs.yr <- 2024

ss.rec <- 3
ss.ssb <- 2

n_om <- 100

first.yr <- 1989
last.obs.yr <- 2024

ss.rec <- 3
ss.ssb <- 2

n_om <- 100

# ============================================================
# 2. LOAD INPUT OBJECTS
# ============================================================
# Load the conditioned biological OMs, the OM-bootstrap mapping and
# the positive BH parameter estimates generated in the previous step.

biols <- readRDS(file.path(rds_boot_dir,"biols_conditioned_100iter.rds"))
iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map.csv"))
bh_summary_pos <- read.csv(file.path(rds_boot_dir,"bh_summary_positive.csv"))

# ============================================================
# 3. INSPECT INPUT STRUCTURE
# ============================================================
# Inspect the classes, dimensions and identifiers of the input objects
# before constructing the stock-recruitment reference.

class(biols)
names(biols)
class(biols$ANE)
dim(biols$ANE@n)

dim(iteration_map)
names(iteration_map)

dim(bh_summary_pos)
names(bh_summary_pos)

head(iteration_map)
head(bh_summary_pos)

# ============================================================
# 4. ALIGN BH PARAMETERS WITH CONDITIONED OMs
# ============================================================
# Match BH estimates to conditioned OMs using both iteration and
# bootstrap identifiers. The explicit ordering by iteration is required
# because the parameter array below is filled by OM iteration.

sr_parameter_map <- iteration_map %>%
  dplyr::left_join(
    bh_summary_pos %>% dplyr::select(iter,bootstrap,convergence,a,b,sigmaR),
    by = c("iter","bootstrap")) %>%
  dplyr::arrange(iter)

# ============================================================
# 5. VALIDATE OM-PARAMETER MAPPING
# ============================================================
# Confirm that all 100 OMs have one complete, converged BH parameter
# set before the parameters are copied into the FLSRsim array.

stopifnot(nrow(sr_parameter_map) == n_om)
stopifnot(dplyr::n_distinct(sr_parameter_map$iter) == n_om)
stopifnot(identical(sr_parameter_map$iter,seq_len(n_om)))
stopifnot(!anyNA(sr_parameter_map$a))
stopifnot(!anyNA(sr_parameter_map$b))
stopifnot(!anyNA(sr_parameter_map$sigmaR))
stopifnot(all(sr_parameter_map$convergence == 0))

# ============================================================
# 6. BUILD STOCK-RECRUITMENT PARAMETER ARRAY
# ============================================================
# Construct a two-parameter array containing OM-specific BH a and b.
# Each parameter is constant across years and seasons within an OM.

yrs <- dimnames(biols$ANE@n)$year
seasons <- dimnames(biols$ANE@n)$season
iters <- dimnames(biols$ANE@n)$iter

length(yrs)
length(seasons)
length(iters)

sr_params <- array(
  NA_real_,
  dim = c(param = 2,year = length(yrs),season = length(seasons),iter = length(iters)),
  dimnames = list(param = c("a","b"),year = yrs,season = seasons,iter = iters))

for(i in seq_len(n_om)) {
  sr_params["a",,,i] <- sr_parameter_map$a[i]
  sr_params["b",,,i] <- sr_parameter_map$b[i]
}

# ============================================================
# 7. BUILD FLSRsim FOR 100 CONDITIONED OMs
# ============================================================
# Use the conditioned FLBiol dimensions as the template for the
# stock-recruitment object. Recruitment is assigned to Q3.
# Recruitment process-error multipliers are not introduced here.

sr_template <- biols$ANE@n[1,,,,,]

dim(sr_template)

sr_rec <- sr_template
sr_ssb <- sr_template
sr_unc <- sr_template
sr_prop <- sr_template

sr_rec[] <- NA_real_
sr_ssb[] <- NA_real_

sr_unc[] <- 1

sr_prop[] <- 0
sr_prop[,,,ac(ss.rec),,] <- 1

SR_ANE <- FLSRsim(
  name = "ANE",
  model = "bevholt",
  rec = sr_rec,
  ssb = sr_ssb,
  params = sr_params,
  uncertainty = sr_unc,
  proportion = sr_prop,
  covar = FLQuants()
)

SR_ANE@timelag["year",] <- 0
SR_ANE@timelag["season",] <- 2

# ============================================================
# 8. DIAGNOSTIC AND WORKING DOCUMENT OUTPUTS
# ============================================================
# Generate the diagnostic showing the 100 OM-specific Beverton-Holt
# relationships and their median relationship over the historical SSB
# range. The corresponding 5th percentile, median and 95th percentile
# envelope is saved as a table for Working Document reporting.

sr_hist <- readRDS(file.path(rds_boot_dir,"SR_historical_bootstrap.rds"))

SSB_grid <- seq(min(sr_hist$SSB,na.rm = TRUE),max(sr_hist$SSB,na.rm = TRUE),length.out = 300)

bh_curves_all <- tidyr::crossing(iter = seq_len(n_om),SSB = SSB_grid) %>%
  dplyr::left_join(sr_parameter_map %>% dplyr::select(iter,bootstrap,a,b),by = "iter") %>%
  dplyr::mutate(R_pred = a * SSB / (b + SSB))

bh_envelope <- bh_curves_all %>%
  dplyr::group_by(SSB) %>%
  dplyr::summarise(
    R_p05 = quantile(R_pred,0.05),
    R_median = median(R_pred),
    R_p95 = quantile(R_pred,0.95),
    .groups = "drop")

p_bh <- ggplot() +
  geom_line(
    data = bh_curves_all,
    aes(x = SSB,y = R_pred,group = iter,colour = factor(iter)),
    alpha = 0.55,
    linewidth = 0.45) +
  geom_line(data = bh_envelope,aes(x = SSB,y = R_median),colour = "black",linewidth = 1.2) +
  labs(x = "Spawning stock biomass (t)",y = "Recruitment",
    title = "Beverton-Holt relationships across bootstrap-conditioned OMs",
    subtitle = "100 bootstrap-specific relationships; black line = median") +
  guides(colour = "none") +
  theme_bw()

ggsave(file.path(fig_dir,"BH_reference_100OMs.png"),p_bh,width = 10,height = 7,dpi = 300)
ggsave(file.path(fig_dir,"BH_reference_100OMs.pdf"),p_bh,width = 10,height = 7)

write.csv(bh_envelope,file.path(tab_dir,"BH_reference_envelope.csv"),row.names = FALSE)
# ============================================================
# 9. BUILD SR REFERENCE MAP
# ============================================================
# Retain the explicit OM-bootstrap correspondence together with the
# BH parameters and sigmaR required by subsequent recruitment steps.

sr_reference_map <- sr_parameter_map %>%
  dplyr::select(iter,bootstrap,a,b,sigmaR)

# ============================================================
# 10. SAVE COMPUTATIONAL OUTPUTS
# ============================================================
# Store the FLSRsim reference objects and the validated parameter map
# required by the subsequent expansion from 100 OMs to 1000 trajectories.

SRs <- list(ANE = SR_ANE)

saveRDS(SRs,file.path(rds_boot_dir,"SRs_reference_conditioned_100iter.rds"))
saveRDS(SR_ANE,file.path(rds_boot_dir,"SR_ANE_reference_conditioned_100iter.rds"))
write.csv(sr_reference_map,file.path(rds_boot_dir,"SR_reference_iteration_map.csv"),row.names = FALSE)


